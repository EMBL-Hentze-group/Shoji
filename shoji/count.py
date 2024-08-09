import multiprocessing as mp
import tempfile
from bisect import bisect_left, bisect_right
from collections import defaultdict
from functools import partial
from pathlib import Path
from shutil import copyfile, rmtree
from typing import Any, Callable, Dict, Generator, List, Optional, Set, Union

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pysam
from loguru import logger
from pyarrow import dataset as pds
from pyarrow.dataset import ParquetFileFragment

from .helpers import BedFeature, check_tabix, set_cores, Crosslinks
from .pyarrow_reader import PartionedParquetReader


class Count:
    def __init__(
        self,
        annotation: str,
        bed: str,
        out: str,
        cores: int,
        sample_name: Optional[str] = None,
        tmp_dir: Optional[str] = None,
    ) -> None:
        """__init__ _summary_

        Args:
            annotation: str, Path to the annotation file in BED format.
            bed: str, Path to the bed file in BED format.
            out: str, Path to the output directory.
            cores: int, Number of cores to use for parallelization
            sample_name: str, sample name to use. If none provided, sample name will be inferred from bed file name.
            tmp_dir: Optional, str, Path to a temporary directory for intermediate files. If not provided, a temporary directory will be automatically created within the output directory

        Raises:
            FileNotFoundError: If the specified tmp_dir does not exist.
        """

        self.annotation = annotation
        self.bed = bed
        self.out = Path(out)
        self.cores: int = set_cores(cores)
        pa.set_cpu_count(self.cores)
        if tmp_dir is None:
            self._tmp = Path(self.out).parent / next(tempfile._get_candidate_names())
        else:
            if not Path(tmp_dir).exists():
                raise FileNotFoundError(f"Directory {tmp_dir} does not exists!")
            self._tmp: Path = Path(tmp_dir) / next(tempfile._get_candidate_names())
        self._tmp.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using {str(self._tmp)} as temp. directory")
        if sample_name is None:
            bpath: Path = Path(bed)
            self.sample_name = bpath.name.replace("".join(bpath.suffixes), "")
            logger.warning(
                f"No sample name provided. Using '{self.sample_name}' inferred from file name '{bed}' as sample name"
            )
        else:
            self.sample_name = sample_name
            logger.info(f"Using '{self.sample_name}' as sample name for '{bed}'")
        # annotation handler and reader
        self._annh: Union[pysam.TabixFile, PartionedParquetReader]
        # boolean
        self._ann_tabix: bool = False
        # annotation reader
        self._ann_reader: Callable
        # annotation paritioned fragments
        self._ann_fragments: Dict[str, ParquetFileFragment] = {}
        # count/bed handler and reader
        self._bedh: Union[pysam.TabixFile, PartionedParquetReader]
        # boolean
        self._bed_tabix: bool = False
        # bed/count reader
        self._bed_reader: Callable
        # bed/count paritioned fragments
        self._bed_fragments: Dict[str, ParquetFileFragment] = {}

    def __enter__(self) -> "Count":
        self._check_out_suffix()
        self._get_annotation_handler()
        self._annh.__enter__()
        self._get_count_handler()
        self._bedh.__enter__()
        return self

    def __exit__(self, except_type, except_val, except_traceback):
        self._annh.__exit__(except_type, except_val, except_traceback)
        self._bedh.__exit__(except_type, except_val, except_traceback)
        rmtree(self._tmp)

    def _check_out_suffix(self):
        """_check_out_suffix Helper function
        Check whether the output suffix is .parquet, otherwise warn
        TODO: perhaps change suffix ?
        """
        if self.out.suffix != ".parquet":
            logger.warning(
                "Output sliding windows count will be written only in Apache Parquet format. File name suffix inconsistency!"
            )
        if self.out.exists():
            logger.warning(f"Rewriting file {str(self.out)}")

    def _get_annotation_handler(self) -> None:
        """_get_annotation_handler Helper function
        Appropriate file handle and associated read function
        depending on file type for annotation data
        """
        if check_tabix(self.annotation):
            logger.debug(f"{self.annotation} is tabix indexed")
            self._annh = pysam.TabixFile(self.annotation)
            self._ann_reader = tabix_reader
            self._ann_tabix = True
        else:
            logger.debug(f"Using pyarrow for {self.annotation}")
            self._annh = PartionedParquetReader(
                file_name=self.annotation,
                fformat="bed6",
                temp_dir=str(self._tmp),
                cores=self.cores,
            )
            self._ann_fragments = self._annh.get_partitioned_fragments()
            self._ann_reader = fragment_reader

    def _get_count_handler(self) -> None:
        """_get_count_handler Helper function
        Appropriate file handle and associated read function
        depending on file type for count data
        """
        if check_tabix(self.bed):
            logger.debug(f"{self.bed} is tabix indexed")
            self._bedh = pysam.TabixFile(self.bed)
            self._bed_reader = tabix_reader
            self._bed_tabix = True
        else:
            logger.debug(f"Using pyarrow for {self.bed}")
            self._bedh = PartionedParquetReader(
                file_name=self.bed,
                fformat="bed6",
                temp_dir=str(self._tmp),
                cores=self.cores,
            )
            self._bed_fragments = self._bedh.get_partitioned_fragments()
            self._bed_reader = fragment_reader

    def _get_ann_fragment(self, chrom: str) -> Union[str, ParquetFileFragment]:
        """_get_ann_fragment Helper function
        If annotation file is tabix indexed, return the original file name.
        Else, return parquet file fragment for `chrom` from self._ann_fragments
        Args:
            chrom: str, chromosome name

        Returns:
            str or ParquetFileFragment
        """
        if self._ann_tabix:
            logger.debug(f"{self.annotation} is tabix indexed")
            return self.annotation
        return self._ann_fragments[chrom]

    def _get_bed_fragment(self, chrom: str) -> Union[str, ParquetFileFragment]:
        """_get_bed_fragment Helper function
        If count file is tabix indexed, return the original file name.
        Else, return parquet file fragment for `chrom` from self._bed_fragments
        Args:
            chrom: str, chromosome name

        Returns:
            str or ParquetFileFragment
        """
        if self._bed_tabix:
            logger.debug(f"{self.bed} is tabix indexed")
            return self.bed
        return self._bed_fragments[chrom]

    def count(self) -> None:
        """count
         Count crosslink positions across common chromosomes found between annotation and count files.

        Returns:
            None

        Raises:
            RuntimeError: If there are no common chromosomes between Annotation and Bed files, suggesting input files may be incorrect or misaligned.
        """
        common_chroms: List[str] = sorted(
            set(self._annh.contigs) & set(self._bedh.contigs)
        )
        if len(common_chroms) == 0:
            raise RuntimeError(
                f"There are no common chromosomes between {self.annotation} and {self.bed}! Check your input files"
            )
        logger.info(
            f"{self.annotation} and {self.bed} have {len(common_chroms)} common chromosomes"
        )
        logger.debug(f"Common chromosomes: {', '.join(common_chroms)}")
        # schema for output file
        schema: pa.Schema = pa.schema(
            {
                "chrom": pa.string(),
                "begin": pa.uint32(),
                "end": pa.uint32(),
                "gene_id": pa.string(),
                "annotation": pa.map_(pa.string(), pa.string()),
                "strand": pa.string(),
                "sample": pa.string(),
                "pos_counts": pa.map_(pa.uint32(), pa.uint32()),
            },
            metadata={
                "chrom": "chromosome name",
                "begin": "window start position",
                "end": "window end position",
                "gene_id": "gene id",
                "annotation": "key: feature names, values: feature attributes",
                "strand": "strand info",
                "sample": "sample name",
                "pos_counts": "keys: positions, values: crosslink counts",
            },
        )
        tmp_dict: Dict[str, Path] = {}  # tmp file dictionary
        with mp.Pool(processes=self.cores) as pool:
            for chrom in common_chroms:
                ann_fn: Callable = partial(
                    self._ann_reader, fragment=self._get_ann_fragment(chrom)
                )
                bed_fn: Callable = partial(
                    self._bed_reader, fragment=self._get_bed_fragment(chrom)
                )
                tmp_file: Path = (
                    self._tmp
                    / f"{chrom}_{next(tempfile._get_candidate_names())}.parquet"
                )
                tmp_dict[chrom] = tmp_file
                pool.apply_async(
                    count_crosslinks,
                    args=(ann_fn, bed_fn, chrom, schema, tmp_file, self.sample_name),
                )
            pool.close()
            pool.join()
        self._combine_outputs(tmp_dict)

    def _combine_outputs(self, tmp_dict: Dict[str, Path]) -> None:
        """
        Combines multiple count files into a single Parquet file with cross-link window counts.

        Args:
            tmp_dict (Dict[str, Path]): A dictionary mapping identifiers to count files' paths.

        Raises:
            RuntimeError:
                - If no valid count files are found in input directories,
                it raises an error indicating that crosslink window counts cannot be created
                due to a lack of data points in both `self.annotation` and `self.bed`.
            RuntimeError:
                - When more than one intermediate Parquet file is generated during the process,
                which indicates multiple inputs were not handled correctly, raising an error explaining the issue.

        Returns:
            None

        """
        data_files: List[Path] = []
        for chrom in sorted(tmp_dict.keys()):
            # so that chromosomes appear in sorted order
            p = tmp_dict[chrom]
            if p.exists():
                data_files.append(p)
        if len(data_files) == 0:
            raise RuntimeError(
                f"Cannot create crosslink window count for windows in {self.annotation} based on counts in {self.bed}. Check the input files!"
            )
        merged = pds.dataset(data_files, format="parquet")
        tmp_out_name = (
            self._tmp / f"{next(tempfile._get_candidate_names())}"
        )  # temp. dir for storing parquet file
        pds.write_dataset(
            merged, tmp_out_name, format="parquet", min_rows_per_group=30000
        )
        pq_path: List[Path] = sorted(tmp_out_name.glob("*.parquet"))
        if len(pq_path) > 1:
            raise RuntimeError(
                "Expected only one intermediate output file, but found multiple intermediate outputs! Cannot handle multiple intermediate output files"
            )
        _ = copyfile(pq_path[0], self.out)


def tabix_reader(fragment: str, chrom: str) -> Generator[pysam.BedProxy, None, None]:
    """
    Creates a generator to iterate over a tabix-indexed file for a specific chromosome.

    Args:
        fragment: str, The path to the tabix-indexed file.
        chrom: str, The chromosome to iterate over.

    Yields:
        A tuple containing the chromosome, start position, end position, and
        annotation for each record in the tabix file for the specified chromosome.
    """
    with pysam.TabixFile(fragment, parser=pysam.asBed()) as th:
        for dat in th.fetch(chrom, multiple_iterators=True):
            yield dat


def fragment_reader(
    fragment: ParquetFileFragment, chrom: str, batch_size: int = 5000
) -> Generator[BedFeature, None, None]:
    """
    Reads BedFeatures from a ParquetFileFragment.

    Args:
        fragment: str, A ParquetFileFragment object.
        chrom: str, The chromosome to read.
        batch_size: int, The maximum row count for scanned record batches

    Yields:
        BedFeature objects.
    """
    for chunk in fragment.to_batches(batch_size=batch_size):
        for dat in chunk.to_pylist():
            yield BedFeature(
                contig=chrom,
                start=dat["chromStart"],
                end=dat["chromEnd"],
                name=dat["name"],
                score=dat["score"],
                strand=dat["strand"],
            )


def _crosslink_counter(bed_fn: Callable, chrom: str) -> Dict[str, Crosslinks]:
    """_crosslink_counter Helper function

    Count crosslink sites within a specified chromosome.

    Args:
        bed_fn (Callable): A callable that accepts a string representing
                            a chromosome and yields BED entry entries for this chromosome.
        chrom (str): The chromosome identifier to analyze crosslink events on.

    Returns:
        Dict[str, Crosslinks]: A dictionary mapping strand orientation ('+', '-') to a Crosslinks
                                object containing the count of crosslink events at each position and their corresponding positions.

    Raises:
        RuntimeError: If no crosslink sites are found for the specified chromosome in the input files,
                        indicating an issue with the input data or its format.

    Note:
        The function relies on a `Crosslinks` dataclass (not shown here) to store the count of crosslink events at each position
        and their corresponding positions as attributes.
    """
    crosslinks_counter: defaultdict = defaultdict(dict)
    crosslinks: Dict[str, Crosslinks] = {}
    for bed in bed_fn(chrom=chrom):
        for pos in range(bed.end, bed.start, -1):
            try:
                crosslinks_counter[bed.strand][pos] += 1
            except KeyError:
                crosslinks_counter[bed.strand][pos] = 1
    if len(crosslinks_counter) == 0:
        raise RuntimeError(
            f"Cannot find crosslink sites for {chrom}. Check input files!"
        )
    for strand, dat in crosslinks_counter.items():
        counts: np.ndarray = np.array(list(dat.items()), dtype=np.uint32)
        counts = counts[counts[:, 0].argsort()]
        crosslinks[strand] = Crosslinks(counts=counts, pos=counts[:, 0].tolist())
        logger.info(
            f"{chrom}, strand: {strand} positions with crosslinks: {len(dat):,}"
        )
    return crosslinks


def count_crosslinks(
    annotation_fn: Callable,
    bed_fn: Callable,
    chrom: str,
    schema: pa.Schema,
    output: Path,
    sample_name: str,
):
    """count_crosslinks count crosslinks for a chromosome
     Count crosslinks within specified regions and annotate them with corresponding gene information.

     Args:
        annotation_fn (Callable): A callable that returns an iterator over annotated genomic regions.
        bed_fn (Callable): A callable that returns an iterator over BED-formatted crosslink positions.
        chrom (str): The chromosome for which to process the data.
        schema (pa.Schema): An Arrow Schema defining the structure of the output table.
        output (Path): The file path or name where the results will be written.
        sample_name (str): The current sample name. Used as a column in output file
        gene_id_index (int, optional): Attribute index within name column at which to find the gene ID. Defaults to 0.
        uniq_id_index (int, optional): Attribute index within the name column at which to find the unique identifier. Defaults to 5.

    Returns:
        None

    Raises:
        RuntimeError: If there is a discrepancy between expected and actual column names in the output schema or
                        if crosslinks data cannot be found for specified chromosome and strand combinations.

    Note:
        This function expects that `annotation_fn` returns an iterator over dictionaries with keys:
        "chrom", "start", "end", "gene_id", "uniq_id", "strand" and any additional required metadata
    """
    expected_cols: Set[str] = {
        "chrom",
        "begin",
        "end",
        "gene_id",
        "annotation",
        "strand",
        "sample",
        "pos_counts",
    }
    if expected_cols != set(schema.names):
        raise RuntimeError(
            f"Output column name discrepancy! expected column names:{','.join(expected_cols)} found column names: {','.join(schema.names)}"
        )
    logger.debug(
        f"Process id {mp.current_process().pid}, chromosome: {chrom}, temp. file: {output}"
    )
    crosslinks: Dict[str, Crosslinks] = _crosslink_counter(bed_fn, chrom)
    # output data
    out_dict: Dict[str, List[Any]] = dict([(name, []) for name in schema.names])
    total, used = 0, 0
    for ann in annotation_fn(chrom=chrom):
        total += 1
        if ann.strand not in crosslinks:
            logger.warning(
                f"Skipping {chrom} window {ann.start}:{ann.end}({ann.strand}). Cannot find crosslink counts for {chrom} {ann.strand}!"
            )
            continue
        start_pos = bisect_left(crosslinks[ann.strand].pos, ann.start)
        if start_pos == len(crosslinks[ann.strand].pos):
            # ann.start is downstream of all crosslink positions
            # skip checking this window
            continue
        end_pos = bisect_right(crosslinks[ann.strand].pos, ann.end, lo=start_pos)
        if (end_pos == 0) or (start_pos == end_pos):
            # ann.end is upstream of all crosslink positions
            # or this window falls in between two crosslink sites, ie no crosslink sites in this window
            # skip checking this window
            continue
        indices: np.ndarray = np.where(
            (crosslinks[ann.strand].counts[:, 0] > ann.start)
            & (crosslinks[ann.strand].counts[:, 0] <= ann.end)
        )[0]
        if len(indices) == 0:
            continue
        name_dat: List[str] = ann.name.split("@")
        nr_region, total_nr_region = name_dat[4].split("/")
        if len(name_dat) < 6:
            logger.warning(
                f"{chrom} Window {ann.start}:{ann.end}({ann.strand}) does not have enough attributes seperated by '@' in name column: {ann.name}. Expected at least 6 attributes!"
            )
            continue
        out_dict["chrom"].append(chrom)
        out_dict["begin"].append(ann.start)
        out_dict["end"].append(ann.end)
        out_dict["gene_id"].append(name_dat[0])
        out_dict["strand"].append(ann.strand)
        out_dict["sample"].append(sample_name)
        out_dict["pos_counts"].append(
            list(map(tuple, crosslinks[ann.strand].counts[indices]))
        )
        annotation: Dict[str, str] = {
            "gene_name": name_dat[1],
            "gene_type": name_dat[2],
            "feature": name_dat[3],
            "nr_of_region": nr_region,
            "total_region": total_nr_region,
            "uniq_id": name_dat[5],
        }
        if len(name_dat) == 7:
            annotation["window_number"] = name_dat[6]
        out_dict["annotation"].append(annotation)
        # out_dict["positions"].append(
        #     (crosslinks[ann.strand].counts[indices, 0]).tolist()
        # )
        # out_dict["counts"].append((crosslinks[ann.strand].counts[indices, 1]).tolist())
        used += 1
    if len(out_dict) > 0:
        logger.info(
            f"{chrom}: total windows: {total:,} windows with crosslink sites: {used:,}"
        )
        out_table = pa.Table.from_pydict(out_dict, schema=schema)
        pq.write_table(out_table, output)
    else:
        logger.warning("f{chrom}: cannot compute crosslink counts within windows!")
