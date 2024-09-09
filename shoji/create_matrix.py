import multiprocessing as mp
import tempfile
from itertools import chain
from pathlib import Path
from shutil import rmtree
from typing import Callable, Dict, List, Optional, Tuple, Set

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
import pyarrow.dataset as ds

# import pyarrow.parquet as pq
from loguru import logger

from .helpers import set_cores
from .output import general_accumulator, output_writer


class CreateMatrix:
    def __init__(
        self,
        in_dir: str,
        annotation: str,
        out: str,
        max_out: Optional[str] = None,
        cores: int = 1,
        prefix: Optional[str] = None,
        suffix: Optional[str] = ".parquet",
        tmp_dir: Optional[str] = None,
    ) -> None:
        self.in_dir = in_dir
        self.annotation = annotation
        self.out = out
        self.max_out = max_out
        self.prefix = prefix
        self.suffix = suffix
        if tmp_dir is None:
            self._tmp = Path(self.out).parent / next(tempfile._get_candidate_names())
        else:
            if not Path(tmp_dir).exists():
                raise FileNotFoundError(f"Directory {tmp_dir} does not exists!")
            self._tmp: Path = Path(tmp_dir) / next(tempfile._get_candidate_names())
        # count files
        self._count_files: List[Path] = []
        # sample names
        self._samples: List[str] = []
        # windowed
        self._is_windowed: bool = False
        # temp. partitioned parquet dataset
        self._partitioned_ds = self._tmp / next(tempfile._get_candidate_names())
        if Path(self.annotation).exists():
            logger.warning(f"Re-writing file {self.annotation}")
        if Path(self.out).exists():
            logger.warning(f"Re-writing file {self.out}")
        if self.max_out is None:
            logger.info("No output for max count matrix.")
        elif Path(self.max_out).exists():
            logger.warning(f"Re-writing file {self.max_out}")
        # set cores
        self.cores: int = set_cores(cores)
        # pa.set_cpu_count(self.cores)
        # divide up cores amont arrow and mp
        pa_cores: int = max(1, int(self.cores / 2))
        pa.set_cpu_count(pa_cores)
        self._rest_cores = max(1, self.cores - pa_cores)
        logger.debug(f"Using {pa_cores} for arrow and {self._rest_cores} for mp")

    def __enter__(self) -> "CreateMatrix":
        self._glob_count_files()
        self._sanity_check()
        self._prepare_dataset()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        rmtree(self._tmp)

    def _glob_count_files(self) -> None:
        """_glob_count_files Helper function
        Generate a list of count files

        Raises:
            FileNotFoundError: If no files match the constructed wildcard pattern in `in_dir`.
        Returns:
            None
        """
        wildcard: str = ""
        if self.prefix is None:
            wildcard = "*" + self.suffix
        elif self.suffix is None:
            wildcard = self.prefix + "*"
        else:
            wildcard = self.prefix + "*" + self.suffix
        self._count_files = sorted(Path(self.in_dir).glob(wildcard))
        if len(self._count_files) == 0:
            raise FileNotFoundError(
                f"Cannot find count files in {self.in_dir} with wildcard {wildcard}"
            )
        logger.info(f"Found {len(self._count_files)} count files in {self.in_dir}")

    def _sanity_check(self) -> None:
        """_sanity_check Helper function
        Checks whether all input files have unique sample names
        checks whether all input files have the same annotation length

        Raises:
            RuntimeError: If the number of unique sample_names and number of files in self._count_files does not match
            RuntimeError: If the number of unique annotation elements are not the same across all samples
        Returns:
            None
        """
        annotation_size: Set[int] = set()
        for c in self._count_files:
            c_table = pq.read_table(c, columns=["annotation", "sample"])
            sample_size = int(np.sqrt(c_table.num_rows))
            # take 50 random rows, or the whole table
            ns = np.random.randint(0, sample_size, min(sample_size, 50))
            logger.info(f"file: {str(c)}, rows: {c_table.num_rows:,}")
            # add sample names
            self._samples.extend(c_table.take(ns)["sample"].unique().to_pylist())
            # annotation
            annotation_size.update(
                {len(a) for a in c_table.take(ns)["annotation"].to_pylist()}
            )
        # check if all files have same annotation schema
        if len(annotation_size) != 1:
            raise RuntimeError(
                "Count files have different annotation schema! Check input files"
            )
        ann_len: int = list(annotation_size)[0]
        # windowed tables have 7 annotation elements, 6 for non windowed
        if ann_len == 7:
            self._is_windowed = True
            logger.info("Count files are windowed")
        self._samples = sorted(self._samples)
        # number of sample names MUST match number of input files!
        if len(self._samples) != len(self._count_files):
            raise RuntimeError(
                f"Mismatch in number of samples! found {len(self._count_files)} in {self.in_dir} but found {len(self._samples)} sample names in merged count files"
            )

    def _prepare_dataset(self) -> None:
        """_prepare_dataset Helper function
        Generate partitioned parquet dataset from the count files
        Returns:
            None
        """
        partition_schema = ds.partitioning(
            schema=pa.schema([("chrom", pa.string()), ("strand", pa.string())]),
            flavor="hive",
        )
        count_ds = ds.dataset(self._count_files, format="parquet")
        self._partitioned_ds.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Writing partitioned parquet file {str(self._partitioned_ds)}")
        ds.write_dataset(
            count_ds,
            base_dir=str(self._partitioned_ds),
            format="parquet",
            partitioning=partition_schema,
        )

    def create_matrices(self, allow_duplicates: bool = False):
        partitioned_ds = ds.dataset(
            str(self._partitioned_ds), format="parquet", partitioning="hive"
        )
        fragments: Dict[Tuple[str, str], ds.ParquetFileFragment] = {}
        for fragment in partitioned_ds.get_fragments():
            partition_dict = ds.get_partition_keys(fragment.partition_expression)
            fragments[(partition_dict["chrom"], partition_dict["strand"])] = fragment
        annotation_suffix: str = Path(self.annotation).suffix
        out_suffix: str = Path(self.out).suffix
        if self.max_out is None:
            max_suffix = None
        else:
            max_suffix = Path(self.out).suffix
        annotation_tmp, out_tmp, max_tmp, first_map = {}, {}, {}, {}
        first: bool = True
        for chrom in sorted(fragments.keys()):
            # iterate over chromosomes and generate output files
            annotation_tmp[f"{chrom[0]}{chrom[1]}"] = str(
                self._tmp
                / f"{next(tempfile._get_candidate_names())}{annotation_suffix}"
            )
            out_tmp[f"{chrom[0]}{chrom[1]}"] = str(
                self._tmp / f"{next(tempfile._get_candidate_names())}{out_suffix}"
            )
            if self.max_out is None:
                max_tmp[f"{chrom[0]}{chrom[1]}"] = None
            else:
                max_tmp[f"{chrom[0]}{chrom[1]}"] = str(
                    self._tmp / f"{next(tempfile._get_candidate_names())}{max_suffix}"
                )
            first_map[chrom] = first
            first = False
        with mp.Pool(self._rest_cores) as pool:
            for chrom in sorted(fragments.keys()):
                ann_file: str = annotation_tmp[f"{chrom[0]}{chrom[1]}"]
                out_file: str = out_tmp[f"{chrom[0]}{chrom[1]}"]
                max_file: str = max_tmp[f"{chrom[0]}{chrom[1]}"]
                first: bool = first_map[chrom]
                pool.apply_async(
                    create_matrices,
                    args=(
                        fragments[chrom],
                        chrom[0],
                        chrom[1],
                        self._samples,
                        ann_file,
                        out_file,
                        max_file,
                        first,
                        self._is_windowed,
                        allow_duplicates,
                    ),
                )
            pool.close()
            pool.join()
        # accumulate and write annotations
        logger.info(f"writing annotations to {self.annotation}")
        general_accumulator(annotation_tmp, self.annotation)
        # accumulate and write window counts
        logger.info(f"writing counts to {self.out}")
        general_accumulator(out_tmp, self.out)
        # accumulate and write max counts
        if self.max_out is not None:
            logger.info(f"writing max counts to {self.max_out}")
            general_accumulator(max_tmp, self.max_out)


class WindowCount:
    def __init__(self) -> None:
        self._sample_counts: Dict[str, List[Tuple[int, int]]] = {}
        self._annotations: Dict[str, str] = {}

    def __hash__(self) -> int:
        keys = sorted(self._sample_counts.keys())
        values = sorted(chain(*self._sample_counts.values()))
        return hash(tuple(keys + values))

    def add(self, sample_name, pos_count):
        if sample_name in self._sample_counts:
            raise RuntimeError(f"{sample_name}  already exists!")
        self._sample_counts[sample_name] = pos_count

    @property
    def samples(self):
        return sorted(self._sample_counts.keys())

    @property
    def annotations(self) -> Dict[str, str]:
        return self._annotations

    @annotations.setter
    def annotations(self, annotations) -> None:
        self._annotations = annotations

    @property
    def window_sum(self) -> Dict[str, int]:
        wsums: Dict[str, int] = {}
        for sample, dat in self._sample_counts.items():
            wsums[sample] = sum([d[1] for d in dat])
        return wsums

    @property
    def window_max(self) -> Dict[str, int]:
        wmax: Dict[str, int] = {}
        for sample, dat in self._sample_counts.items():
            wmax[sample] = max([d[1] for d in dat])
        return wmax


def create_matrices(
    fragment: ds.ParquetFileFragment,
    chrom: str,
    strand: str,
    samples: List[str],
    ann: str,
    sums: str,
    maxs: Optional[str] = None,
    first: bool = False,
    is_windowed: bool = False,
    allow_duplicates: bool = False,
):
    ann_header: List[str] = [
        "unique_id",
        "chromosome",
        "begin",
        "end",
        "strand",
        "gene_id",
        "gene_name",
        "gene_type",
        "gene_region",
        "Nr_of_region",
        "Total_nr_of_region",
    ]
    if is_windowed:
        ann_header.append("window_number")
    # annotation handler
    ann_writer: Callable = output_writer(out=ann, use_tabix=False)
    # sum handler
    sum_writer: Callable = output_writer(out=sums, use_tabix=False)
    # counts table
    counts: pa.Table = fragment.to_table().sort_by("begin")
    # genes
    genes: pa.StringArray = counts["gene_id"].unique()
    logger.debug(f"{fragment}")
    logger.info(f"Found {len(genes):,} genes in {chrom} {strand}")
    # max counts
    max_counts: List[List[str]] = []
    if allow_duplicates:
        row_fn: Callable = all_rows
    else:
        row_fn: Callable = pick_rows
    with ann_writer(ann) as aw, sum_writer(sums) as sw:
        if first:
            aw.write("\t".join(ann_header) + "\n")
            sw.write("\t".join(["unique_id"] + samples) + "\n")
        for gene in genes:
            row_map: Dict[Tuple[int, int], WindowCount] = {}
            for row in counts.filter(pc.field("gene_id") == gene).to_pylist():
                uid = (row["begin"], row["end"])
                try:
                    row_map[uid].add(row["sample"], row["pos_counts"])
                except KeyError:
                    row_map[uid] = WindowCount()
                    row_map[uid].add(row["sample"], row["pos_counts"])
                    row_map[uid].annotations = dict(row["annotation"])
            uniq_rows: List[Tuple[int, int]] = row_fn(row_map=row_map, strand=strand)
            for ur in uniq_rows:
                annotation: List[str] = [
                    row_map[ur].annotations["uniq_id"],
                    chrom,
                    str(ur[0]),
                    str(ur[1]),
                    strand,
                    str(gene),
                    row_map[ur].annotations["gene_name"],
                    row_map[ur].annotations["gene_type"],
                    row_map[ur].annotations["feature"],
                    row_map[ur].annotations["nr_of_region"],
                    row_map[ur].annotations["total_region"],
                ]
                if is_windowed:
                    annotation.append(row_map[ur].annotations["window_number"])
                aw.write("\t".join(annotation) + "\n")
                window_sums: List[str] = [row_map[ur].annotations["uniq_id"]]
                window_maxs: List[str] = [row_map[ur].annotations["uniq_id"]]
                for sample in samples:
                    if sample not in row_map[ur].window_sum:
                        window_sums.append("0")
                        window_maxs.append("0")
                        continue
                    window_sums.append(str(row_map[ur].window_sum[sample]))
                    window_maxs.append(str(row_map[ur].window_max[sample]))
                sw.write("\t".join(window_sums) + "\n")
                max_counts.append(window_maxs)
    if maxs is not None:
        max_writer: Callable = output_writer(out=maxs, use_tabix=False)
        with max_writer(maxs) as mw:
            if first:
                mw.write("\t".join(["unique_id"] + samples) + "\n")
            for m in max_counts:
                mw.write("\t".join(m) + "\n")


def all_rows(
    row_map: Dict[Tuple[int, int], WindowCount], strand: str
) -> List[Tuple[int, int]]:
    """all_rows Helper function
    Return all rows. Proxy function for pick_rows
    Args:
        row_map (Dict[Tuple[int, int], WindowCount]): Intervals with count data
        strand (str): gene strand

    Returns:
        List[Tuple[int, int]]: list of non overlapping intervals
    """
    return sorted(row_map.keys())


def pick_rows(
    row_map: Dict[Tuple[int, int], WindowCount], strand: str
) -> List[Tuple[int, int]]:
    """pick_rows Helper function
    From overlapping intervals with identical crosslink data, pick one
    Args:
        row_map (Dict[Tuple[int, int], WindowCount]): Intervals with count data
        strand (str): gene strand

    Returns:
        List[Tuple[int, int]]: list of non overlapping intervals
    """
    rows: List[Tuple[int, int]] = []
    if strand == "-":
        _fn: Callable = max
    else:
        _fn: Callable = min
    hash_map: Dict[int, List[Tuple[int, int]]] = {}
    for row, dat in row_map.items():
        try:
            hash_map[hash(dat)].append(row)
        except KeyError:
            hash_map[hash(dat)] = [row]
    for pos in hash_map.values():
        if len(pos) == 1:
            rows.append(pos[0])
        else:
            rows.extend(_row_picker(pos, _fn))
    return rows


def _row_picker(
    rows: List[Tuple[int, int]], picker_fn: Callable
) -> List[Tuple[int, int]]:
    """_row_picker Helper function to pick rows from a list of tuples
    From a list intervals, for the overlapping set of intervals pick
    the most 5' interval depending on the strand.
    Args:
        rows (List[Tuple[int, int]]): list of intervals
        picker_fn (Callable): min or max function

    Returns:
        List[Tuple[int, int]]: list of non overlapping intervals
    """
    uniq_rows: List[Tuple[int, int]] = []
    irows = iter(sorted(rows))
    try:
        start, end = next(irows)
        for nstart, nend in irows:
            if nstart > end:
                uniq_rows.append((start, end))
                start, end = nstart, nend
            else:
                start = picker_fn(start, nstart)
                end = picker_fn(end, nend)
        uniq_rows.append((start, end))
    except StopIteration:
        return
    return uniq_rows
