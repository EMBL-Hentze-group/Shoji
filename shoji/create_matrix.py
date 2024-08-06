import tempfile
from itertools import chain
from pathlib import Path
from shutil import rmtree
from typing import Callable, Dict, List, Optional, Set, Tuple

import pyarrow as pa
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
        cores: int,
        prefix: Optional[str] = None,
        suffix: Optional[str] = ".parquet",
        tmp_dir: Optional[str] = None,
    ) -> None:
        self.in_dir = in_dir
        self.annotation = annotation
        self.out = out
        self.prefix = prefix
        self.suffix = suffix
        self.cores: int = set_cores(cores)
        pa.set_cpu_count(self.cores)
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

    def __enter__(self) -> "CreateMatrix":
        self._glob_count_files()
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

    def _prepare_dataset(self) -> None:
        """_prepare_dataset Helper function
        Generate partitioned parquet dataset from the count files

        Raises:
            RuntimeError: If the number of unique sample_names and number of files in self._count_files does not match
        Returns:
            None
        """
        count_ds = ds.dataset(self._count_files, format="parquet").sort_by(
            [("chrom", "ascending"), ("start", "ascending")]
        )
        # sample names
        self._samples: List[str] = sorted(
            count_ds.to_table()["sample"].unique().to_pylist()
        )
        if len(self._samples) != len(self._count_files):
            raise RuntimeError(
                f"Mismatch in number of samples! found {len(self._count_files)} in {self.in_dir} but found {len(self._samples)} sample names in merged count files"
            )
        # windowed
        # TODO find a better way to do this!
        annotation = (count_ds.to_table()["annotation"][:100]).to_pylist()
        alen = list(set([len(a) for a in annotation]))
        if len(alen) != 1:
            raise RuntimeError(
                "Count files have different annotation schema! Check input files"
            )
        if alen[0] == 7:
            self._is_windowed = True
            logger.info("Count files are windowed")
        partition_schema = ds.partitioning(
            schema=pa.schema([("chrom", pa.string()), ("strand", pa.string())]),
            flavor="hive",
        )
        self._partitioned_ds.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Writing partitioned parquet file {str(self._partitioned_ds)}")
        ds.write_dataset(
            count_ds,
            base_dir=str(self._partitioned_ds),
            format="parquet",
            partitioning=partition_schema,
        )

    def create_matrices(self):
        partitioned_ds = ds.dataset(
            str(self._partitioned_ds), format="parquet", partitioning="hive"
        )
        fragments: Dict[Tuple[str, str], ds.ParquetFileFragment] = {}
        for fragment in partitioned_ds.get_fragments():
            partition_dict = ds.get_partition_keys(fragment.partition_expression)
            fragments[(partition_dict["chrom"], partition_dict["strand"])] = fragment
        annotation_suffix: str = Path(self.annotation).suffix
        out_suffix: str = Path(self.out).suffix
        annotation_tmp, out_tmp, max_tmp = {}, {}, {}
        first: bool = True
        for chrom in sorted(fragments.keys()):
            ann_file: str = str(
                self._tmp
                / f"{next(tempfile._get_candidate_names())}{annotation_suffix}"
            )
            out_file: str = str(
                self._tmp / f"{next(tempfile._get_candidate_names())}{out_suffix}"
            )
            annotation_tmp[f"{chrom[0]}{chrom[1]}"] = ann_file
            out_tmp[f"{chrom[0]}{chrom[1]}"] = out_file
            create_matrices(
                fragment=fragments[chrom],
                chrom=chrom[0],
                strand=chrom[1],
                samples=self._samples,
                ann=ann_file,
                sums=out_file,
                first=first,
            )
            first = False
        # accumulate and write annotations
        logger.info(f"writing annotations to {self.annotation}")
        general_accumulator(annotation_tmp, self.annotation)
        # accumulate and write window counts
        logger.info(f"writing counts to {self.out}")
        general_accumulator(out_tmp, self.out)


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
):
    ann_header: List[str] = [
        "unique_id",
        "chromosome",
        "start",
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
    # sum handlser
    sum_writer: Callable = output_writer(out=sums, use_tabix=False)
    # counts table
    counts: pa.Table = fragment.to_table()
    # genes
    genes: pa.StringArray = counts["gene_id"].unique()
    logger.debug(f"{fragment}")
    logger.info(f"Found {len(genes):,} genes in {chrom} {strand}")
    with ann_writer(ann) as aw, sum_writer(sums) as sw:
        if first:
            aw.write("\t".join(ann_header) + "\n")
            sw.write("\t".join(["unique_id"] + samples) + "\n")
        for gene in genes:
            per_gene = counts.filter(pc.field("gene_id") == gene)
            row_map: Dict[Tuple[int, int], WindowCount] = {}
            for row in per_gene.to_pylist():
                uid = (row["start"], row["end"])
                try:
                    row_map[uid].add(row["sample"], row["pos_counts"])
                except KeyError:
                    row_map[uid] = WindowCount()
                    row_map[uid].add(row["sample"], row["pos_counts"])
                    row_map[uid].annotations = dict(row["annotation"])
            uniq_rows: List[Tuple[int, int]] = pick_rows(row_map=row_map, strand=strand)
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
                aw.write("\t".join(annotation) + "\n")
                window_sums: List[str] = [row_map[ur].annotations["uniq_id"]]
                for sample in samples:
                    if sample not in row_map[ur].window_sum:
                        window_sums.append("0")
                        continue
                    window_sums.append(str(row_map[ur].window_sum[sample]))
                sw.write("\t".join(window_sums) + "\n")


def pick_rows(
    row_map: Dict[Tuple[int, int], WindowCount], strand: str
) -> List[Tuple[int, int]]:
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
