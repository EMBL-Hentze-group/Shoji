import multiprocessing as mp
import tempfile
from gzip import open as gzopen
from itertools import chain
from pathlib import Path
from shutil import rmtree
from typing import Any, Callable, Dict, List, Optional, Set, Tuple

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from loguru import logger
from pyarrow.csv import CSVWriter, WriteOptions

from .helpers import set_cores
from .schemas import count_schema, get_annotation_schema


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
        format: str = "csv",
        tmp_dir: Optional[str] = None,
    ) -> None:
        """CreateMatrix Class
        Args:
            in_dir (str): Input directory with count files
            annotation (str): Output file for the annotation matrix
            out (str): Output file for the sum matrix
            max_out (Optional[str], optional): Output file for the max count matrix. Defaults to None.
            cores (int, optional): Number of cores to use. Defaults to 1.
            prefix (Optional[str], optional): Prefix for the input files. Defaults to None.
            suffix (Optional[str], optional): Suffix for the input files. Defaults to ".parquet".
            format (str, optional): Format of the output files. Defaults to "csv", could be either "csv" or "parquet".
            tmp_dir (Optional[str], optional): Temporary directory for intermediate files. Defaults to None.
        """
        self.in_dir = in_dir
        self.annotation = annotation
        self.out = out
        self.max_out = max_out
        self.prefix = prefix
        self.suffix = suffix
        self.format = format
        if tmp_dir is None:
            self._tmp = Path(self.out).parent / next(
                tempfile._get_candidate_names()  # type: ignore
            )
        else:
            if not Path(tmp_dir).exists():
                raise FileNotFoundError(f"Directory {tmp_dir} does not exists!")
            self._tmp: Path = Path(tmp_dir) / next(
                tempfile._get_candidate_names()  # type: ignore
            )
        # count files
        self._count_files: List[Path] = []
        # sample names
        self._samples: List[str] = []
        # windowed
        self._is_windowed: bool = False
        # temp. partitioned parquet dataset
        self._partitioned_ds = self._tmp / next(tempfile._get_candidate_names())  # type: ignore
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
        # @TODO: is there a better way to do this?
        pa_cores: int = max(1, self.cores // 2)
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
            wildcard = f"*{self.suffix}"
        elif self.suffix is None:
            wildcard = f"{self.prefix}*"
        else:
            wildcard = f"{self.prefix}*{self.suffix}"
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
        window_col: Set[bool] = set()
        for c in self._count_files:
            c_table = pq.read_table(
                c, schema=count_schema, columns=["window_number", "sample"]
            )
            # add sample names
            self._samples.extend(c_table.column("sample").unique().to_pylist())
            # check if the window number column is null
            window_col.add(
                c_table.column("window_number").unique().is_null()[0].as_py()
            )
        # check if all files have same annotation schema
        if len(window_col) != 1:
            raise RuntimeError(
                "Count files appear to be a mix of windowed and non windowed data! Check input files"
            )
        if not window_col.pop():
            # if first element, and the only element is False, then all are windowed
            # since checking is_null()
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
        count_ds = ds.dataset(self._count_files, format="parquet", schema=count_schema)
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
        ann_files, out_files, max_files = [], [], []
        with mp.Pool(self._rest_cores) as pool:
            for chrom in sorted(fragments.keys()):
                tmp_name: str = next(tempfile._get_candidate_names())  # type: ignore
                ann_file: str = str(self._tmp / f"ann_{chrom[0]}_{tmp_name}.parquet")
                ann_files.append(ann_file)
                out_file: str = str(self._tmp / f"count_{chrom[0]}_{tmp_name}.parquet")
                out_files.append(out_file)
                max_file: Optional[str] = None
                if self.max_out is not None:
                    max_file = str(self._tmp / f"maxs_{chrom[0]}_{tmp_name}.parquet")
                    max_files.append(max_file)
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
                        self._is_windowed,
                        allow_duplicates,
                    ),
                )
            pool.close()
            pool.join()
        # write out files
        if self.format == "parquet":
            self._parquet_writer(ann_files, self.annotation)
            self._parquet_writer(out_files, self.out)
            if self.max_out is not None:
                self._parquet_writer(max_files, self.max_out)
        elif self.format == "csv":
            self._csv_writer(ann_files, self.annotation)
            self._csv_writer(out_files, self.out)
            if self.max_out is not None:
                self._csv_writer(max_files, self.max_out)

    def _parquet_writer(self, files: List[str], out: str) -> None:
        """_parquet_writer Helper function
        Write a list of parquet files to a single parquet file
        Args:
            files (List[str]): List of parquet files
            out (str): Output file path

        Returns:
            None
        """
        if Path(out).suffix != ".parquet":
            logger.warning(
                f"Output file {out} does not have .parquet suffix altough output format is set to parquet. Fix this inconsistency!"
            )
        fschema: pa.Schema = pq.read_schema(files[0])
        with pq.ParquetWriter(out, schema=fschema) as writer:
            for f in files:
                writer.write_table(pq.read_table(f, schema=fschema))

    def _csv_writer(self, files: List[str], out: str) -> None:
        """_csv_writer Helper function
        Write a list of csv files to a single csv file
        Args:
            files (List[str]): List of csv files
            out (str): Output file path

        Returns:
            None
        """
        gz_suffixes: Set[str] = {".gz", ".gzip", ".bgz"}
        if Path(out).suffix.lower() in gz_suffixes:
            handler: Callable = gzopen
        else:
            handler = open
        schema: pa.Schema = pq.read_schema(files[0])
        write_opts = WriteOptions(
            include_header=False, delimiter="\t", quoting_style="none"
        )
        # The following is a hack until this issue in pyarrow is resolved:
        # https://github.com/apache/arrow/issues/41239
        with handler(out, "wb") as wh:
            wh.write(("\t".join(schema.names) + "\n").encode("utf-8"))
            with CSVWriter(wh, schema=schema, write_options=write_opts) as csvh:
                for f in files:
                    csvh.write_table(pq.read_table(f, schema=schema))


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


def _build_row_map(
    gene_df: pa.Table, ann_cols: Set[str]
) -> Dict[Tuple[int, int], WindowCount]:
    """_build_row_map Helper function
    Build a dictionary of rows for a given gene
    Args:
        gene_df (pa.Table): gene table with crosslink count per window

    Returns:
        Dict[Tuple[int, int], WindowCount]: Dictionary of rows
    """
    row_map: Dict[Tuple[int, int], WindowCount] = {}
    for row in gene_df.to_pylist():
        uid = (row["begin"], row["end"])
        try:
            row_map[uid].add(row["sample"], row["pos_counts"])
        except KeyError:
            row_map[uid] = WindowCount()
            row_map[uid].add(row["sample"], row["pos_counts"])
            annotations: Dict[str, str] = {}
            for an in ann_cols:
                annotations[an] = row[an]
            row_map[uid].annotations = annotations
    return row_map


def _generate_count_schema(
    samples: List[str], unique_id: str = "unique_id"
) -> pa.Schema:
    """_generate_count_schema Helper function
    Generate a schema for the count file
    Args:
        samples (List[str]): List of sample names
        unique_id (str): Unique id column name for the count file

    Returns:
        pa.Schema: Schema for the count file
    """
    return pa.schema(
        [
            pa.field(unique_id, pa.string(), nullable=False),
            *[pa.field(s, pa.uint32(), nullable=False) for s in samples],
        ]
    )


def create_matrices(
    fragment: ds.ParquetFileFragment,
    chrom: str,
    strand: str,
    samples: List[str],
    ann: str,
    sums: str,
    maxs: Optional[str] = None,
    is_windowed: bool = False,
    allow_duplicates: bool = False,
):
    """
    Processes a Parquet fragment to generate annotation, sum, and max count matrices for genomic intervals.

    Args:
        fragment (pyarrow.dataset.ParquetFileFragment): The Parquet fragment containing count data.
        chrom (str): Chromosome name.
        strand (str): Strand information ('+' or '-').
        samples (List[str]): List of sample names.
        ann (str): Output file path for the annotation matrix.
        sums (str): Output file path for the sum matrix.
        maxs (Optional[str], optional): Output file path for the max matrix. Defaults to None.
        is_windowed (bool, optional): Whether the data is windowed. Defaults to False.
        allow_duplicates (bool, optional): Whether to allow duplicate intervals. Defaults to False.

    Returns:
        None
    """
    annotation_cols: Set[str] = set(
        [
            "uniq_id",
            "gene_id",
            "gene_name",
            "gene_type",
            "feature",
            "nr_of_region",
            "total_region",
            "window_number",
        ]
    )
    diff_cols: Set[str] = annotation_cols - set(count_schema.names)
    if len(diff_cols) > 0:
        missing = ", ".join(diff_cols)
        raise RuntimeError(
            f"Missing columns in count files: {missing}!. Please check the count file schema."
        )
    # how to handle duplicates
    if allow_duplicates:
        row_fn: Callable = all_rows
    else:
        row_fn: Callable = pick_rows
    # annotation schema and dictionary
    ann_schema: pa.Schema = get_annotation_schema(windowed=is_windowed)
    ann_dict: Dict[str, List[Any]] = dict([(name, list()) for name in ann_schema.names])
    # count and max count output schema
    out_schema: pa.Schema = _generate_count_schema(
        samples=samples, unique_id="unique_id"
    )
    count_dict: Dict[str, List[str | int]] = dict(
        [(name, list()) for name in out_schema.names]
    )
    # max count dictionary
    max_dict: Dict[str, List[str | int]] = dict(
        [(name, list()) for name in out_schema.names]
    )
    # counts table
    counts: pa.Table = fragment.to_table().sort_by("begin")
    # genes
    genes: List[str] = counts["gene_id"].unique().to_pylist()
    for gene in genes:
        gene_df: pa.Table = counts.filter(pc.field("gene_id") == gene)
        gene_name: str = gene_df.column("gene_name").unique()[0].as_py()
        gene_type: str = gene_df.column("gene_type").unique()[0].as_py()
        row_map: Dict[Tuple[int, int], WindowCount] = _build_row_map(
            gene_df=gene_df,
            ann_cols={
                "uniq_id",
                "feature",
                "nr_of_region",
                "total_region",
                "window_number",
            },
        )
        uniq_rows: List[Tuple[int, int]] = row_fn(row_map=row_map, strand=strand)
        for ur in uniq_rows:
            # fill annotation dictionary
            ann_dict["unique_id"].append(row_map[ur].annotations["uniq_id"])
            ann_dict["chromosome"].append(chrom)
            ann_dict["begin"].append(ur[0])
            ann_dict["end"].append(ur[1])
            ann_dict["strand"].append(strand)
            ann_dict["gene_id"].append(gene)
            ann_dict["gene_name"].append(gene_name)
            ann_dict["gene_type"].append(gene_type)
            ann_dict["gene_region"].append(row_map[ur].annotations["feature"])
            ann_dict["Nr_of_region"].append(row_map[ur].annotations["nr_of_region"])
            ann_dict["Total_nr_of_region"].append(
                row_map[ur].annotations["total_region"]
            )
            if is_windowed:
                ann_dict["window_number"].append(
                    row_map[ur].annotations["window_number"]
                )
            # fill count dictionary and max count dictionary
            cx_sum = row_map[ur].window_sum
            cx_max = row_map[ur].window_max
            count_dict["unique_id"].append(row_map[ur].annotations["uniq_id"])
            max_dict["unique_id"].append(row_map[ur].annotations["uniq_id"])
            for sample in samples:
                if sample in cx_sum:
                    count_dict[sample].append(cx_sum[sample])
                    max_dict[sample].append(cx_max[sample])
                else:
                    count_dict[sample].append(0)
                    max_dict[sample].append(0)
    # write annotation table
    pq.write_table(pa.table(ann_dict, schema=ann_schema), ann)
    # write count table
    pq.write_table(pa.table(count_dict, schema=out_schema), sums)
    # write max count table if provided
    if maxs is not None:
        pq.write_table(pa.table(max_dict, schema=out_schema), maxs)


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
