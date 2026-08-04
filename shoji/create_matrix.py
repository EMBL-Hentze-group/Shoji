import tempfile
from collections import defaultdict
from collections.abc import Callable, Generator
from dataclasses import dataclass, field
from pathlib import Path
from shutil import rmtree
from typing import Self

import polars as pl
import polars.selectors as cs
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from loguru import logger

from .schemas import count_schema, get_annotation_schema


@dataclass
class Annotations:
    gene_ids: set[str] = field(default_factory=set)
    uniq_ids: set[str] = field(default_factory=set)

    @property
    def gene_id_count(self) -> int:
        return len(self.gene_ids)

    @property
    def uniq_id_count(self) -> int:
        return len(self.uniq_ids)


class CreateMatrix:
    def __init__(
        self,
        in_dir: str,
        annotation: str,
        out: str,
        max_out: str | None = None,
        cores: int = 1,
        prefix: str | None = None,
        suffix: str | None = ".parquet",
        format: str = "csv",
        tmp_dir: str | None = None,
    ) -> None:
        """CreateMatrix Class
        Args:
            in_dir (str): Input directory with count files
            annotation (str): Output file for the annotation matrix
            out (str): Output file for the sum matrix
            max_out (str | None, optional): Output file for the max count matrix. Defaults to None.
            cores (int, optional): Number of cores to use. Defaults to 1.
            prefix (str | None, optional): Prefix for the input files. Defaults to None.
            suffix (str | None, optional): Suffix for the input files. Defaults to ".parquet".
            format (str, optional): Format of the output files. Defaults to "csv", could be either "csv" or "parquet".
            tmp_dir (str | None, optional): Temporary directory for intermediate files. Defaults to None.
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
        ####
        # Place holders for various attributes
        ####
        self._count_files: list[Path] = []  # count files
        self._samples: list[str] = []  # sample names
        self._count_ds: ds.Dataset = None  # placeholder for count dataset
        self._is_windowed: bool = False  #  windowed flag
        self._gene_id_map: dict[tuple[str, str], Annotations] = defaultdict(
            Annotations
        )  # gene id map [chrom, strand] = Annotations(gene_ids=set(), uniq_ids=set())
        self._partition_schema = ds.partitioning(
            schema=pa.schema([("chrom", pa.string()), ("strand", pa.string())]),
            flavor="hive",
        )  # partition schema for the dataset
        ####
        # Place holders for dictionary encoded data
        ####
        self._new_schema: pa.Schema = (
            None  # new schema for count dataset with dictionary encoded columns
        )

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
        pa.set_cpu_count(cores)
        pa.set_io_thread_count(cores)

    def __enter__(self) -> Self:
        self._glob_count_files()
        self._validate_count_schema()
        self._validate_data()
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
        logger.info(
            f"Found {len(self._count_files)} count files in {self.in_dir}"
        )

    def _validate_count_schema(self):
        """_validate_count_schema Helper function
        Validates that all count files have the same schema as the expected count_schema
        Raises:
            RuntimeError: If any of the count files have a schema that does not match the expected count_schema
        """
        for count_file in self._count_files:
            file_schema = pq.read_schema(count_file)
            if not file_schema.equals(count_schema):
                raise RuntimeError(
                    f"Count file {count_file} has incompatible schema.\n Expected \n{count_schema}\n found \n{file_schema}\n Check input files!"
                )
        logger.info(
            "Finished validating count file schema. All files are compatible with expected schema."
        )

    def _validate_data(self) -> None:
        """_validate_data Helper function
        Checks whether all input files have unique sample names
        checks whether all input files have the same annotation length
        Parse gene id and unique id information from the count files and store in self._gene_id_map

        Raises:
            RuntimeError: If the number of unique sample_names and number of files in self._count_files does not match
            RuntimeError: If the number of unique annotation elements are not the same across all samples
        Returns:
            None
        """

        samples: set[str] = set()
        has_any_null: bool = False

        self._count_ds = ds.dataset(
            self._count_files, format="parquet", schema=count_schema
        )
        logger.info(
            "Validating sample consistency and windowed annotation for all count files. This may take a while for large datasets..."
        )
        meta_scanner: ds.Scanner = self._count_ds.scanner(
            columns=[
                "chrom",
                "strand",
                "gene_id",
                "uniq_id",
                "window_number",
                "sample",
            ]
        )
        for batch in meta_scanner.to_batches():
            samples.update(pc.unique(batch["sample"]).to_pylist())
            has_any_null = (
                has_any_null or pc.any(batch["window_number"].is_null()).as_py()
            )
            # gene_ids.update(pc.unique(batch["gene_id"]).to_pylist())
            unique_dat = (
                pa.Table.from_batches([batch])
                .group_by(["chrom", "strand"])
                .aggregate([("gene_id", "distinct"), ("uniq_id", "distinct")])
                .to_pylist()
            )
            for ud in unique_dat:
                cid = (ud["chrom"], ud["strand"])
                self._gene_id_map[cid].gene_ids.update(ud["gene_id_distinct"])
                self._gene_id_map[cid].uniq_ids.update(ud["uniq_id_distinct"])
        # check if all files have same annotation schema
        if has_any_null:
            raise RuntimeError(
                "Count files appear to be a mix of windowed and non windowed data! Check input files"
            )
        # is windowed if any of the files have a non-null window_number column
        self._is_windowed = True
        logger.info("Count files are windowed")
        self._samples = sorted(samples)
        # number of sample names MUST match number of input files!
        if len(self._samples) != len(self._count_files):
            raise RuntimeError(
                f"Mismatch in number of samples! found {len(self._count_files)} in {self.in_dir} but found {len(self._samples)} sample names in merged count files"
            )
        # check for gene id map
        if len(self._gene_id_map) == 0:
            raise RuntimeError(
                f"Cannot parse chromosome name, strand and gene ids from {self._count_files}! Check input files"
            )

    def _generate_schema_and_columns(self) -> None:
        """_generate_schema_and_columns Helper function
        Generate the new schema and column expressions for the count dataset with dictionary encoded columns.

        Returns:
            None
        """
        # @TODO: This is a bit of a hack, but it works for now. Refactor this to be more elegant
        encoding_map: dict[str, pa.DataType] = {
            "sample": pa.dictionary(
                index_type=self._get_pyarrow_type(len(self._samples)),
                value_type=pa.string(),
            ),
            "gene_id": pa.dictionary(
                index_type=self._get_pyarrow_type(
                    max([v.gene_id_count for v in self._gene_id_map.values()])
                ),
                value_type=pa.string(),
            ),
            "uniq_id": pa.dictionary(
                index_type=self._get_pyarrow_type(
                    max([v.uniq_id_count for v in self._gene_id_map.values()])
                ),
                value_type=pa.string(),
            ),
        }
        mod_fields: list[pa.Field] = []
        for ds_field in self._count_ds.schema:
            if ds_field.name in encoding_map:
                mod_fields.append(
                    pa.field(
                        ds_field.name,
                        encoding_map[ds_field.name],
                        nullable=ds_field.nullable,
                    )
                )
            else:
                mod_fields.append(ds_field)

        self._new_schema = pa.schema(mod_fields)

    def _get_pyarrow_type(self, val: int) -> pa.DataType:
        """_get_pyarrow_type Helper function
        Get the appropriate pyarrow type for a given integer value

        Args:
            val (int): The integer value to determine the pyarrow type for

        Returns:
            pa.DataType: The corresponding pyarrow data type
        """
        if not isinstance(val, int):
            raise TypeError(f"Expected int, got {type(val)}")
        if val <= 0:
            raise ValueError(f"Expected positive int, got {val}")
        if val < 2**8 - 1:
            return pa.uint8()
        elif val < 2**16 - 1:
            return pa.uint16()
        elif val < 2**32 - 1:
            return pa.uint32()
        else:
            return pa.uint64()

    def _prepare_dataset(self) -> None:
        """_prepare_dataset Helper function
        Generate partitioned parquet dataset from the count files
        Returns:
            None
        """
        logger.debug(
            f"Writing dataset partitioned by chrom and strand to {self._partitioned_ds}"
        )
        self._generate_schema_and_columns()
        mapping: dict[str, pa.Array] = {
            "sample": pa.array(self._samples, type=pa.string()),
        }
        # partition per chromosome and strand
        for (chrom, strand), annotations in self._gene_id_map.items():
            ds_filter: pc.Expression = (pc.field("chrom") == chrom) & (
                pc.field("strand") == strand
            )  # filter expression for the current chromosome and strand
            mapping["gene_id"] = pa.array(
                sorted(annotations.gene_ids), type=pa.string()
            )  # gene id mapping
            mapping["uniq_id"] = pa.array(
                sorted(annotations.uniq_ids), type=pa.string()
            )  # unique id mapping
            # stream and encode the projected dataset
            ds.write_dataset(
                data=self._stream_and_project(ds_filter, mapping),
                base_dir=self._partitioned_ds,
                format="ipc",
                schema=self._new_schema,
                existing_data_behavior="overwrite_or_ignore",
                partitioning=self._partition_schema,
            )

    def _stream_and_project(
        self,
        filter_expr: pc.Expression,
        mapping: dict[str, pa.Array],
    ) -> Generator[pa.RecordBatch, None, None]:
        projections: dict[str, pc.Expression] = {}  # projections
        for ds_field in self._count_ds.schema:
            if ds_field.name in mapping:
                projections[ds_field.name] = pc.index_in(  # type: ignore
                    pc.field(ds_field.name), value_set=mapping[ds_field.name]
                )
            else:
                projections[ds_field.name] = pc.field(ds_field.name)
        scanner: ds.Scanner = self._count_ds.scanner(
            filter=filter_expr, columns=projections
        )
        for batch in scanner.to_batches():
            new_columns: list[pa.Array] = []
            for name, col_data in zip(batch.schema.names, batch.columns):
                if name in mapping:
                    encoded_col: pa.DictionaryArray = (
                        pa.DictionaryArray.from_arrays(
                            indices=col_data, dictionary=mapping[name]
                        )
                    )
                    new_columns.append(encoded_col)
                else:
                    new_columns.append(col_data)
            yield pa.RecordBatch.from_arrays(
                new_columns, schema=self._new_schema
            )

    def create_matrices(self, allow_duplicates: bool = False) -> None:
        logger.info(
            f"Processing count files per chromosome and strand. Filter overlapping windows with identical crosslink counts: {not allow_duplicates}"
        )
        partitioned_ds = ds.dataset(
            str(self._partitioned_ds), format="ipc", partitioning="hive"
        )
        fragments: dict[tuple[str, str], str] = self._get_fragment_map(
            partitioned_ds
        )
        ann_files, out_files, max_files = [], [], []
        for chrom, strand in sorted(fragments.keys()):
            tmp_name: str = next(tempfile._get_candidate_names())  # type: ignore
            ann_file: str = str(self._tmp / f"ann_{chrom}_{tmp_name}.ipc")
            ann_files.append(ann_file)
            out_file: str = str(self._tmp / f"count_{chrom}_{tmp_name}.ipc")
            out_files.append(out_file)
            max_file: str | None = None
            if self.max_out is not None:
                max_file = str(self._tmp / f"maxs_{chrom}_{tmp_name}.ipc")
                max_files.append(max_file)
            create_matrices(
                fragments[(chrom, strand)],
                chrom,
                strand,
                self._samples,
                ann_file,
                out_file,
                max_file,
                self._is_windowed,
                allow_duplicates,
            )
            n_genes: int = self._gene_id_map[(chrom, strand)].gene_id_count
            n_uniq: int = self._gene_id_map[(chrom, strand)].uniq_id_count
            logger.info(
                f"Finished chromosome {chrom} strand {strand}, genes: {n_genes:,}  windows: {n_uniq:,}"
            )
        logger.info("Writing output files")
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

    def _get_fragment_map(
        self, partitioned_ds: ds.Dataset
    ) -> dict[tuple[str, str], str]:
        """_get_fragment_map Helper function
        Get a mapping of (chrom, strand) to fragment file path from the partitioned dataset

        Args:
            partitioned_ds (ds.Dataset): The partitioned dataset
        Returns:
            dict[tuple[str, str], str]: A mapping of (chrom, strand) to fragment file path
        """
        fragments: dict[tuple[str, str], str] = {}
        for fragment in partitioned_ds.get_fragments():
            partition_dict = ds.get_partition_keys(
                fragment.partition_expression
            )
            fragments[(partition_dict["chrom"], partition_dict["strand"])] = (
                fragment.path
            )
        return fragments

    def _parquet_writer(self, files: list[str], out: str) -> None:
        """_parquet_writer Helper function
        Write a list of parquet files to a single parquet file
        Args:
            files (list[str]): List of parquet files
            out (str): Output file path

        Returns:
            None
        """
        if Path(out).suffix != ".parquet":
            logger.warning(
                f"Output file {out} does not have .parquet suffix altough output format is set to parquet. Fix this inconsistency!"
            )
        pl.scan_ipc(files).sink_parquet(out)

    def _csv_writer(self, files: list[str], out: str) -> None:
        """_csv_writer Helper function
        Write a list of csv files to a single csv file
        Args:
            files (list[str]): List of csv files
            out (str): Output file path

        Returns:
            None
        """
        gz_suffixes: set[str] = {".gz", ".gzip", ".bgz"}
        if Path(out).suffix.lower() in gz_suffixes:
            compression: str = "gzip"
        else:
            compression: str = "uncompressed"
        pl.scan_ipc(files).sink_csv(
            out, separator="\t", compression=compression
        )


def _strand_specific_row_picker(
    strand: str,
) -> Callable[[pl.LazyFrame, list[str], str], pl.LazyFrame]:
    """
    Return a callable to pick the most 5' window from a set of overlapping winodws with same crosslink count over samples

    Args:
        strand (str): Strand information ('+' or '-').
    Returns:
        Callable[[pl.LazyFrame, list[str], str], pl.LazyFrame]: A callable function that takes a LazyFrame, a list of data columns, and a picker column, and returns a deduplicated LazyFrame with the most 5' window selected based on the strand information
    """

    def _plus_strand_picker(
        uid_agg: pl.LazyFrame,
        data_cols: list[str] = ["gene_id", "sample", "pos_counts"],
        picker_col: str = "begin",
    ) -> pl.LazyFrame:
        """
        Function to pick most 5' prime interval from overlapping intervals with identical crosslink data for the plus strand for plus strand oriented windows.
        Args:
            uid_agg (pl.LazyFrame): polars LazyFrame containing overlapping count data
            data_cols (list[str], optional): columns to check for duplicate data. Defaults to ["sample", "pos_counts"].
            picker_col (str, optional): Column to select most 5' prime value from. Defaults to "begin".

        Returns:
            pl.LazyFrame: Deduplicated LazyFrame
        """
        return (
            uid_agg.with_columns(grp_id=pl.struct(data_cols).rank("dense"))
            .filter(
                pl.col(picker_col) == pl.col(picker_col).min().over("grp_id")
            )
            .drop("grp_id")
        )

    def _minus_strand_picker(
        uid_agg: pl.LazyFrame,
        data_cols: list[str] = ["gene_id", "sample", "pos_counts"],
        picker_col: str = "begin",
    ) -> pl.LazyFrame:
        """
        Function to pick most 5' prime interval from overlapping intervals with identical crosslink data for the minus strand for minus strand oriented windows.
        Args:
            uid_agg (pl.LazyFrame): polars LazyFrame containing overlapping count data
            data_cols (list[str], optional): columns to check for duplicate data. Defaults to ["sample", "pos_counts"].
            picker_col (str, optional): Column to select most 5' prime value from. Defaults to "begin".

        Returns:
            pl.LazyFrame: Deduplicated LazyFrame
        """
        return (
            uid_agg.with_columns(grp_id=pl.struct(data_cols).rank("dense"))
            .filter(
                pl.col(picker_col) == pl.col(picker_col).max().over("grp_id")
            )
            .drop("grp_id")
        )

    if strand == "-":
        return _minus_strand_picker
    else:
        return _plus_strand_picker


def _write_count_matrix(
    count_matrix: pl.LazyFrame,
    out_file: str,
    rename_schema: dict[str, str] | None = None,
) -> None:
    """
    Helper function to write count matrix or max count matrix to disk in parquet format, with optional renaming of columns and sorting.
    Args:
        count_matrix (pl.LazyFrame): Count matrix to be written to disk.
        out_file (str): Output file path for the count matrix
        rename_schema (dict[str, str] | None, optional): Optional schema to rename columns. Defaults to None.
    """
    if rename_schema is not None:
        count_matrix = count_matrix.rename(rename_schema)
    count_matrix.sink_ipc(out_file)


def _write_annotation_matrix(
    counts: pl.LazyFrame,
    out_file: str,
    chrom: str,
    strand: str,
    unique_id_col: str,
    uniq_ids: pl.Series,
    annotation_cols: list[str],
    sort_col: str,
    rename_schema: dict[str, str] | None,
    annotation_order: list[str] | None,
) -> None:
    """
    Helper function to write the annotation matrix to disk in parquet format.

    Args:
        counts (pl.LazyFrame): The counts LazyFrame containing the annotation data.
        out_file (str): Output file path for the annotation matrix.
        uniq_ids (pl.Series): Series of unique IDs to filter the counts.
        annotation_cols (list[str]): List of columns to include in the annotation matrix.
        rename_schema (dict[str, str]): Schema to rename columns in the annotation matrix.
        annotation_order (list[str]): Order of columns in the final annotation matrix.
        chrom (str): Chromosome name.
        strand (str): Strand information ('+' or '-').

    Returns:
        None
    """
    annot: pl.LazyFrame = (
        counts.filter(pl.col(unique_id_col).is_in(uniq_ids.implode()))
        .select(annotation_cols)
        .unique()
        .with_columns(chromosome=pl.lit(chrom), strand=pl.lit(strand))
    )
    if rename_schema is not None:
        annot = annot.rename(rename_schema)
    if annotation_order is not None:
        annot = annot.select(annotation_order)
    annot.sort(sort_col).with_columns(
        cs.categorical().cast(pl.String)
    ).sink_ipc(out_file)


def create_matrices(
    fragment: str,
    chrom: str,
    strand: str,
    samples: list[str],
    ann: str,
    sums: str,
    maxs: str | None = None,
    is_windowed: bool = False,
    allow_duplicates: bool = False,
) -> None:
    """
    Processes a Parquet fragment to generate annotation, sum, and max count matrices for genomic intervals.

    Args:
        fragment (str): The Parquet fragment file path containing count data.
        chrom (str): Chromosome name.
        strand (str): Strand information ('+' or '-').
        samples (list[str]): List of sample names.
        ann (str): Output file path for the annotation matrix.
        sums (str): Output file path for the sum matrix.
        maxs (str | None, optional): Output file path for the max matrix. Defaults to None.
        is_windowed (bool, optional): Whether the data is windowed. Defaults to False.
        allow_duplicates (bool, optional): Whether to allow duplicate intervals. Defaults to False.

    Returns:
        None
    """
    if is_windowed:
        annotation_cols: set[str] = {
            "begin",
            "end",
            "uniq_id",
            "gene_id",
            "gene_name",
            "gene_type",
            "feature",
            "nr_of_region",
            "total_region",
            "window_number",
        }
    else:
        annotation_cols: set[str] = {
            "begin",
            "end",
            "uniq_id",
            "gene_id",
            "gene_name",
            "gene_type",
            "feature",
            "nr_of_region",
            "total_region",
        }
    rename_schema: dict[str, str] = {
        "uniq_id": "unique_id",
        "feature": "gene_region",
        "nr_of_region": "Nr_of_region",
        "total_region": "Total_nr_of_region",
    }
    counts_rename_schema: dict[str, str] = {
        "uniq_id": "unique_id",
    }
    annotation_order: list[str] = get_annotation_schema(
        windowed=is_windowed
    ).names
    # for count and max count matrices
    target_type: pl.DataType = pl.List(
        pl.Struct([pl.Field("key", pl.UInt32), pl.Field("value", pl.UInt32)])
    )
    # compute sum
    xlink_sum: pl.Expr = (
        cs.by_dtype(target_type)
        .list.eval(pl.element().struct.field("value").sum())
        .list.get(0)
        .fill_null(0)
    )
    # compute max count
    xlink_max: pl.Expr = (
        cs.by_dtype(target_type)
        .list.eval(pl.element().struct.field("value").max())
        .list.get(0)
        .fill_null(0)
    )
    # row picker function based on strand
    row_picker: Callable = _strand_specific_row_picker(strand)
    # counts table
    counts: pl.LazyFrame = pl.scan_ipc(fragment)
    diff_cols: set[str] = annotation_cols - set(counts.collect_schema().names())
    if len(diff_cols) > 0:
        # sanity check failed, missing columns in count files
        missing = ", ".join(diff_cols)
        raise RuntimeError(
            f"Missing columns in count files: {missing}!. Please check the count file schema."
        )
    per_gene_query: pl.LazyFrame = counts.group_by(
        ["begin", "end", "gene_id", "uniq_id"]
    ).agg(pl.col("sample").sort(), pl.col("pos_counts").sort_by("sample"))
    if not allow_duplicates:
        # no duplicated info windows allowed, pick the most 5' of overlapping windows
        per_gene_query = row_picker(per_gene_query)
    per_sample_query: pl.LazyFrame = (
        per_gene_query.explode(["sample", "pos_counts"], empty_as_null=True)
        .pivot(
            on="sample",
            on_columns=samples,
            index=["uniq_id", "begin"],
            values="pos_counts",
            aggregate_function="first",
        )
        .sort("begin")
        .drop("begin")
    )
    # write count matrix
    _write_count_matrix(
        count_matrix=per_sample_query.with_columns(xlink_sum),
        out_file=sums,
        rename_schema=counts_rename_schema,
    )
    if maxs is not None:
        # write max count matrix
        _write_count_matrix(
            count_matrix=per_sample_query.with_columns(xlink_max),
            out_file=maxs,
            rename_schema=counts_rename_schema,
        )
    # write annotation matrix
    uniq_ids: pl.Series = (
        per_sample_query.select("uniq_id").sort("uniq_id").collect()["uniq_id"]
    )
    _write_annotation_matrix(
        counts=counts,
        out_file=ann,
        chrom=chrom,
        strand=strand,
        unique_id_col="uniq_id",
        uniq_ids=uniq_ids,
        sort_col="begin",
        annotation_cols=list(annotation_cols),
        rename_schema=rename_schema,
        annotation_order=annotation_order,
    )
