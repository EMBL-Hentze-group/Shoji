import tempfile
from pathlib import Path
from shutil import rmtree
from typing import Dict, List, Optional

import pyarrow as pa
import pyarrow.csv as pv
import pyarrow.dataset as pds
import pyarrow.parquet as pq
from loguru import logger


def skip_comment(row) -> str:
    """skip_comment comment skipper
    Skip lines starting with "#" in gff and bed formatted files
    Args:
        row: pyarrow row

    Returns:
        "skip" or error
    """
    # comment skipper
    if row.text.startswith("#") or row.text.startswith("track "):
        return "skip"
    return "error"


class Reader:
    """
    Class to read GFF and bed formatted files using pyarrow
    """

    def __init__(self, file_name: str, tmp_dir: Optional[str] = None) -> None:
        """__init__

        Args:
            file_name: BED/GFF file name
            tmp_dir: Temp directory for intermediate files. If None, use system temp. dir.
        """
        self.file_name = file_name
        self.tmp_dir = tmp_dir
        # read options
        self._gff_read = pv.ReadOptions(
            # source: https://genome.ucsc.edu/FAQ/FAQformat.html#format3
            column_names=[
                "seqname",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attributes",
            ]
        )
        self._bed6_read = pv.ReadOptions(
            # source: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
            column_names=["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
        )
        # parse options
        self._general_parse = pv.ParseOptions(
            delimiter="\t", invalid_row_handler=skip_comment
        )
        # convert options
        self._gff_convert = pv.ConvertOptions(
            # "score" column is omitted for the moment as it can also be "."
            column_types={
                "seqname": pa.string(),  # to be compatible with ENSEMBL chromosome names
                "start": pa.uint32(),
                "end": pa.uint32(),
                "score": pa.string(),
                "frame": pa.string(),
            }
        )
        self._bed6_convert = pv.ConvertOptions(
            column_types={
                "chrom": pa.string(),  # to be compatible with ENSEMBL chromosome names
                "chromStart": pa.uint32(),
                "chromEnd": pa.uint32(),
                "score": pa.float32(),
            }
        )
        if not Path(self.file_name).exists():
            raise FileNotFoundError(
                f"Cannot find {file_name}! No such file or directory"
            )

    def _to_partitioned_parquet(
        self,
        root_path: str,
        read_opts: pa._csv.ReadOptions,  # type: ignore
        parse_opts: pa._csv.ParseOptions,  # type: ignore
        convert_opts: pa._csv.ConvertOptions,  # type: ignore
        partition_opts,
    ) -> None:
        """_to_partitioned_parquet
        Helper function
        convert the file to partitioned parquet format
        Args:
            root_path: base file name for the new file
            read_opts: pyarrrow.csv.ReadOptions
            parse_opts: pyarrow.csv.ParseOptions
            convert_opts: pyarrow.csv.ConvertOptions
            partition_opts: pyarrow.dataset.partition options
        """
        with tempfile.NamedTemporaryFile(suffix=".parquet", dir=self.tmp_dir) as tmpq:
            logger.debug(f"Tmp file: {tmpq.name}")
            writer: Optional[pq.core.ParquetWriter] = None
            with pv.open_csv(
                self.file_name,
                read_options=read_opts,
                parse_options=parse_opts,
                convert_options=convert_opts,
            ) as reader:
                for next_chunk in reader:
                    if next_chunk is None:
                        break
                    if writer is None:
                        writer = pq.ParquetWriter(tmpq.name, next_chunk.schema)
                    writer.write_table(pa.Table.from_batches([next_chunk]))
            writer.close()  # type: ignore
            tmpds = pq.ParquetDataset(tmpq, memory_map=True)
            pq.write_to_dataset(
                tmpds.read(), root_path=root_path, partitioning=partition_opts
            )

    def gff(self):
        """gff gff reader
        Read gff file using pyarrow and return a pyarrow table
        Returns:
            pyarrow table
        """
        # TODO: add return type for pyarrow table
        return pv.read_csv(
            self.file_name,
            read_options=self._gff_read,
            parse_options=self._general_parse,
            convert_options=self._gff_convert,
        )

    def gff_to_partitioned_parquet(
        self, root_path: str, flavor: Optional[str] = None
    ) -> None:
        """gff_to_partitioned_parquet gff to partitioned parquet
        Convert GFF file to partitioned parquet format
        Args:
            root_path:  root directory of the new partitioned dataset
            flavor: str, partitioning flavor. If None, use directory partitioning
        """
        # parition options
        # for now the only available paritiioning will be based on chromosome
        gff_partition = pds.partitioning(
            pa.schema([("seqname", pa.string())]), flavor=flavor
        )
        self._to_partitioned_parquet(
            root_path=root_path,
            read_opts=self._gff_read,
            parse_opts=self._general_parse,
            convert_opts=self._gff_convert,
            partition_opts=gff_partition,
        )

    def bed6(self):
        """bed6 bed6 reader
        Read bed6 formatted file using pyarrow

        Returns:
            pyarrow table
        """
        return pv.read_csv(
            self.file_name,
            read_options=self._bed6_read,
            parse_options=self._general_parse,
            convert_options=self._bed6_convert,
        )

    def bed6_to_partitioned_parquet(
        self, root_path: str, flavor: Optional[str] = None
    ) -> None:
        """bed6_to_partitioned_parquet bed6 to partitioned parquet
        Convert BED6 file to partitioned parquet format
        Args:
            root_path:  root directory of the new partitioned dataset
            flavor: str, partitoning flavor. If None, use directory partitioning
        """
        bed_partition = pds.partitioning(
            pa.schema([("chrom", pa.string())]), flavor=flavor
        )
        self._to_partitioned_parquet(
            root_path=root_path,
            read_opts=self._bed6_read,
            parse_opts=self._general_parse,
            convert_opts=self._bed6_convert,
            partition_opts=bed_partition,
        )


class PartionedParquetReader:
    """
    Class to handle GFF/BED6 files in partitoned parquet format
    Given a GFF3/BED6 file, convert the file to partitioned parquet,
    partitiotioning on chromosome (default) and read individual files
    """

    def __init__(self, file_name: str, fformat: str, temp_dir: str, cores: int) -> None:
        self.file_name: str = file_name
        # set number of cpus to use
        pa.set_cpu_count(cores)
        self._pr = Reader(self.file_name, temp_dir)
        self._tmp_pq: str = self._touch_tmp_pq(temp_dir)
        self._is_gff: bool = False
        self._is_bed6: bool = False
        self._is_supported_format(fformat)
        # only "hive" partitioning supports pds.get_partition_keys(fragment.partition_expression)
        self._flavor: str = "hive"
        # partitioned fragment dictionary
        self._fragments: Dict[str, pds.ParquetFileFragment] = {}
        # partitioned dataset
        self._ppq: pds.Dataset = None

    def __enter__(self) -> "PartionedParquetReader":
        logger.debug(
            f"Temp. parquet file: {self._tmp_pq}, partitioning flavor: {self._flavor}"
        )
        if self._is_gff:
            self._pr.gff_to_partitioned_parquet(self._tmp_pq, flavor=self._flavor)
        elif self._is_bed6:
            self._pr.bed6_to_partitioned_parquet(self._tmp_pq, flavor=self._flavor)
        self._partition_data()
        return self

    def __exit__(self, except_type, except_val, except_traceback):
        rmtree(self._tmp_pq)

    def _is_supported_format(self, fformat: str) -> None:
        """_is_supported_format Helper function
        Check if the input file is one of the supported formats
        Args:
            format: str, one of gff, gff3, bed6

        Raises:
            NoeImplementedErro: if the format is not one of the supported ones
        """
        gff_formats = set(
            [
                "gff",
                "gff3",
            ]
        )
        bed6_format = set(["bed6"])
        fformat = fformat.lower()
        if fformat in gff_formats:
            self._is_gff = True
        elif fformat in bed6_format:
            self._is_bed6 = True
        else:
            format_str = ", ".join(sorted(gff_formats | bed6_format))
            raise NotImplementedError(
                f"Format not supported! 'format' MUST BE one of {format_str}, found {fformat}"
            )

    def _touch_tmp_pq(self, temp_dir) -> str:
        """_tmp_pq Helper function
        Create tmp. parquet file
        Args:
            temp_dir: str, tmp. dir name. MUST exist

        Returns:
            temp. parquet file name
        """
        second_tmp_dir = tempfile.mkdtemp(dir=temp_dir)
        return second_tmp_dir  # type: ignore

    def _partition_data(self) -> None:
        """_partition_data Helper function
        For GFF/BED files paritioned on chromosome, partition file based on chromosome name
        Raises:
            RuntimeError:Raise runtime error if there are more than one partition key per fragment
        """
        self._ppq = pds.dataset(
            self._tmp_pq, format="parquet", partitioning=self._flavor
        )
        for frag in self._ppq.get_fragments():
            partition_dict = pds.get_partition_keys(frag.partition_expression)
            if len(partition_dict) != 1:  # type: ignore
                raise RuntimeError(
                    f"Supports only one partition key at the moment. {self._tmp_pq} either the partition keys are empty or there are more than one partition keys!"
                )
            self._fragments[list(partition_dict.values())[0]] = frag

    def get_partitioned_fragments(self) -> Dict[str, pds.ParquetFileFragment]:
        """get_partitioned_fragments partitioned fragment

        Returns:
            Dict[str, pds.ParquetFileFragment], key: chromosome name, value: parquet file fragment for the chromosome
        """
        return self._fragments

    @property
    def contigs(self) -> List[str]:
        """contigs
        List of chromosome names
        Returns:
            List[str]
        """
        return sorted(self._fragments.keys())
