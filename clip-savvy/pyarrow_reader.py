import pyarrow as pa
from pyarrow import csv


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

    def __init__(self) -> None:
        self._gff_read = csv.ReadOptions(
            column_names=[
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "attributes",
            ]
        )
        self._bed6_read = csv.ReadOptions(
            column_names=["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
        )
        self._general_parse = csv.ParseOptions(
            delimiter="\t", invalid_row_handler=skip_comment
        )
        self._gff_convert = csv.ConvertOptions(
            column_types={"start": pa.uint32(), "end": pa.uint32()}
        )
        self._bed6_convert = csv.ConvertOptions(
            column_types={
                "chromStart": pa.uint32(),
                "chromEnd": pa.uint32(),
                "score": pa.float32(),
            }
        )

    def gff(self, file_name: str):
        """gff gff reader
        Read gff file using pyarrow and return a pyarrow table
        Args:
            file_name: str, gff file name, can be plain text, .gz or .bz2

        Returns:
            pyarrow table
        """
        # TODO: add return type for pyarrow table
        return csv.read_csv(
            file_name,
            read_options=self._gff_read,
            parse_options=self._general_parse,
            convert_options=self._gff_convert,
        )

    def bed6(self, file_name: str):
        """bed6 bed6 reader
        Read bed6 formatted file using pyarrow
        Args:
            file_name: str, bed6 file name, can be plain text, .gz or .bz2

        Returns:
            pyarrow table
        """
        return csv.read_csv(
            file_name,
            read_options=self._bed6_read,
            parse_options=self._general_parse,
            convert_options=self._bed6_convert,
        )
