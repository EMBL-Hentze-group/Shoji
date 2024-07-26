import gzip
from functools import partial
from typing import Callable, List, Set

from loguru import logger

from .output import output_writer


class MapToId:
    """
    extract "name" column from the annotation file and map the entries to unique id and print out in tab separated format
    """

    def __init__(self, annotation: str, out: str) -> None:
        self.annotation = annotation
        self.out = out
        # parser
        self._reader = self._parser(annotation)
        # check if input is windowed
        self._is_windowed: bool = False
        self._check_windows()
        # output writer
        self._ow = output_writer(self.out, use_tabix=False, preset="bed")

    @staticmethod
    def _parser(fname: str) -> Callable:
        """_parser_ Helper function
        Return appropriate parser
        Args:
            fname: str, input file name

        Returns:
            Callable, either a gzip file handle or a plain text file handle
        """
        with open(fname, "rb") as _rb:
            if _rb.read(2) == b"\x1f\x8b":
                return partial(gzip.open, mode="rt")
            return partial(open, mode="r")

    def _check_windows(self) -> None:
        """_check_windows Helper function
        Check first 250 lines of the input file to to see whether the name column is windowed or not
        Raises:
            RuntimeError: If the name column is a mix of windowed and non windowed entries
        """
        limit: int = 250
        name_components: Set[int] = set()
        with self._reader(self.annotation) as _ah:
            i: int = 0
            for l in _ah:
                if l[0] == "#":
                    continue
                ldat = l.strip().split("\t")
                name_components.add(len(ldat[3].split("@")))
                i += 1
                if i >= limit:
                    break
        ncomp = sorted(name_components)
        if len(ncomp) > 1:
            raise RuntimeError(
                f"'Name' column in file:{self.annotation} has varying numbers of elements separated by '@'. Check this input file!"
            )
        elif len(ncomp) == 1 and ncomp[0] == 6:
            logger.info(f"{self.annotation} is not windowed")
            self._is_windowed = False
        elif len(ncomp) == 1 and ncomp[0] == 7:
            logger.info(f"{self.annotation} is windowed")
            self._is_windowed = True

    def map_to_id(self) -> None:
        """map_to_id
        extract "name" column from the annotation file and map the entries to unique id and print out in tab separated format
        """

        def last_col(name_dat: List[str]) -> str:
            return name_dat[-1]

        def second_last_col(name_dat: List[str]) -> str:
            return name_dat[-2]

        def add_window_number(out_list: List[str], window_number: str) -> List[str]:
            return out_list + [window_number]

        def add_nothing(out_list: List[str], window_number: str) -> List[str]:
            return out_list

        header: List[str] = [
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
        uniq_id_fn: Callable = last_col
        window_number: Callable = add_nothing
        if self._is_windowed:
            header.append("window_number")
            uniq_id_fn = second_last_col
            window_number = add_window_number
        with self._reader(self.annotation) as _ah, self._ow(self.out) as oh:
            oh.write("\t".join(header) + "\n")
            for l in _ah:
                if l[0] == "#":
                    continue
                ldat: List[str] = l.strip().split("\t")
                name_dat: List[str] = ldat[3].split("@")
                nr_region, total_nr_region = name_dat[4].split("/")
                out_list: List[str] = (
                    [
                        uniq_id_fn(name_dat),
                    ]
                    + ldat[0:3]
                    + [ldat[-1]]
                    + name_dat[0:4]
                    + [nr_region, total_nr_region]
                )
                out_list = window_number(out_list, name_dat[-1])
                oh.write("\t".join(out_list) + "\n")
