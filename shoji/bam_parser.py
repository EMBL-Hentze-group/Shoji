import logging
import tempfile
from decimal import ROUND_HALF_UP, Decimal
from itertools import chain
from pathlib import Path
from shutil import rmtree
from typing import Callable, List, Optional, Set, Tuple

import numpy as np
import pysam
from sortedcontainers import SortedList

from .helpers import set_cores
from .output import general_accumulator, output_writer, tabix_accumulator

logger = logging.getLogger(__name__)


class BamParser:
    """Class to parse Bam files.
    Bam file must be co-ordinate sorted and indexed
    """

    def __init__(
        self,
        bam: str,
        output: str,
        use_tabix: bool,
        cores: int,
        tmp_dir: Optional[str] = None,
    ) -> None:
        """__init__ _summary_

        Args:
            bam: bam file to parse. Must be co-ordinate sorted and indexed
            output: Output file name (bed format)
            use_tabix: boolean, if True, use tabix to index the output file
            cores: int, number of cores to use
            tmp_dir: Tmp. directory to store intermediate outputs. Defaults to None.
        """
        self.bam: str = bam
        self.output: str = output
        self.use_tabix: bool = use_tabix
        self.cores: int = set_cores(cores)
        if tmp_dir is None:
            self.tmp: Path = Path(self.output).parent / next(
                tempfile._get_candidate_names()
            )
        else:
            if not Path(tmp_dir).exists():
                raise FileNotFoundError(f"Directory {tmp_dir} does not exist")
            self.tmp: Path = Path(tmp_dir) / next(tempfile._get_candidate_names())
        logger.info("Using %s as temp. directory", str(self.tmp))
        # list of chromosomes in the bam file
        self._chroms: List[str] = self._get_chromosomes()
        logger.info("Found %i chromosomes in %s", len(self._chroms), self.bam)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        rmtree(self.tmp)  # clean up
        if exc_type:
            logger.exception(traceback)

    def _get_chromosomes(self) -> List[str]:
        """_get_chromosomes Helper function
        Check if the input bam file is coordinate sorted and indexed.
        Collect chromosome names from bam header
        """
        with pysam.AlignmentFile(
            self.bam, mode="rb", check_sq=True, require_index=True
        ) as _bam:
            header = dict(_bam.header)  # type: ignore
            if "SQ" not in header:  # type: ignore
                msg = f"Cannot find @SQ header lines in file {self.bam}!"
                logger.error(msg)
                raise LookupError(msg)
        return sorted(map(lambda s: s["SN"], header["SQ"]))  # type: ignore

    def extract(
        self,
        site: str,
        mate: int,
        offset: int,
        min_qual: int,
        min_len: int,
        max_len: int,
        max_interval_len: int,
        primary: bool,
        ignore_PCR: bool,
    ) -> None:
        """extract Extract crosslink sites
        Extract crosslink sites from bam file based on user defined parameters
        and map it to "site"
        Args:
            site: str, Crosslink site choices, must be one of ["s", "i", "d", "m", "e"]
            mate: int, Select the read/mate to extract the crosslink sites from. Must be one of [1,2]
            offset: int, Number of nucleotides to offset for crosslink sites
            min_qual: int, Minimum alignment quality
            min_len: int, Minimum read length
            max_len: int, Maximum read length
            max_interval_len: int, Maximum interval length, for paired end reads, splice length otherwise
            primary: bool, Flag to extract only primary alignments
            ignore_PCR: bool, Flag to ignore PCR duplicates (only if bam file has PCR duplicate flag in alignment)
        """
        pass


def extract_single_site(
    bam: str,
    chrom: str,
    output: str,
    site: str,
    mate: int,
    offset: int,
    min_qual: int,
    min_len: int,
    max_len: int,
    max_interval_len: int,
    primary: bool,
    ignore_PCR: bool,
) -> None:
    """extract_single_site
    Extract one crosslink site/event per read.
    Could be either start site, middle site or end site
    Insertion and deletion sites are not considered here as there can be more than one such event per read
    Args:
        bam: str, BAM file to parse, must be co-ordinate sorted and indexed
        chrom: str, Chromosome name to extract reads
        output: str, tmp. output file name
        site: str, site to extract, could be either s (start), m (middle) or e (end)
        mate: int, mate to extract the crosslink sites from. Must be one of [1,2]
        offset: int, offset start and end sites by "offset" base pairs
        min_qual: int, minimum alignment quality
        min_len: int, minimum read length
        max_len: int, maximum read length
        max_interval_len: int, maximum interval length, for paired end reads, splice length otherwise
        primary: bool, flag to extract only primary alignments
        ignore_PCR: bool, flag to ignore PCR duplicates (only if bam file has PCR duplicate flag in alignment)
    """
    extract_fn = _site_ops(site)
    positions: SortedList = SortedList()
    with pysam.AlignmentFile(bam, mode="rb") as _bam:
        for aln in _bam.fetch(chrom, multiple_iterators=True):
            if _discard_read(
                aln,
                mate,
                min_qual,
                min_len,
                max_len,
                max_interval_len,
                primary,
                ignore_PCR,
            ):
                continue
            start, end = extract_fn(aln, offset)
            strand: str = "-" if aln.is_reverse else "+"
            try:
                yb = aln.get_tag("YB")
            except KeyError:
                yb = 1
            positions.add(
                (start, end, f"{aln.query_name}|{aln.query_length}", yb, strand)
            )
    _tmp_output_writer(output, chrom, positions)


def extract_multiple_sites(
    bam: str,
    chrom: str,
    output: str,
    site: str,
    mate: int,
    offset: int,
    min_qual: int,
    min_len: int,
    max_len: int,
    max_interval_len: int,
    primary: bool,
    ignore_PCR: bool,
) -> None:
    extract_fn = _site_ops(site)
    positions: SortedList = SortedList()
    with pysam.AlignmentFile(bam, mode="rb") as _bam:
        for aln in _bam.fetch(chrom, multiple_iterators=True):
            if _discard_read(
                aln,
                mate,
                min_qual,
                min_len,
                max_len,
                max_interval_len,
                primary,
                ignore_PCR,
            ):
                continue
            strand: str = "-" if aln.is_reverse else "+"
            try:
                yb = aln.get_tag("YB")
            except KeyError:
                yb = 1
            pos_list: List[Tuple[int, int]] = extract_fn(aln)
            for start, end in pos_list:
                positions.add(
                    (start, end, f"{aln.query_name}|{aln.query_length}", yb, strand)
                )
    _tmp_output_writer(output, chrom, positions)


def _discard_read(
    aln: pysam.AlignedSegment,
    mate: int,
    qual: int,
    min_len: int,
    max_len: int,
    max_interval_len: int,
    primary: bool,
    ignore_PCR: bool,
) -> bool:
    """_discard_read Helper function
    Discard reads based on input criteria

    Args:
        aln: pysam.AlignedSegment, Aligned read
        mate: int, Select the read/mate to extract the crosslink sites from. Must be one of [1,2]
        qual: int, Minimum alignment quality
        min_len: int, Minimum read length
        max_len: int, Maximum read length
        max_interval_len: int, Maximum interval length, for paired end reads, splice length otherwise
        primary: bool, Flag to extract only primary alignments
        ignore_PCR: bool, Flag to ignore PCR duplicates (only if bam file has PCR duplicate flag in alignment)

    Returns:
        bool
    """
    if (
        aln.is_unmapped
        or aln.is_qcfail
        or aln.mapping_quality < qual
        or aln.query_length < min_len
        or aln.query_length > max_len
        or aln.reference_length > max_interval_len
    ):
        return True
    if (primary and aln.is_secondary) or (ignore_PCR and aln.is_duplicate):
        return True
    if mate == 1 and aln.is_read2:
        return True
    elif mate == 2 and aln.is_read1:
        return True
    return False


def _start(aln: pysam.AlignedSegment, offset: int) -> Tuple[int, int]:
    """_start Helper function
    Get stranded genomic start position of the read and offset it by "offset" base pairs
    Args:
        aln: pysam.AlignedSegment, Aligned read
        offset: int, offset to add to the start position

    Returns:
        Tuple[int, int], offset start coordinates
    """
    if aln.is_reverse:
        start: int = aln.reference_end - offset
        return start - 1, start
    # aln.reference_start: 0-based leftmost coordinate
    begin0 = aln.reference_start + offset
    return begin0, begin0 + 1


def _end(aln: pysam.AlignedSegment, offset: int) -> Tuple[int, int]:
    """_end Helper function
    Get stranded genomic end position of the read and offset it by "offset" base pairs
    Args:
        aln: pysam.AlignedSegment, Aligned read
        offset: int, offset to add to the start position

    Returns:
        Tuple[int, int], offset end coordinates
    """
    if aln.is_reverse:
        # aln.reference_start: 0-based leftmost coordinate
        begin0 = aln.reference_start + offset
        return begin0, begin0 + 1
    start: int = aln.reference_end - offset
    return start - 1, start


def _middle(aln: pysam.AlignedSegment, offset: int) -> Tuple[int, int]:
    """_middle Helper function
    Get stranded genomic middle position of the aligned part of the read
    offset is simply a placeholder for offset position, not used
    Args:
        aln: pysam.AlignedSegment, Aligned read
        offset: int, offset placeholder

    Returns:
        Tuple[int,int], middle site position
    """
    match_ops: Set[int] = {0, 4, 5}
    cigars: Set[int] = set([x[0] for x in aln.cigartuples])
    mid: int = int(Decimal(aln.query_alignment_length / 2).quantize(0, ROUND_HALF_UP))

    if (len(cigars - match_ops) == 0) and (aln.is_reverse):
        # only match, soft clip, hard clip, negative strand
        middle: int = aln.reference_end - mid
        return middle - 1, middle
    elif (len(cigars - match_ops) == 0) and (not aln.is_reverse):
        # only match, soft clip, hard clip, positive strand
        middle: int = aln.reference_start + mid
        return middle - 1, middle
    else:
        # insertion, deletion or splice operations
        pairs: np.ndarray = np.array(aln.get_aligned_pairs())
        try:
            qmid: np.ndarray = np.where(pairs[:, 0] == mid)[0]
        except IndexError:
            return -1, -1
        if len(qmid) > 1:
            logger.warning("Multiple mid points found for read %s!", aln.query_name)
        rmid = pairs[qmid[0], 1]
        if rmid is None:
            logger.warning("Mid point not found for read %s!", aln.query_name)
            return -1, -1
        if aln.is_reverse:
            return rmid - 2, rmid - 1
        return rmid - 1, rmid


def _insertion(aln: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """_insertion Helper function
    Get insertion points from the read
    Args:
        aln: pysam.AlignedSegment, Aligned read

    Returns:
        List[Tuple[int, int]], List of insertion points (start, end)
    """
    cigars: Set[int] = set([x[0] for x in aln.cigartuples])
    if 1 not in cigars:
        # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
        # No insertion found
        return []
    return _insertion_deletion_points(aln.cigartuples, aln.get_aligned_pairs(), 1)


def _deletion(aln: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """_deletion Helper function
    Get deletion points from the read
    Args:
        aln: pysam.AlignedSegment, Aligned read

    Returns:
        List[Tuple[int, int]], List of deletion points (start, end)
    """
    cigars: Set[int] = set([x[0] for x in aln.cigartuples])  # type: ignore
    if 2 not in cigars:
        # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
        # No deletion found
        return []
    return _insertion_deletion_points(aln.cigartuples, aln.get_aligned_pairs(), 2)  # type: ignore


def _insertion_deletion_points(
    cigartuples: List[Tuple[int, int]], pairs: List[Tuple[int, int]], ops: int
) -> List[Tuple[int, int]]:
    """_insertion_deletion_points Helper function
    Return list of insertion or deletion points from the read
    Args:
        cigartuples: List[Tuple[int, int]], see aln.cigartuples
        pairs: List[Tuple[int, int]], see aln.get_aligned_pairs()
        ops: int, 1 for insertion, 2 for deletion

    Returns:
        List[Tuple[int, int]], List of insertion or deletion points (start, end)
    """
    lcigars: List[int] = list(chain(*[[op] * c for op, c in cigartuples]))
    cigar_ref: np.ndarray = np.array(list(zip(lcigars, [i[1] for i in pairs])))
    ops_pos: np.ndarray = np.where(cigar_ref[:, 0] == ops)[0]
    ops_locations: List[Tuple[int, int]] = []
    for ops_index in np.split(ops_pos, np.where(np.diff(ops_pos) != 1)[0] + 1):
        if len(ops_index) == 0:
            continue
        start_index = np.min(ops_index) - 1
        end_index = np.max(ops_index) + 1
        ops_locations.append((cigar_ref[start_index, 1], cigar_ref[end_index, 1]))  # type: ignore
    return ops_locations


def _site_ops(site: str) -> Callable:
    """_pos_ops Helper function
    Return the appropriate function based on the choice of crosslink site
    Args:
        site: str, one of ["s", "i", "d", "m", "e"]

    Raises:
        NotImplementedError: If site is not one of ["s", "i", "d", "m", "e"]

    Returns:
        Callable, appropriate function based on the site
    """
    sites: Set[str] = {"s", "i", "d", "m", "e"}
    if site not in sites:
        raise NotImplementedError(f"Site must be one of {site}, but found {sites}")
    if site == "s":
        return _start
    elif site == "e":
        return _end
    elif site == "m":
        return _middle
    elif site == "i":
        return _insertion
    else:
        return _deletion


def _tmp_output_writer(
    output: str,
    chrom: str,
    positions: SortedList,
) -> None:
    """_tmp_output_writer Helper function
    Write extracted crosslink sites to a temporary file

    Args:
        output: str, output file name
        chrom: str, chromosome name
        positions: SortedList, crosslink site details
    """
    owriter = output_writer(output, use_tabix=False, preset="bed")
    with owriter(output) as _ow:
        for spos in positions:
            _ow.write(
                f"{chrom}\t{spos[0]}\t{spos[1]}\t{spos[2]}\t{spos[3]}\t{spos[4]}\n"
            )


# def chen_fox_lyndon_factorization(s):
#     n = len(s)  # length of the string
#     i = 0  # start of the string
#     factorization = []  # list to store the factorization
#     while i < n:  # iterate over the string
#         j, k = (
#             i + 1,
#             i,
#         )  # j is the next character index, k is the current character index
#         while (
#             j < n and s[k] <= s[j]
#         ):  # while j is less than length of the string and current character is less than or equal to next character
#             if s[k] < s[j]:  # if current character is less than next character
#                 print(k, s[k], j, s[j], "current character < next character, k=", k)
#                 k = i  # set k to i
#             else:
#                 print(
#                     k,
#                     s[k],
#                     j,
#                     s[j],
#                     "current character == next character, so k=",
#                     k + 1,
#                 )
#                 k += 1  # else increment k by 1
#             j += 1  # increment j by 1
#             print(f"while j<{len(s)} j=", j)
#             # print("i=", i, "j=", j, "k=", k)
#         while i <= k:
#             print("i=", i, "j=", j, "k=", k)
#             factorization.append(s[i : i + j - k])
#             i += j - k
#     assert "".join(factorization) == s
#     return factorization
