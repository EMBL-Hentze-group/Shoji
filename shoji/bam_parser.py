import multiprocessing as mp
import tempfile
from decimal import ROUND_HALF_UP, Decimal
from functools import partial
from itertools import chain
from pathlib import Path
from shutil import rmtree
from typing import Callable, Dict, List, Optional, Set, Tuple

import numpy as np
import pysam
from loguru import logger
from sortedcontainers import SortedList

from .helpers import set_cores
from .output import general_accumulator, output_writer, tabix_accumulator


class BamParser:
    """Class to parse Bam files.
    Bam file must be co-ordinate sorted and indexed
    """

    def __init__(
        self,
        bam: str,
        index: str,
        out: str,
        use_tabix: bool,
        cores: int,
        tmp_dir: Optional[str] = None,
    ) -> None:
        """__init__ _summary_

        Args:
            bam: bam file to parse. Must be co-ordinate sorted and indexed
            index: BAM index file name
            out: Output file name (bed format)
            use_tabix: boolean, if True, use tabix to index the output file
            cores: int, number of cores to use
            tmp_dir: Tmp. directory to store intermediate outputs. Defaults to None.
        """
        self.bam: str = bam
        self.index: str = index
        self.out: str = out
        self.use_tabix: bool = use_tabix
        # check if the output file suffix is in the list of tabix supported formats
        self._check_suffix()
        self.cores: int = set_cores(cores)
        if tmp_dir is None:
            self._tmp: Path = Path(self.out).parent / next(
                tempfile._get_candidate_names()
            )
        else:
            if not Path(tmp_dir).exists():
                raise FileNotFoundError(f"Directory {tmp_dir} does not exist")
            self._tmp: Path = Path(tmp_dir) / next(
                tempfile._get_candidate_names()
            )
        self._tmp.mkdir(parents=True, exist_ok=True)
        logger.debug(
            f"Using {str(self._tmp)} as temp. directory",
        )
        # list of chromosomes in the bam file
        self._chroms: List[str] = self._get_chromosomes()
        logger.info(f"Found {len(self._chroms)} chromosomes in {self.bam}")
        logger.debug(f"chromosomes: {','.join(self._chroms)}")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        rmtree(self._tmp)  # clean up
        if exc_type:
            logger.exception(traceback)

    def _check_suffix(self) -> None:
        """_check_suffix
        Helper function
        If the suffix of the output file is not in the list of tabix supported formats, raise an error
        Raises:
            NotImplementedError: _description_
        """
        tb_suffix: Set[str] = set([".gz", ".gzip", ".bgz", ".bgzip"])
        out_suffix = Path(self.out).suffix.lower()
        if self.use_tabix and (out_suffix not in tb_suffix):
            raise NotImplementedError(
                f"Cannot use tabix index with '{out_suffix}' suffix for output file {self.out}"
            )

    def _get_chromosomes(self) -> List[str]:
        """_get_chromosomes Helper function
        Check if the input bam file is coordinate sorted and indexed.
        Collect chromosome names from bam header
        """
        with pysam.AlignmentFile(
            self.bam,
            mode="rb",
            check_sq=True,
            require_index=True,
            index_filename=self.index,
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
        min_aln_len: int,
        aln_frac: float,
        mismatch_frac: float,
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
            min_aln_len: int, Minimum aligned read length
            aln_frac: float, Minimum fraction of aligned bases in the read
            mismatch_frac: float, Maximum fraction of mismatches to alignment length, (needs tag "NM" in the bam file)
            max_interval_len: int, Maximum interval length, for paired end reads, splice length otherwise
            primary: bool, Flag to extract only primary alignments
            ignore_PCR: bool, Flag to ignore PCR duplicates (only if bam file has PCR duplicate flag in alignment)
        """
        extract_fn: Callable = self._site_ops(site)
        suffix: str = Path(self.out).suffix
        if self.use_tabix:
            suffix = ".bed"
        temp_dict: Dict[str, str] = {}
        for chrom in self._chroms:
            temp_dict[chrom] = str(
                self._tmp / f"{next(tempfile._get_candidate_names())}{suffix}"  # type: ignore
            )
        with mp.Pool(self.cores) as pool:
            for chrom, temp_file in temp_dict.items():
                pool.apply_async(
                    extract_fn,
                    args=(
                        self.bam,
                        self.index,
                        chrom,
                        temp_file,
                        mate,
                        offset,
                        min_qual,
                        min_len,
                        max_len,
                        min_aln_len,
                        aln_frac,
                        mismatch_frac,
                        max_interval_len,
                        primary,
                        ignore_PCR,
                    ),
                )
            pool.close()
            pool.join()
        if self.use_tabix:
            tabix_accumulator(temp_dict, str(self._tmp), self.out, "bed")
        else:
            general_accumulator(temp_dict, self.out)

    def _site_ops(self, site: str) -> Callable:
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
            raise NotImplementedError(
                f"Site must be one of {sites}, but found {site}"
            )
        if site == "s":
            return partial(extract_single_site, extract_fn=_start)
        elif site == "e":
            return partial(extract_single_site, extract_fn=_end)
        elif site == "m":
            return partial(extract_single_site, extract_fn=_middle)
        elif site == "i":
            return partial(extract_multiple_sites, extract_fn=_insertion)
        else:
            return partial(extract_multiple_sites, extract_fn=_deletion)


def extract_single_site(
    bam: str,
    index: str,
    chrom: str,
    output: str,
    mate: int,
    offset: int,
    min_qual: int,
    min_len: int,
    max_len: int,
    min_aln_len: int,
    aln_frac: float,
    mismatch_frac: float,
    max_interval_len: int,
    primary: bool,
    ignore_PCR: bool,
    extract_fn: Callable,
) -> None:
    """extract_single_site
    Extract one crosslink site/event per read.
    Could be either start site, middle site or end site
    Insertion and deletion sites are not considered here as there can be more than one such event per read
    Args:
        bam: str, BAM file to parse, must be co-ordinate sorted and indexed
        index: str, BAM index file name
        chrom: str, Chromosome name to extract reads
        output: str, tmp. output file name
        site: str, site to extract, could be either s (start), m (middle) or e (end)
        mate: int, mate to extract the crosslink sites from. Must be one of [1,2]
        offset: int, offset start and end sites by "offset" base pairs
        min_qual: int, minimum alignment quality
        min_len: int, minimum read length
        max_len: int, maximum read length
        min_aln_len: int, minimum aligned read length
        aln_frac: float, minimum fraction of aligned bases in the read
        mismatch_frac: float, Maximum fraction of mismatches to alignment length, (needs tag "NM" in the bam file)
        max_interval_len: int, maximum interval length, for paired end reads, splice length otherwise
        primary: bool, flag to extract only primary alignments
        ignore_PCR: bool, flag to ignore PCR duplicates (only if bam file has PCR duplicate flag in alignment)
        extract_fn: Callable, function to extract appropriate crosslink event
    """
    logger.debug(
        f"Process id: {mp.current_process().pid}, chromosome: {chrom}, temp. file: {output}"
    )
    used, discarded = 0, 0
    positions: SortedList = SortedList()
    with pysam.AlignmentFile(bam, mode="rb", index_filename=index) as _bam:
        for aln in _bam.fetch(chrom, multiple_iterators=True):
            if _discard_read(
                aln,
                mate,
                min_qual,
                min_len,
                max_len,
                min_aln_len,
                aln_frac,
                mismatch_frac,
                max_interval_len,
                primary,
                ignore_PCR,
            ):
                discarded += 1
                continue
            start, end = extract_fn(aln, offset)
            if start < 0:
                logger.warning(
                    f"{aln.reference_name}: crosslink site position for read {aln.query_name} is {start}! skipping..."
                )
                discarded += 1
                continue
            strand: str = "-" if aln.is_reverse else "+"
            try:
                yb = aln.get_tag("YB")
            except KeyError:
                yb = 1
            positions.add(
                (start, end, f"{aln.query_name}|{aln.query_length}", yb, strand)
            )
            used += 1
    total = used + discarded
    logger.info(
        f"Finished {chrom}. total reads: {total:,} used reads: {used:,} discarded reads: {discarded:,}"
    )
    _tmp_output_writer(output, chrom, positions)


def extract_multiple_sites(
    bam: str,
    index: str,
    chrom: str,
    output: str,
    mate: int,
    offset: int,
    min_qual: int,
    min_len: int,
    max_len: int,
    min_aln_len: int,
    aln_frac: float,
    mismatch_frac: float,
    max_interval_len: int,
    primary: bool,
    ignore_PCR: bool,
    extract_fn: Callable,
) -> None:
    """extract_multiple_sites
    Extract multiple crosslink sites/events (insertion, deletion) per read.

    Args:
        bam: str, bam file to parse, must be co-ordinate sorted and indexed
        index: str, index file name
        chrom: str, chromosome name
        output: str, output file name
        mate: int, mate to extract the crosslink sites from. Must be one of [1,2]
        offset: int, offset start and end sites by "offset" base pairs, not used here
        min_qual: int, minimum alignment quality
        min_len: int, minimum read length
        max_len: int, maximum read length
        min_aln_len: int, minimum aligned read length
        aln_frac: float, minimum fraction of aligned bases in the read
        mismatch_frac: float, Maximum fraction of mismatches to alignment length, (needs tag "NM" in the bam file)
        max_interval_len: int, maximum interval length, for paired end reads, splice length otherwise
        primary: bool, flag to extract only primary alignments
        ignore_PCR: bool, flag to ignore PCR duplicates (only if bam file has PCR duplicate flag in alignment)
        extract_fn: Callable, extract crosslink sites function
    """
    logger.debug(
        f"Process id: {mp.current_process().pid}, chromosome: {chrom}, temp. file: {output}"
    )
    used, discarded = 0, 0
    positions: SortedList = SortedList()
    with pysam.AlignmentFile(bam, mode="rb", index_filename=index) as _bam:
        for aln in _bam.fetch(chrom, multiple_iterators=True):
            if _discard_read(
                aln,
                mate,
                min_qual,
                min_len,
                max_len,
                min_aln_len,
                aln_frac,
                mismatch_frac,
                max_interval_len,
                primary,
                ignore_PCR,
            ):
                discarded += 1
                continue
            strand: str = "-" if aln.is_reverse else "+"
            try:
                yb = aln.get_tag("YB")
            except KeyError:
                yb = 1
            pos_list: List[Tuple[int, int]] = extract_fn(aln)
            for start, end in pos_list:
                if start < 0:
                    logger.warning(
                        f"{aln.reference_name}: crosslink site position for read {aln.query_name} is {start}! skipping..."
                    )
                    continue
                positions.add(
                    (
                        start,
                        end,
                        f"{aln.query_name}|{aln.query_length}",
                        yb,
                        strand,
                    )
                )
            used += 1
    total = used + discarded
    logger.info(
        f"Finished {chrom}. total reads: {total:,} used reads: {used:,} discarded reads: {discarded:,}"
    )
    _tmp_output_writer(output, chrom, positions)


def _discard_read(
    aln: pysam.AlignedSegment,
    mate: int,
    qual: int,
    min_len: int,
    max_len: int,
    min_aln_len: int,
    aln_frac: float,
    mismatch_frac: float,
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
        min_aln_len: int, Minimum aligned read length
        aln_frac: float, Minimum fraction of aligned bases in the read
        aln_frac: float, Minimum fraction of aligned bases in the read
        mismatch_frac: float, maximum fraction of mismatches to alignment length (needs tag "NM" in the bam file)
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
        or aln.query_alignment_length < min_aln_len
        or (aln.query_alignment_length / aln.query_length) < aln_frac
        or aln.reference_length  # type: ignore
        > max_interval_len  # @TODO: fix this to consider only spliced segments
    ):
        return True
    if (primary and aln.is_secondary) or (ignore_PCR and aln.is_duplicate):
        return True
    if mate == 1 and aln.is_read2:
        return True
    elif mate == 2 and aln.is_read1:
        return True
    if aln.has_tag("NM"):
        # only check for mismatches if the tag "NM" is present
        if (aln.get_tag("NM") / aln.query_alignment_length) > mismatch_frac:  # type: ignore
            return True
        else:
            return False
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
        start: int = aln.reference_end - offset  # type: ignore
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
    start: int = aln.reference_end - offset  # type: ignore
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
    cigars: Set[int] = set([x[0] for x in aln.cigartuples])  # type: ignore
    mid: int = int(
        Decimal(aln.query_alignment_length / 2).quantize(0, ROUND_HALF_UP)
    )

    if (len(cigars - match_ops) == 0) and (aln.is_reverse):
        # only match, soft clip, hard clip, negative strand
        middle: int = aln.reference_end - mid  # type: ignore
        return middle - 1, middle
    elif (len(cigars - match_ops) == 0) and (not aln.is_reverse):
        # only match, soft clip, hard clip, positive strand
        middle: int = aln.reference_start + mid
        return middle - 1, middle
    else:
        # insertion, deletion or splice operations
        pairs: np.ndarray = np.array(aln.get_aligned_pairs())
        try:
            aln_indx: np.ndarray = np.where(
                (pairs[:, 0] != None) & (pairs[:, 1] != None)
            )[0]
        except IndexError:
            return -1, -1
        if aln_indx.shape[0] <= 1:
            return -1, -1
        logger.debug(
            f"{aln.query_name} cigar operations: {','.join([str(c) for c in cigars])}"
        )
        # redefine mid position
        mid = int(Decimal(aln_indx.shape[0] / 2).quantize(0, ROUND_HALF_UP))
        pairs = pairs[aln_indx]
        mid_end = pairs[mid, 1]
        mid_start = pairs[mid - 1, 1]
        if mid_end - mid_start > 1:
            if 3 in cigars:
                logger.warning(
                    f"{aln.reference_name}: {aln.query_name} {mid_start}-{mid_end} middle of the read could be a splice junction!"
                )
            if aln.is_reverse:
                return _get_reverse_mid_pos(pairs, mid)
            else:
                return _get_forward_mid_pos(pairs, mid)
        return mid_start, mid_end


def _get_forward_mid_pos(pairs: np.ndarray, mid: int) -> Tuple[int, int]:
    """_get_forward_mid_pos Helper function

    Args:
        pairs(np.ndarray): A 2D numpy array with aligned pair positions.
                          The first element of each pair is considered for finding reverse positions.
        mid (int): current middle of the aligned read.

    Returns:
        Tuple[int, int]: New middle of the read in BED format
    """
    while mid > 0:
        if pairs[mid, 1] - pairs[mid - 1, 1] == 1:
            return pairs[mid - 1, 1], pairs[mid, 1]
        else:
            mid -= 1
    return -1, -1


def _get_reverse_mid_pos(pairs: np.ndarray, mid: int) -> Tuple[int, int]:
    """_get_reverse_mid_pos Helper function

    Args:
        pairs(np.ndarray): A 2D numpy array with aligned pair positions.
                          The first element of each pair is considered for finding reverse positions.
        mid (int): current middle of the aligned read.

    Returns:
        Tuple[int, int]: New middle of the read in BED format
    """
    while mid < pairs.shape[0]:
        if pairs[mid, 1] - pairs[mid - 1, 1] == 1:
            return pairs[mid - 1, 1], pairs[mid, 1]
        else:
            mid += 1
    return -1, -1


def _insertion(aln: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    """_insertion Helper function
    Get insertion points from the read
    Args:
        aln: pysam.AlignedSegment, Aligned read

    Returns:
        List[Tuple[int, int]], List of insertion points (start, end)
    """
    cigars: Set[int] = set([x[0] for x in aln.cigartuples])  # type: ignore
    if 1 not in cigars:
        # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
        # No insertion found
        return []
    return _insertion_deletion_points(
        aln.cigartuples,
        aln.get_aligned_pairs(),
        1,  # type: ignore
    )


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
    return _insertion_deletion_points(
        aln.cigartuples,
        aln.get_aligned_pairs(),
        2,  # type: ignore
    )


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
        ops_locations.append(
            (cigar_ref[start_index, 1], cigar_ref[end_index, 1])
        )  # type: ignore
    return ops_locations


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
