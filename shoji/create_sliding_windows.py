import multiprocessing as mp
import tempfile
from pathlib import Path
from shutil import rmtree
from typing import Dict, List, Tuple

import pysam
from loguru import logger
from pyarrow.dataset import ParquetFileFragment
from sortedcontainers import SortedList

from . import pyarrow_reader as pr
from .helpers import set_cores, check_tabix
from .output import general_accumulator, output_writer, tabix_accumulator


class SlidingWindows:
    """
    Class to create sliding windows from processed annotation
    """

    def __init__(
        self,
        annotation: str,
        out: str,
        cores: int,
    ) -> None:
        self.annotation = annotation
        self.out = out
        self.cores = set_cores(cores)
        self._temp_dir = tempfile.mkdtemp(dir=Path(out).parent)
        logger.debug(f"Temp. dir {self._temp_dir}")

    def __enter__(self) -> "SlidingWindows":
        return self

    def __exit__(self, except_type, except_val, except_traceback):
        rmtree(self._temp_dir)  # clean up temp dir
        if except_type:
            logger.exception(except_val)

    def generate_sliding_windows(self, step: int, size: int, use_tabix: bool) -> None:
        """
        Generate sliding windows of a given size and step.

        Args:
            step (int): The distance between two adjacent windows.
            size (int): The size each window.
            use_tabix (bool): Whether to use Tabix or not. If True, and the output is gzip, the file will be tabix indexed.
        """
        if check_tabix(self.annotation):
            self._tabix_sliding_windows(
                step=step,
                size=size,
                use_tabix=use_tabix,
            )
        else:
            self._parquet_sliding_windows(
                step=step,
                size=size,
                use_tabix=use_tabix,
            )

    def _tabix_sliding_windows(self, step: int, size: int, use_tabix: bool) -> None:
        """Helper function
        Create sliding windows from an annotation file.

        Args:
            step (int): Size of the window step.
            size (int): Size of the window.
            use_tabix (bool): Whether to use tabix compression
        """
        suffix: str = Path(self.out).suffix
        if use_tabix:
            suffix = ".bed"
        sw_dict: Dict[str, str] = {}
        temp_dir = Path(self._temp_dir)
        with pysam.TabixFile(self.annotation) as _annh:
            # create temp file per chromosome
            for chrom in sorted(_annh.contigs):
                sw_dict[chrom] = str(
                    temp_dir / f"{next(tempfile._get_candidate_names())}{suffix}"  # type: ignore
                )
        with mp.Pool(processes=self.cores) as pool:
            # generate and write position sorted sliding windows per chromosome
            for chrom, temp_file in sw_dict.items():
                pool.apply_async(
                    tabix_sw_worker,
                    args=(self.annotation, temp_file, chrom, step, size),
                )
            pool.close()
            pool.join()

        if use_tabix:
            # gather position sorted data from all chromosomes and write to tabix compressed, indexed file
            tabix_accumulator(sw_dict, self._temp_dir, self.out, "bed")
        else:
            general_accumulator(sw_dict, self.out)

    def _parquet_sliding_windows(self, step: int, size: int, use_tabix: bool) -> None:
        """Helper function
        Generates position-sorted sliding windows from a partitioned Parquet file.

        Args:
            step (int): The step size for the sliding window.
            size (int): The size of the sliding window.
            use_tabix (bool): Whether to write the output as a tabixed compressed, indexed file.
        """
        suffix: str = Path(self.out).suffix
        if use_tabix:
            suffix = ".bed"
        sw_dict: Dict[str, str] = {}
        temp_dir = Path(self._temp_dir)
        with pr.PartionedParquetReader(
            file_name=self.annotation,
            fformat="bed6",
            temp_dir=self._temp_dir,
            cores=self.cores,
        ) as ppq:
            fragments = ppq.get_partitioned_fragments()
            with mp.Pool(processes=self.cores) as pool:
                # generate and write position sorted sliding windows per chromosome
                for chrom, fragment in fragments.items():
                    temp_file = str(
                        temp_dir / f"{next(tempfile._get_candidate_names())}{suffix}"  # type: ignore
                    )
                    sw_dict[chrom] = temp_file
                    pool.apply_async(
                        parquet_sw_worker,
                        args=(fragment, temp_file, chrom, step, size),
                    )
                pool.close()
                pool.join()
        if use_tabix:
            # gather position sorted data from all chromosomes and write to tabix compressed, indexed file
            tabix_accumulator(sw_dict, self._temp_dir, self.out, "bed")
        else:
            general_accumulator(sw_dict, self.out)


def tabix_sw_worker(
    annotation: str, out: str, chrom: str, step: int, size: int
) -> None:
    """tabix_sw_worker tabix worker function
    Function to generate sliding windows from tabix compressed and indexed bed file
    Args:
        annotation: str, Tabix compressed and indexed input file
        out: str, output file name
        chrom: str, chromosome name to extract features
        step: int, step size for sliding windows
        size: int, window size
    """
    logger.debug(
        f"Process id: {mp.current_process().pid}, chromosome: {chrom}, temp. file {out}"
    )
    fc, wc = 0, 0
    wwriter = output_writer(out, use_tabix=False, preset="bed")
    heap = SortedList()
    with pysam.TabixFile(annotation) as _annw, wwriter(out) as _ow:
        for feature in _annw.fetch(chrom, parser=pysam.asBed()):
            fc += 1
            sliding_windows: List[Tuple[int, int]] = _sliding_windows(
                feature.start, feature.end, step, size
            )
            if len(sliding_windows) == 0:
                logger.warning(f"{chrom}:{feature.start}-{feature.end} has no windows!")
                continue
            for i, wi in enumerate(
                _window_indexer(feature.strand, len(sliding_windows))
            ):
                heap.add(
                    (
                        sliding_windows[i][0],
                        sliding_windows[i][1],
                        feature.name + f"W{wi:05}@{wi}",
                        int(feature.score),
                        feature.strand,
                    )
                )
                wc += 1
        for dat in heap:
            # chromsome, start, stop, name, score, strand
            _ow.write(f"{chrom}\t{dat[0]}\t{dat[1]}\t{dat[2]}\t{dat[3]}\t{dat[4]}\n")
    logger.info(f"Finished {chrom}: # features: {fc:,} # windows: {wc:,}")


def parquet_sw_worker(
    fragment: ParquetFileFragment, out: str, chrom: str, step: int, size: int
) -> None:
    """
    Process a Parquet file fragment into sliding windows.

    Args:
        fragment (ParquetFileFragment): Input Parquet file fragment.
        out (str): Output file path.
        chrom (str): Chromosome name.
        step (int): Step size for sliding window.
        size (int): Window size.

    Writes output to the specified file in bed format. Each line
    represents a feature with its start, stop, name, score, and strand information.
    The features are grouped by chromosome and ordered by their start position.
    """
    logger.debug(
        f"Process id: {mp.current_process().pid}, chromosome: {chrom}, temp. file {out}"
    )
    wwriter = output_writer(out, use_tabix=False, preset="bed")
    heap = SortedList()
    fc, wc = 0, 0
    with wwriter(out) as _ow:
        for feature in fragment.to_table().to_pylist():
            fc += 1  # feature count
            sliding_windows: List[Tuple[int, int]] = _sliding_windows(
                feature["chromStart"],
                feature["chromEnd"],
                step,
                size,
            )
            if len(sliding_windows) == 0:
                logger.warning(
                    f"{chrom}:{feature['chromStart']}-{feature['chromEnd']} has no windows!"
                )
                continue
            indexer: range = _window_indexer(feature["strand"], len(sliding_windows))
            for i, wi in enumerate(indexer):
                heap.add(
                    (
                        sliding_windows[i][0],
                        sliding_windows[i][1],
                        feature["name"] + f"W{wi:05}@{wi}",
                        int(feature["score"]),
                        feature["strand"],
                    )
                )
                wc += 1  # window count
        for dat in heap:
            # chromsome, start, stop, name, score, strand
            _ow.write(f"{chrom}\t{dat[0]}\t{dat[1]}\t{dat[2]}\t{dat[3]}\t{dat[4]}\n")
    logger.info(f"Finished {chrom}: # features: {fc:,} # windows: {wc:,}")


def _sliding_windows(
    start: int, end: int, step: int, size: int
) -> List[Tuple[int, int]]:
    """_sliding_windows generate sliding windows
    Given a start and end postion for a feature, a step value and a window size,
    generate a list of sliding windows:
        [(start, start+size), (start+step, start+step+size,)...]
    Args:
        start: int, 0 based start position of the feature
        end: int, end position of the feature
        step: int, step to take from start for the begin position of the next window
        size: int, size of the windw

    Returns:
        List[Tuple[int,int]]
    """
    windows: List[Tuple[int, int]] = []
    wstart: int = start
    wend: int = min(wstart + size, end)
    while wend <= end:
        windows.append((wstart, wend))
        wstart = wstart + step
        wend = wstart + size
    if wstart < end and windows[-1][1] < end:
        windows.append((wstart, end))
    elif windows[-1][1] < end:
        windows.append((windows[-1][1], end))
    return windows


def _window_indexer(strand: str, window_length: int) -> range:
    """_window_indexer Helper function
    Generate window indexing corresponing to window length, strand
    Args:
        strand: str, feature strand
        window_length: int, length of total windows in this feature

    Returns:
        Callable
    """
    if strand == "-":
        return range(window_length, 0, -1)
    else:
        return range(1, window_length + 1, 1)
