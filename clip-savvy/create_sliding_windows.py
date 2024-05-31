import heapq
import logging
from pathlib import Path
from typing import List, Tuple

import pysam

# from .
import pyarrow_reader as pr
from output import output_writer

logger = logging.getLogger(__name__)


class SlidingWindows:
    """
    Class to create sliding windows from processed annotation
    """

    def __init__(
        self, annotation: str, out: str, step: int, size: int, use_tabix: bool
    ) -> None:
        self.annotation = annotation
        self.out = out
        self.step = step
        self.size = size
        self.use_tabix = use_tabix

    def generate_sliding_windows(self) -> None:
        if self._check_tabix():
            _tabix_sliding_windows(
                annotation=self.annotation,
                out=self.out,
                step=self.step,
                size=self.size,
                use_tabix=self.use_tabix,
            )
        else:
            _parquet_sliding_windows(
                annotation=self.annotation,
                out=self.out,
                step=self.step,
                size=self.size,
                use_tabix=self.use_tabix,
            )
        logger.info("Output wrote to %s", self.out)

    def _check_tabix(self) -> bool:
        """_check_tabix check for tabix indices
        Helper function to check for tabx indices (.csi or .tbi)
        Returns:
            bool
        """
        annpath = Path(self.annotation)
        if (annpath.parent / (annpath.name + ".tbi")).exists() or (
            annpath.parent / (annpath.name + ".csi")
        ).exists():
            logger.info("%s is tabix indexed", self.annotation)
            return True
        return False


def _tabix_sliding_windows(
    annotation: str, out: str, step: int, size: int, use_tabix: bool
) -> None:
    """_tabix_sliding_windows generate sliding windows
    Generate sliding windows from tabix indexed bed.gz file

    Args:
        annotation: str, tabix indexed bed.gz file
        out: str, output file name
        step: int, step to take from start position of each feature for next window
        size: int, size of each window
        use_tabix: bool, If True and output has .gz suffix, tabix compress and index output
    """
    owriter = output_writer(out, use_tabix=use_tabix, preset="bed")
    with pysam.TabixFile(annotation) as annh, owriter(out) as oh:  # type: ignore
        chroms: List[str] = sorted(annh.contigs)
        for chrom in chroms:
            # heap to store windows
            heap = []
            fc: int = 0
            for feature in annh.fetch(chrom, parser=pysam.asBed()):
                fc += 1
                for wi, (wstart, wend) in enumerate(
                    _sliding_windows(feature.start, feature.end, step, size)
                ):
                    heapq.heappush(
                        heap,
                        (
                            wstart,
                            wend,
                            feature.name + f"W{wi+1:05}@{wi+1}",
                            int(feature.score),
                            feature.strand,
                        ),
                    )
            logging.info("Chromosome: %s, # features: %s", chrom, f"{fc:,}")
            _heap_windows_writer(heap, chrom, oh)


def _parquet_sliding_windows(
    annotation: str, out: str, step: int, size: int, use_tabix: bool
) -> None:
    tmp_dir = str(Path(out).parent)
    owriter = output_writer(out=out, use_tabix=use_tabix, preset="bed")
    with pr.PartionedParquetReader(
        file_name=annotation, fformat="bed6", temp_dir=tmp_dir
    ) as ppq, owriter(out) as oh:  # type: ignore
        fragments = ppq.get_partitioned_fragments()
        chroms = sorted(fragments.keys())
        for chrom in chroms:
            # heap to store windows
            heap = []
            fc: int = 0
            for feature in fragments[chrom].to_table().to_pylist():
                fc += 1
                for wi, (wstart, wend) in enumerate(
                    _sliding_windows(
                        feature["chromStart"],
                        feature["chromEnd"],
                        step,
                        size,
                    )
                ):
                    heapq.heappush(
                        heap,
                        (
                            wstart,
                            wend,
                            feature["name"] + f"W{wi+1:05}@{wi+1}",
                            int(feature["score"]),
                            feature["strand"],
                        ),
                    )
            logging.info("Chromosome: %s, # features: %s", chrom, f"{fc:,}")
            _heap_windows_writer(heap, chrom, oh)


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


def _heap_windows_writer(heap, chrom, out_handle) -> None:
    """_heap_windows_writer heap writer
    pop items from the given list and write to the given output handle
    Args:
        heap: heap of windows_
        chrom: str, chromosome name
        out_handle: output file writer
    """
    while True:
        try:
            (start, end, name, score, strand) = heapq.heappop(heap)
            out_handle.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")
        except IndexError:
            out_handle.flush()
            break
