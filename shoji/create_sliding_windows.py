import logging
import multiprocessing as mp
import tempfile
from os import cpu_count
from pathlib import Path
from shutil import rmtree
from typing import Dict, List, Tuple

import pysam
from pyarrow.dataset import ParquetFileFragment
from sortedcontainers import SortedList

# import pyarrow_reader as pr
# from output import general_accumulator, output_writer, tabix_accumulator

from . import pyarrow_reader as pr
from .output import general_accumulator, output_writer, tabix_accumulator

logger = logging.getLogger(__name__)


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
        self.cores = cores
        self._temp_dir = tempfile.mkdtemp(dir=Path(out).parent)
        logger.debug("Temp. dir %s", self._temp_dir)

    def __enter__(self) -> "SlidingWindows":
        self._set_cores()
        return self

    def __exit__(self, except_type, except_val, except_traceback):
        rmtree(self._temp_dir)  # clean up temp dir
        if except_type:
            logging.exception(except_val)

    def _set_cores(self) -> None:
        """_set_cores Helper function
        Sanity check self.cores and total number of cores available,
        reset number of available cores if self.cores > total cores
        """
        allcores = cpu_count()
        if (self.cores > allcores) and (allcores > 1):  # type: ignore
            setcores = max(allcores - 1, 1)  # type: ignore
            logger.warning(
                "Give number of cores %i > number of cores detected %i. Setting cores to %i",
                self.cores,
                allcores,
                setcores,
            )
            self.cores = setcores
        elif allcores == 1:
            logger.warning(
                "Available # cores: 1, resetting cores parameter from %i to 1",
                self.cores,
            )
            self.cores = 1
        else:
            logger.info("Using %i cores out of %i...", self.cores, allcores)

    def generate_sliding_windows(self, step: int, size: int, use_tabix: bool) -> None:
        if self._check_tabix():
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

    def _tabix_sliding_windows(self, step: int, size: int, use_tabix: bool) -> None:
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
        suffix: str = Path(self.out).suffix
        if use_tabix:
            suffix = ".bed"
        sw_dict: Dict[str, str] = {}
        temp_dir = Path(self._temp_dir)
        with pr.PartionedParquetReader(
            file_name=self.annotation, fformat="bed6", temp_dir=self._temp_dir
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
    wwriter = output_writer(out, use_tabix=False, preset="bed")
    heap = SortedList()
    with pysam.TabixFile(annotation) as _annw, wwriter(out) as _ow:
        for feature in _annw.fetch(chrom, parser=pysam.asBed()):
            for wi, (wstart, wend) in enumerate(
                _sliding_windows(feature.start, feature.end, step, size)
            ):
                heap.add(
                    (
                        wstart,
                        wend,
                        feature.name + f"W{wi+1:05}@{wi+1}",
                        int(feature.score),
                        feature.strand,
                    )
                )
        for dat in heap:
            # chromsome, start, stop, name, score, strand
            _ow.write(f"{chrom}\t{dat[0]}\t{dat[1]}\t{dat[2]}\t{dat[3]}\t{dat[4]}\n")


def parquet_sw_worker(
    fragment: ParquetFileFragment, out: str, chrom: str, step: int, size: int
) -> None:
    wwriter = output_writer(out, use_tabix=False, preset="bed")
    heap = SortedList()
    with wwriter(out) as _ow:
        for feature in fragment.to_table().to_pylist():
            for wi, (wstart, wend) in enumerate(
                _sliding_windows(
                    feature["chromStart"],
                    feature["chromEnd"],
                    step,
                    size,
                )
            ):
                heap.add(
                    (
                        wstart,
                        wend,
                        feature["name"] + f"W{wi+1:05}@{wi+1}",
                        int(feature["score"]),
                        feature["strand"],
                    )
                )
        for dat in heap:
            # chromsome, start, stop, name, score, strand
            _ow.write(f"{chrom}\t{dat[0]}\t{dat[1]}\t{dat[2]}\t{dat[3]}\t{dat[4]}\n")


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


# import sys


# def main():
#     root = logging.getLogger()
#     root.setLevel(logging.DEBUG)

#     handler = logging.StreamHandler(sys.stdout)
#     handler.setLevel(logging.DEBUG)
#     root.addHandler(handler)
#     tabix_bed = "/workspaces/clip_savvy/test_data/baseline.gencode.v42.tabix.gz"
#     window = 100
#     step = 20
#     tabix_sw = f"/workspaces/clip_savvy/test_data/mp_tabix2plain.{window}_{step}.bed"
#     # sw_tabix = SlidingWindows(annotation=tabix_bed, out=tabix_sw, cores=6)
#     with SlidingWindows(annotation=tabix_bed, out=tabix_sw, cores=6) as swx:
#         swx.generate_sliding_windows(step=step, size=window, use_tabix=False)
#     # sw_tabix.generate_sliding_windows()


# if __name__ == "__main__":
#     main()
