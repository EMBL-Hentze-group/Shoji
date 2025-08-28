from dataclasses import dataclass
from os import cpu_count
from pathlib import Path
from typing import List, NamedTuple

import numpy as np
from loguru import logger

"""_summary_
A collection of general helper functions and modules
"""


def set_cores(ncores: int) -> int:
    """_set_cores Helper function
    Sanity check self.cores and total number of cores available,
    reset number of available cores if self.cores > total cores
    """
    allcores = cpu_count()
    if (ncores > allcores) and (allcores > 1):  # type: ignore
        setcores = max(allcores - 1, 1)  # type: ignore
        logger.warning(
            f"Give number of cores {ncores} > number of cores detected {allcores}. Setting cores to {setcores}"
        )
        return setcores
    if allcores == 1:
        logger.warning(
            f"Available # cores: 1, resetting cores parameter from {ncores} to 1"
        )
        return 1
    else:
        logger.info(f"Using {ncores} cores out of {allcores}...")
        return ncores


def check_tabix(annotation) -> bool:
    """_check_tabix check for tabix indices
    Helper function to check for tabx indices (.csi or .tbi)
    Returns:
        bool
    """
    annpath = Path(annotation)
    if (annpath.parent / (annpath.name + ".tbi")).exists() or (
        annpath.parent / (annpath.name + ".csi")
    ).exists():
        logger.info(f"{annotation} is tabix indexed")
        return True
    return False


def check_bam_index(bam: str) -> str:
    """_check_bam_index check for BAM indices
    Helper function to check for BAM indices (.bai)
    Args:
        bam (str): Path to the BAM file.
    Raises:
        FileNotFoundError: If the BAM index file does not exist.
    Returns:
        str, index file name
    """
    bp: Path = Path(bam)
    bai: Path = bp.with_name(f"{bp.name}.bai")
    csi: Path = bp.with_name(f"{bp.name}.csi")
    if bai.exists():
        logger.debug(f"Found BAM index file: {str(bai)}")
        return str(bai)
    if csi.exists():
        logger.debug(f"Found BAM index file: {str(csi)}")
        return str(csi)
    else:
        # @TODO: improve this!
        raise FileNotFoundError(
            f"BAM index file not found for {bam}. Please run `samtools index {bam}` first!"
        )


class BedFeature(NamedTuple):
    """_summary_
    A named tuple representing a BED feature.
    Args:
        contig: The chromosome or contig where the bed feature is located.
        start: The starting position of the bed feature.
        end: The ending position of the bed feature.
        name: The name of the bed feature.
        score: The score of the bed feature.
        strand: The strand of the bed feature.

    """

    contig: str
    start: int
    end: int
    name: str
    score: int
    strand: int


@dataclass
class Crosslinks:
    """
    A dataclass to store the strand specific crosslink counts for a given chromosome

    Attributes:
       counts (np.ndarray): A NumPy array containing two columns where the first column represents
                            crosslink positions, and the second column contains corresponding counts
                            of those positions. The structure is assumed to be a 2D array with shape-like [N, 2].
       pos (List[int]): A sorted list of unique integers representing the distinct crosslink positions.
    """

    counts: (
        np.ndarray
    )  # first column crosslink pos, second column crosslink counts
    pos: List[int]  # sorted list of crosslink positions
