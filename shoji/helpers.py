from os import cpu_count
from pathlib import Path
from loguru import logger
from typing import NamedTuple

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
