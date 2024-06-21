import logging
from os import cpu_count

logger = logging.getLogger(__name__)
"""_summary_
A collection of general helper functions
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
            "Give number of cores %i > number of cores detected %i. Setting cores to %i",
            ncores,
            allcores,
            setcores,
        )
        return setcores
    if allcores == 1:
        logger.warning(
            "Available # cores: 1, resetting cores parameter from %i to 1",
            ncores,
        )
        return 1
    else:
        logger.info("Using %i cores out of %i...", ncores, allcores)
        return ncores
