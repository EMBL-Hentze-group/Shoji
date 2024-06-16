import logging
import tempfile
from contextlib import contextmanager
from functools import partial
from pathlib import Path
from shutil import copyfileobj
from typing import Callable, Dict, Generator, Optional, Set
from gzip import open as gzopen
from pysam import tabix_compress, tabix_index

logger = logging.getLogger(__file__)


def output_writer(out: str, use_tabix: bool, preset: str) -> Callable:
    """output_writer output writer
    Output writer for plain text, tabix indexed gzip or compressions supported by xopen

    Args:
        out: str, output file name
        use_tabix: boolean, if True, use tabix to compress and index gzip file
        preset: str, preset to use for tabix indexing

    Raises:
        NotImplementedError: if use_tabix is True and file output is not gz

    Returns:
        Callable, appropriate file writer
    """
    tb_suffix: Set[str] = set([".gz", ".gzip", ".bgz", ".bgzip"])
    gz_suffix: Set[str] = set([".gz", ".gzip"])
    out_suffix = Path(out).suffix.lower()
    if Path(out).exists():
        logging.warning("Re-writing file %s", out)
    if use_tabix and (out_suffix not in tb_suffix):
        raise NotImplementedError(
            f"Cannot use tabix index with {out_suffix} for output file {out}"
        )
    if out_suffix in tb_suffix and use_tabix:
        logger.info("%s format: bgzip, create tabix index: True", out)
        return partial(tabix_writer, preset=preset)
    elif out_suffix in gz_suffix:
        logger.info("%s format: %s", out, out_suffix)
        return partial(gzopen, mode="wt")
    else:
        logger.info("%s format: plain text", out)
        return partial(open, mode="w")


@contextmanager
def tabix_writer(
    file_name: str, preset: str, temp_dir: Optional[str] = None
) -> Generator[tempfile._TemporaryFileWrapper, None, None]:
    """tabix_writer tabix writer
    Write plain text file to a temp. file and use tabix to compress and index file
    Args:
        file_name: str, file name for final output file
        preset: str, preset format to index tabix file
        temp_dir: str, temp. dir to write the temp. file

    Yields:
        Generator[tempfile._TemporaryFileWrapper, None, None]
    """
    tmpfile = tempfile.NamedTemporaryFile(mode="wt", delete=False, dir=temp_dir)
    logger.debug("Tmp file: %s", tmpfile.name)
    try:
        yield tmpfile
        tmpfile.close()
    finally:
        if tmpfile.closed:
            logger.info("Writing bgzip file %s", file_name)
            tabix_compress(tmpfile.name, file_name, force=True)
            logger.info("Indexing bgzip file %s", file_name)
            tabix_index(file_name, force=True, preset=preset)  # type: ignore
        Path(tmpfile.name).unlink(missing_ok=True)


def tabix_accumulator(
    temp_files: Dict[str, str], temp_dir: str, file_name: str, preset: str = "bed"
) -> None:
    """tabix_accumulator gather and write tabix file
    Function to gather position sorted feature/sliding window data from multiple files
    and write to a temp. file, compress the tmp. file to given output file_name and tabix index the file
    Args:
        temp_files: Dict[str,str], key: chromosome name, value: temp. file storing co-ordinate sorted data from the chromosome
        temp_dir: str, temp. directory, directory to write the temp. file
        file_name: str, output file name
        preset: str, preset to index tabix file. Defaults to "bed".
    """
    if Path(file_name).exists():
        logger.warning("Re-writing file %s", file_name)
    with tempfile.NamedTemporaryFile(dir=temp_dir, mode="wb") as tmpsw:
        for chrom in sorted(temp_files.keys()):
            copyfileobj(open(temp_files[chrom], "rb"), tmpsw)
        tmpsw.flush()
        logger.info("Writing bgzip file %s from temp file %s", file_name, tmpsw.name)
        tabix_compress(tmpsw.name, file_name, force=True)
        logger.info("Indexing bgzip file %s", file_name)
        tabix_index(file_name, force=True, preset=preset)  # type: ignore


def general_accumulator(temp_files: Dict[str, str], file_name: str) -> None:
    """general_accumulator general file accumulator
    Function to gather files feature/sliding window data from multiple files of the same format
    and write to a single file of the same format
    Args:
        temp_files: Dict[str,str], key: chromosome name, value: temp. file storing data from the chromosome
        file_name: str, output file name
    """
    if Path(file_name).exists():
        logger.warning("Re-writing file %s", file_name)
    with open(file_name, "wb") as fh:
        for chrom in sorted(temp_files.keys()):
            copyfileobj(open(temp_files[chrom], "rb"), fh)
