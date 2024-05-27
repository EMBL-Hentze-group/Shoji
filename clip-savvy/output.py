import logging
import tempfile
from contextlib import contextmanager
from functools import partial
from pathlib import Path
from typing import Callable, Generator, Set

from pysam import tabix_compress, tabix_index
from xopen import xopen

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
    tabix_suffix: Set[str] = set([".gz", ".gzip", ".bgz", ".bgzip"])
    xopen_suffix: Set[str] = set([".gz", ".bz2", ".xz"])
    out_suffix = Path(out).suffix.lower()
    if Path(out).exists():
        logging.warning("Re-writing file %s", out)
    if use_tabix and (out_suffix not in tabix_suffix):
        raise NotImplementedError(
            f"Cannot use tabix index with {out_suffix} for output file {out}"
        )
    elif out_suffix in tabix_suffix and use_tabix:
        logger.info("%s format: bgzip, create tabix index: True", out)
        return partial(tabix_writer, preset=preset)
    elif out_suffix in xopen_suffix:
        logger.info("%s format: %s", out, out_suffix)
        return partial(xopen, mode="wt", compresslevel=None)
    else:
        logger.info("%s format: plain text", out)
        return partial(open, mode="w")


@contextmanager
def tabix_writer(
    file_name, preset
) -> Generator[tempfile._TemporaryFileWrapper, None, None]:
    """tabix_writer tabix writer
    Write plain text file to a temp. file and use tabix to compress and index file
    Args:
        file_name: str, file name for final output file
        preset: str, preset format to index tabix file

    Yields:
        Generator[tempfile._TemporaryFileWrapper, None, None]
    """
    tmpfile = tempfile.NamedTemporaryFile(mode="wt", delete=False)
    logging.debug("Tmp file: %s", tmpfile.name)
    try:
        yield tmpfile
        tmpfile.close()
    finally:
        if tmpfile.closed:
            logging.info("Writing bgzip file %s", file_name)
            tabix_compress(tmpfile.name, file_name, force=True)
            logging.info("Indexing bgzip file %s", file_name)
            tabix_index(file_name, force=True, preset=preset)
        Path(tmpfile.name).unlink(missing_ok=True)
