import tempfile
from gzip import open as gzopen
from pathlib import Path, PurePath

import pytest


######## general fixtures and functions for tests ########
@pytest.fixture(scope="session")
def tmp_dir():
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmpd:
        yield tmpd


@pytest.fixture
def output_file(tmp_dir):
    """Create a temporary output file path."""
    tmp_out = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)
    tmp_out.close()  # Close it so that the parser can write to it
    yield tmp_out.name
    Path(tmp_out.name).unlink(missing_ok=True)


@pytest.fixture
def output_gz_file(tmp_dir):
    """Create a temporary output file path."""
    tmp_out = tempfile.NamedTemporaryFile(
        dir=tmp_dir, delete=False, suffix=".gz"
    )
    tmp_out.close()  # Close it so that the parser can write to it
    yield tmp_out.name
    Path(tmp_out.name).unlink(missing_ok=True)


@pytest.fixture
def output_parquet_file(tmp_dir):
    """Create a temporary output file path."""
    tmp_out = tempfile.NamedTemporaryFile(
        dir=tmp_dir, delete=False, suffix=".parquet"
    )
    tmp_out.close()  # Close it so that the parser can write to it
    yield tmp_out.name
    Path(tmp_out.name).unlink(missing_ok=True)


def binary_reader(fname, format="plain", mode="rb"):
    """helper function to read files
    Args:
        fname (str): filename to read.
        format (str, optional): File format. Defaults to "plain".
        mode (str, optional): file read mode. Defaults to "rb".
    """
    if format == "plain":
        fn = open
    elif format == "gzip":
        fn = gzopen
    else:
        raise NotImplementedError(
            f"File format {format} not implemented. Use 'plain' or 'gzip'."
        )
    with fn(fname, mode=mode) as fh:
        content = fh.read()
    return content


######## annoation and sliding windows ########


@pytest.fixture
def gff3_basic():
    gff3 = PurePath("tests", "data", "annotation", "basic.gff3.gz")
    return str(gff3)


@pytest.fixture
def basic_annotation():
    """
    expected content of the basic output parser
    also used as input for some tests
    TODO: add test names using this fixture as input
    """
    basic_file = PurePath("tests", "data", "annotation", "basic_output.bed.gz")
    return str(basic_file)


@pytest.fixture
def split_annotation():
    """
    return intron split output
    """
    split_output = PurePath(
        "tests", "data", "annotation", "intron_split.bed.gz"
    )
    return str(split_output)


@pytest.fixture
def expected_bed():
    """
    expected sliding windows output, with tabix index,
    windows=100, step=20
    """
    expected_bed = PurePath(
        "tests", "data", "annotation", "basic_windows_100_20.gz"
    )
    return str(expected_bed)


@pytest.fixture
def expected_index():
    """
    expected sliding windows output, with tabix index,
    windows=100, step=20
    """
    expected_index = PurePath(
        "tests", "data", "annotation", "basic_windows_100_20.gz.tbi"
    )
    return str(expected_index)


@pytest.fixture
def expected_split_bed():
    """
    expected sliding windows output, with tabix index,
    windows=50, step=5
    """
    expected_split_bed = PurePath(
        "tests", "data", "annotation", "split_windows_50_5.gz"
    )
    return str(expected_split_bed)


######## bam extract function fixtures ########


@pytest.fixture
def expected_start_sites():
    """
    expected start sites output
    """
    sites = PurePath(
        "tests", "data", "xlinks", "ENCFF511HSJ_histones.m1_ss_g-1.bed.gz"
    )
    return str(sites)


@pytest.fixture
def expected_middle_sites():
    """
    expected start sites output
    """
    sites = PurePath(
        "tests", "data", "xlinks", "ENCFF511HSJ_histones.m1_sm_g0.bed.gz"
    )
    return str(sites)


@pytest.fixture
def expected_end_sites():
    """
    expected start sites output
    """
    sites = PurePath(
        "tests", "data", "xlinks", "ENCFF511HSJ_histones.m1_se_g0.bed.gz"
    )
    return str(sites)


@pytest.fixture
def expected_deletion_sites():
    """
    expected start sites output
    """
    sites = PurePath(
        "tests", "data", "xlinks", "ENCFF511HSJ_histones.m1_sd_g0.bed.gz"
    )
    return str(sites)


######## count function fixtures ########


@pytest.fixture
def expected_start_counts():
    """
    expected region counts output
    """
    start_counts = PurePath(
        "tests",
        "data",
        "counts",
        "start_regions.parquet",
    )
    return str(start_counts)


@pytest.fixture
def expected_middle_counts():
    """
    expected region counts output
    """
    middle_counts = PurePath(
        "tests",
        "data",
        "counts",
        "middle_regions.parquet",
    )
    return str(middle_counts)


@pytest.fixture
def expected_end_counts():
    """
    expected region counts output
    """
    end_counts = PurePath(
        "tests",
        "data",
        "counts",
        "end_regions.parquet",
    )
    return str(end_counts)


@pytest.fixture
def expected_deletion_counts():
    """
    expected region counts output
    """
    del_counts = PurePath(
        "tests",
        "data",
        "counts",
        "deletion_regions.parquet",
    )
    return str(del_counts)
