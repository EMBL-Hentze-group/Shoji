from pathlib import PurePath

import pytest

from shoji.bam_parser import BamParser
from tests.conftest import binary_reader


@pytest.fixture(scope="module")
def bam():
    """
    return a bam file for testing
    """
    bam_file = PurePath("tests", "data", "xlinks", "ENCFF511HSJ_histones.bam")
    return str(bam_file)


@pytest.fixture(scope="module")
def index():
    """
    return a bam index file for testing
    """
    index_file = PurePath(
        "tests", "data", "xlinks", "ENCFF511HSJ_histones.bam.bai"
    )
    return str(index_file)


@pytest.mark.parametrize(
    "in_bam,index_file,out_file,temp_dir,expected_sites,site,offset,success",
    [
        (
            "bam",
            "index",
            "output_gz_file",
            "tmp_dir",
            "expected_start_sites",
            "s",
            -1,
            True,
        ),
        (
            "bam",
            "index",
            "output_gz_file",
            "tmp_dir",
            "expected_middle_sites",
            "m",
            0,
            True,
        ),
        (
            "bam",
            "index",
            "output_gz_file",
            "tmp_dir",
            "expected_end_sites",
            "e",
            0,
            True,
        ),
        (
            "bam",
            "index",
            "output_gz_file",
            "tmp_dir",
            "expected_deletion_sites",
            "d",
            0,
            True,
        ),
        (
            "bam",
            "index",
            "output_gz_file",
            "tmp_dir",
            "expected_deletion_sites",
            "i",
            0,
            False,  # Invalid site type
        ),
    ],
)
def test_start_sites(
    in_bam,
    index_file,
    out_file,
    temp_dir,
    expected_sites,
    site,
    offset,
    success,
    request,
):
    """
    Parameters:
        in_bam (str): The fixture name for the input BAM file.
        index_file (str): The fixture name for the BAM index file.
        out_file (str): The fixture name for the output file.
        temp_dir (str): The fixture name for the temporary directory.
        expected_sites (str): The fixture name for the expected output sites file.
        site (str): The genomic site to extract.
        offset (int): The offset to apply when extracting the site.
        success (bool): Whether the extraction is expected to succeed.
        request (FixtureRequest): Pytest fixture request object for accessing fixtures.

    """

    with BamParser(
        bam=request.getfixturevalue(in_bam),
        index=request.getfixturevalue(index_file),
        out=request.getfixturevalue(out_file),
        tmp_dir=request.getfixturevalue(temp_dir),
        use_tabix=False,
        cores=1,
    ) as bp:
        bp.extract(
            site=site,
            mate=1,
            offset=offset,
            primary=True,
            min_qual=0,
            min_len=0,
            aln_frac=0.0,
            mismatch_frac=1.0,
            max_len=10000,
            min_aln_len=0,
            max_interval_len=10000,
            ignore_PCR=False,
        )
    # expected bed
    expected_res = binary_reader(
        request.getfixturevalue(expected_sites), "gzip", "rt"
    )
    output_res = binary_reader(request.getfixturevalue(out_file), "gzip", "rt")
    if success:
        assert expected_res == output_res
    else:
        assert expected_res != output_res, (
            f"Expected failure for {site} extraction"
        )
