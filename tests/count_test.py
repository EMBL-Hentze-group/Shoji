import tempfile
from pathlib import Path

import pyarrow.parquet as pq
import pytest

from shoji.count import Count
from shoji.schemas import count_schema


def test_count_schema(
    basic_annotation,
    expected_start_sites,
    output_parquet_file,
    tmp_dir,
    sample_name="test",
):
    """_summary_
    annotation: str,
    bed: str,
    out: str,
    sample_name: str,
    cores: int,
    tmp_dir: str,
    """
    with Count(
        annotation=basic_annotation,
        bed=expected_start_sites,
        out=output_parquet_file,
        cores=1,
        tmp_dir=tmp_dir,
        sample_name=sample_name,
    ) as count:
        count.count()
    assert pq.read_schema(output_parquet_file).equals(count_schema)


@pytest.mark.parametrize(
    "annotation,sites,output_pq_file,temp_dir,expected_out,sample_name",
    [
        (
            "expected_bed",
            "expected_start_sites",
            "output_parquet_file",
            "tmp_dir",
            "expected_start_counts",
            "start_regions",
        ),
        (
            "expected_bed",
            "expected_middle_sites",
            "output_parquet_file",
            "tmp_dir",
            "expected_middle_counts",
            "middle_regions",
        ),
        (
            "expected_bed",
            "expected_end_sites",
            "output_parquet_file",
            "tmp_dir",
            "expected_end_counts",
            "end_regions",
        ),
        (
            "expected_bed",
            "expected_deletion_sites",
            "output_parquet_file",
            "tmp_dir",
            "expected_deletion_counts",
            "deletion_regions",
        ),
    ],
)
def test_count_sites(
    annotation,
    sites,
    output_pq_file,
    temp_dir,
    expected_out,
    sample_name,
    request,
):
    """
    Test the count sites functionality.
    """
    with Count(
        annotation=request.getfixturevalue(annotation),
        bed=request.getfixturevalue(sites),
        out=request.getfixturevalue(output_pq_file),
        cores=1,
        tmp_dir=request.getfixturevalue(temp_dir),
        sample_name=sample_name,
    ) as count:
        count.count()
    assert pq.read_table(request.getfixturevalue(output_pq_file)).equals(
        pq.read_table(request.getfixturevalue(expected_out)),
        check_metadata=True,
    )


def test_start_vs_middle_sites(
    expected_bed,
    expected_start_sites,
    expected_middle_sites,
    tmp_dir,
):
    """
    Test that start sites are not equal to middle sites.
    """
    start_out = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)
    start_out.close()
    middle_out = tempfile.NamedTemporaryFile(dir=tmp_dir, delete=False)
    middle_out.close()

    with Count(
        annotation=expected_bed,
        bed=expected_start_sites,
        out=start_out.name,
        cores=1,
        tmp_dir=tmp_dir,
        sample_name="start_sites",
    ) as start_count:
        start_count.count()

    with Count(
        annotation=expected_bed,
        bed=expected_middle_sites,
        out=middle_out.name,
        cores=1,
        tmp_dir=tmp_dir,
        sample_name="middle_sites",
    ) as middle_count:
        middle_count.count()

    assert not pq.read_table(start_out.name).equals(
        pq.read_table(middle_out.name), check_metadata=True
    )
    Path(start_out.name).unlink(missing_ok=True)
    Path(middle_out.name).unlink(missing_ok=True)
