import os
import tempfile

import polars as pl
import pytest

from shoji.create_matrix import CreateMatrix


def test_bad_schema(count_dir, tmp_dir):
    annotation_file = tempfile.NamedTemporaryFile(
        dir=tmp_dir, delete=False, suffix=".parquet"
    )
    annotation_file.close()
    max_file = tempfile.NamedTemporaryFile(
        dir=tmp_dir, delete=False, suffix=".parquet"
    )
    max_file.close()
    count_file = tempfile.NamedTemporaryFile(
        dir=tmp_dir, delete=False, suffix=".parquet"
    )
    count_file.close()

    with (
        pytest.raises(RuntimeError) as r_error,
        CreateMatrix(
            in_dir=count_dir,
            annotation=annotation_file.name,
            out=count_file.name,
            max_out=max_file.name,
            prefix="bad_schema",
            suffix=".parquet",
            format="parquet",
            tmp_dir=tmp_dir,
            cores=2,
        ) as cm,
    ):
        cm.create_matrices(allow_duplicates=False)


@pytest.mark.parametrize(
    "input_dir,expected_annotation,expected_count,expected_max,temp_dir,prefix,suffix,allow_duplicates",
    [
        (
            "count_dir",
            "start_annotation",
            "start_count",
            "start_max",
            "tmp_dir",
            "start",
            ".parquet",
            False,
        ),
        (
            "count_dir",
            "start_annotation_dup",
            "start_count_dup",
            "start_max_dup",
            "tmp_dir",
            "start",
            ".parquet",
            True,
        ),
        (
            "count_dir",
            "middle_annotation",
            "middle_count",
            "middle_max",
            "tmp_dir",
            "middle",
            ".parquet",
            False,
        ),
    ],
)
def test_create_matrix(
    input_dir,
    expected_annotation,
    expected_count,
    expected_max,
    temp_dir,
    prefix,
    suffix,
    allow_duplicates,
    request,
):
    out_dir = request.getfixturevalue(temp_dir)
    annotation_file = tempfile.NamedTemporaryFile(
        dir=out_dir, delete=False, suffix=".parquet"
    )
    annotation_file.close()
    max_file = tempfile.NamedTemporaryFile(
        dir=out_dir, delete=False, suffix=".parquet"
    )
    max_file.close()
    count_file = tempfile.NamedTemporaryFile(
        dir=out_dir, delete=False, suffix=".parquet"
    )
    count_file.close()

    os.environ.setdefault("POLARS_MAX_THREADS", str(1))
    with CreateMatrix(
        in_dir=request.getfixturevalue(input_dir),
        annotation=annotation_file.name,
        out=count_file.name,
        max_out=max_file.name,
        prefix=prefix,
        suffix=suffix,
        format="parquet",
        tmp_dir=str(out_dir),
        cores=2,
    ) as cm:
        cm.create_matrices(allow_duplicates=allow_duplicates)

    # annotations
    test_annot = pl.read_parquet(annotation_file.name).sort("unique_id")
    expected_annot = pl.read_parquet(
        request.getfixturevalue(expected_annotation)
    ).sort("unique_id")
    assert expected_annot.equals(test_annot)
    # counts
    test_count = pl.read_parquet(count_file.name).sort("unique_id")
    expected_count = pl.read_parquet(
        request.getfixturevalue(expected_count)
    ).sort("unique_id")
    assert expected_count.equals(test_count)
    # max counts
    test_max = pl.read_parquet(max_file.name).sort("unique_id")
    expected_max = pl.read_parquet(request.getfixturevalue(expected_max)).sort(
        "unique_id"
    )
    assert expected_max.equals(test_max)
