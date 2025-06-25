from pathlib import Path

from shoji.create_sliding_windows import SlidingWindows
from tests.conftest import binary_reader


def test_split_intron_windows(
    split_annotation, output_gz_file, expected_split_bed, windows=50, step=5
):
    with SlidingWindows(annotation=split_annotation, out=output_gz_file, cores=1) as sw:
        sw.generate_sliding_windows(step=step, size=windows, use_tabix=False)
    # expected bed
    expected_res = binary_reader(expected_split_bed, "gzip", "rt")
    assert expected_res == binary_reader(output_gz_file, "gzip", "rt")


def test_tabix_out(
    basic_annotation,
    output_gz_file,
    expected_bed,
    expected_index,
    windows=100,
    step=20,
):
    with SlidingWindows(annotation=basic_annotation, out=output_gz_file, cores=1) as sw:
        sw.generate_sliding_windows(step=step, size=windows, use_tabix=True)
    # expected bed
    expected_res = binary_reader(expected_bed)
    # expected index
    expected_indx = binary_reader(expected_index)
    assert expected_res == binary_reader(
        output_gz_file
    ) and expected_indx == binary_reader(output_gz_file + ".tbi")
    # also clean up the index file
    Path(output_gz_file + ".tbi").unlink(missing_ok=True)
