import pytest

from shoji.gff3_parser import GFF3parser
from tests.conftest import binary_reader


def test_gff3parser_init():
    parser = GFF3parser(
        gff="test.gff3",
        out="output.bed",
        gene_like_features={"tRNA", "rRNA"},
        idx_id="ID",
        parent_id="Parent",
        gene_id="gene_id",
        gene_name="gene_name",
        gene_type="gene_type",
        use_tabix=False,
        split_intron=False,
    )
    assert parser.gff == "test.gff3"
    assert parser.out == "output.bed"
    assert parser.gene_like_features == {"trna", "rrna"}  # Lowercased
    assert parser.idx_id == "ID"
    assert parser.parent_id == "Parent"
    assert parser.gene_id == "gene_id"
    assert parser.gene_name == "gene_name"
    assert parser.gene_type == "gene_type"
    assert parser.use_tabix is False
    assert parser.split_intron is False


def test_interval_appender():
    parser = GFF3parser(
        gff="test.gff3",
        out="output.bed",
        gene_like_features=None,
        use_tabix=False,
    )
    assert callable(parser._chrom_appender)


@pytest.mark.parametrize(
    "gff3_test_file,output_test_file,split_intron,expected_output",
    [
        ("gff3_basic", "output_file", False, "basic_annotation"),
        ("gff3_basic", "output_file", True, "split_annotation"),
    ],
)
def test_basic_gff3parser(
    gff3_test_file, output_test_file, split_intron, expected_output, request
):
    """
    Test GFF3parser using file input and file output.
    fixtures with parameterization help: https://miguendes.me/how-to-use-fixtures-as-arguments-in-pytestmarkparametrize
    """
    out = request.getfixturevalue(output_test_file)
    parser = GFF3parser(
        gff=request.getfixturevalue(gff3_test_file),
        out=out,
        gene_like_features=None,
        use_tabix=False,
        split_intron=split_intron,
    )
    parser.process(feature_type="exon")
    # Verify output file exists and contains expected data
    expected_out = binary_reader(
        request.getfixturevalue(expected_output), "gzip", "rt"
    )
    current_out = binary_reader(out, "plain", "r")
    assert expected_out == current_out
