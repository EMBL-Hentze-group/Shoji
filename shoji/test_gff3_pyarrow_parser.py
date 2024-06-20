import pytest
from gff3_pyarrow_parser import GFF3parser


def test_GFF3parser():
    # Create an instance of GFF3parser
    gff_parser = GFF3parser(
        gff="input.gff",
        out="output.bed",
        gene_like_features={"tRNA"},
        idx_id="ID",
        parent_id="Parent",
        gene_id="gene_id",
        gene_name="gene_name",
        gene_type="gene_type",
        use_tabix=True,
        split_intron=False,
    )

    # Test the initialization of the GFF3parser instance
    assert gff_parser.gff == "input.gff"
    assert gff_parser.out == "output.bed"
    assert gff_parser.gene_like_features == {"tRNA"}
    assert gff_parser.idx_id == "ID"
    assert gff_parser.parent_id == "Parent"
    assert gff_parser.gene_id == "gene_id"
    assert gff_parser.gene_name == "gene_name"
    assert gff_parser.gene_type == "gene_type"
    assert gff_parser.use_tabix == True
    assert gff_parser.split_intron == False

    # Test the process() method
    with pytest.raises(RuntimeError):
        gff_parser.process(feature_type="exon")

    # Test the _read_gff3() method
    gff3_pa = gff_parser._read_gff3()
    assert isinstance(gff3_pa, pa.Table)

    # Test the _gene_feature_dependancy() method
    gff_parser._gene_feature_dependancy()
    assert isinstance(gff_parser._gene_graph, nx.DiGraph)

    # Test the _parse_gene_features() method
    gff_parser._parse_gene_features(feature_type="exon")
    assert isinstance(gff_parser._gene_feature_map, dict)

    # Test the _fill_strand_intervals() method
    gff_parser._fill_strand_intervals()
    assert isinstance(gff_parser._strand_intervals, dict)

    # Test the _sort_n_filter_intervals() method
    gff_parser._sort_n_filter_intervals()
    assert isinstance(gff_parser._heap, list)

    # Test the _write() method
    with open("output.bed", "w") as fh:
        gff_parser._write(fh, chrom="chr1")
    # Add assertions for the output file if needed
