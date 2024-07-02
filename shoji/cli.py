import click
import logging
from .gff3_parser import GFF3parser
from .create_sliding_windows import SlidingWindows

logger = logging.getLogger(__name__)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group()
@click.version_option()
def run() -> None:
    """
    \b
    Shoji: A flexible toolset for the analysis of iCLIP and eCLIP sequencing data.
    Intended replacement for htseq-clip
    For help on each sub command run:
        shoji <subcommand> -h
    """


@run.command("annotation", context_settings=CONTEXT_SETTINGS)
@click.option(
    "-g",
    "--gff3",
    "gff",
    type=click.Path(exists=True),
    required=True,
    help="GFF3 file to parse (supports .gz files)",
)
@click.option(
    "-o",
    "--out",
    "out",
    type=click.Path(exists=False),
    required=True,
    help="Output bed file name (supports .gz file, tabix indexing)",
)
@click.option(
    "-i",
    "--id",
    "idx_id",
    type=str,
    default="ID",
    show_default=True,
    help="'ID' tag in GFF3 attribute column",
)
@click.option(
    "-p",
    "--parent",
    "parent_id",
    type=str,
    default="Parent",
    show_default=True,
    help="'Parent' tag in GFF3 attribute column",
)
@click.option(
    "-g",
    "--gene_id",
    "gene_id",
    type=str,
    default="gene_id",
    show_default=True,
    help="'Gene id' tag in GFF3 attribute column",
)
@click.option(
    "-n",
    "--gene_name",
    "gene_name",
    type=str,
    default="gene_name",
    show_default=True,
    help="'Gene name' tag in GFF3 attribute column",
)
@click.option(
    "-t",
    "--gene_type",
    "gene_type",
    type=str,
    default="gene_type",
    show_default=True,
    help="'Gene type' tag in GFF3 attribute column",
)
@click.option(
    "-f",
    "--feature",
    "feature",
    type=str,
    default="exon",
    show_default=True,
    help="Gene feature to extract from GFF3 file (from GFF3 3rd column)",
)
@click.option(
    "-x",
    "--gene_like_features",
    "gene_like_features",
    multiple=True,
    default=["tRNA"],
    show_default=True,
    help="'Gene' like featuers to parse from GFF3 (based on GFF3 3rd column). Multiple values can be passed as -x tRNA -x rRNA...",
)
@click.option(
    "--tabix",
    is_flag=True,
    default=False,
    show_default=True,
    help="If the output suffix is '.gz', use this flag to index the output bed file using tabix",
)
@click.option(
    "--split-intron",
    is_flag=True,
    default=False,
    show_default=True,
    help="If an intron overlaps exon of another genes, split this introns into separate bed entries after removing exon overlap",
)
def annotate(
    gff: str,
    out: str,
    idx_id: str,
    parent_id: str,
    gene_id: str,
    gene_name: str,
    gene_type: str,
    feature: str,
    gene_like_features: list,
    tabix: bool,
    split_intron: bool,
) -> None:
    """
    \b
    Parse gff3 file and extract features to bed format
    See `shoji createSlidingWindows -h` for output file use
    """
    parser = GFF3parser(
        gff=gff,
        out=out,
        idx_id=idx_id,
        parent_id=parent_id,
        gene_id=gene_id,
        gene_name=gene_name,
        gene_type=gene_type,
        gene_like_features=set(gene_like_features),
        use_tabix=tabix,
        split_intron=split_intron,
    )
    parser.process(feature_type=feature)


@run.command("createSlidingWindows", context_settings=CONTEXT_SETTINGS)
@click.option(
    "-a",
    "--annotation",
    "annotation",
    type=click.Path(exists=True),
    required=True,
    help="Input annotation file [see 'shoji annotate -h'] to create sliding windows (supports .gz files, tabix indexed files)",
)
@click.option(
    "-o",
    "--out",
    "out",
    type=click.Path(exists=False),
    required=True,
    help="Output sliding windows bed file name (supports .gz file, tabix indexing)",
)
@click.option(
    "-w",
    "--size",
    "size",
    type=int,
    default=100,
    show_default=True,
    help="Size of the sliding window (in bp)",
)
@click.option(
    "-s",
    "--step",
    "step",
    type=int,
    default=20,
    show_default=True,
    help="Step/slide (in bp), from beginning of the previous window to the beginning of the current window",
)
@click.option(
    "--tabix",
    is_flag=True,
    default=False,
    show_default=True,
    help="If the output suffix is '.gz', use this flag to index the output bed file using tabix",
)
@click.option(
    "-c",
    "--cpus",
    "cores",
    type=int,
    default=6,
    show_default=True,
    help="Number of cores to use for parallel processing",
)
def create_sliding_windows(
    annotation: str,
    out: str,
    size: int,
    step: int,
    tabix: bool,
    cores: int,
) -> None:
    """
    \b
    Create sliding windows from flattened annotation
    See `shoji annotate -h` for annotation file creation
    """
    with SlidingWindows(annotation=annotation, out=out, cores=cores) as sw:
        sw.generate_sliding_windows(step=step, size=size, use_tabix=tabix)


if __name__ == "__main__":
    run()
