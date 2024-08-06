from email.policy import default
import rich_click as click
from loguru import logger
import sys
from .bam_parser import BamParser
from .create_sliding_windows import SlidingWindows
from .gff3_parser import GFF3parser
from .map_to_id import MapToId
from .count import Count

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

# rich_click styling options
# see https://github.com/ewels/rich-click/tree/main/examples
# config options: https://ewels.github.io/rich-click/documentation/configuration/#configuration-options
click.rich_click.USE_MARKDOWN = True
click.rich_click.STYLE_SWITCH = "bold"
click.rich_click.STYLE_OPTIONS_TABLE_LEADING = 1
click.rich_click.STYLE_OPTIONS_TABLE_BOX = "SIMPLE"
click.rich_click.STYLE_COMMANDS_TABLE_PAD_EDGE = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.COMMAND_GROUPS = {
    "shoji": [
        {
            "name": "Annotation",
            "commands": ["annotation", "createSlidingWindows", "mapToId"],
        },
        {
            "name": "Extraction",
            "commands": ["extract"],
        },
        {
            "name": "Counting",
            "commands": ["count"],
        },
    ]
}


def setup_logger(verbose) -> None:
    logger.remove(0)
    logger.add(
        sys.stderr,
        level=verbose,
        format="|<level>{level: <8}</level>| {message}",
        catch=False,
        enqueue=True,
    )


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option()
def run() -> None:
    """
    Shoji: A flexible toolset for the analysis of iCLIP and eCLIP sequencing data.

    Intended replacement for htseq-clip.

    For help on each sub command run `shoji subcommand -h`
    """


# A common set of options used by multiple sub commands
verbose_option = click.option(
    "-v",
    "--verbose",
    "verbose",
    type=click.Choice(["DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"]),
    default="INFO",
    show_default=True,
    help="Verbose level",
)
cores_option = click.option(
    "-c",
    "--cpus",
    "cores",
    type=int,
    default=6,
    show_default=True,
    help="Number of cores to use for parallel processing",
)
tabix_option = click.option(
    "--tabix",
    is_flag=True,
    default=False,
    show_default=True,
    help="If the output suffix is *.gz*, use this flag to index the output bed file using tabix. It is recommended to use *.gz* suffix and this flag for output files",
)
tmp_option = click.option(
    "-t",
    "--tmp",
    "tmp_dir",
    type=str,
    default=None,
    show_default=True,
    help="Temp. directory to save intermediate outputs. If not provided, creates and uses a temporary directory in --out parent directory",
)


# annotation subcommand
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
    help="Output bed file name (supports .gz compression, tabix indexing)",
)
@click.option(
    "-i",
    "--id",
    "idx_id",
    type=str,
    default="ID",
    show_default=True,
    help="*ID* tag in GFF3 attribute column",
)
@click.option(
    "-p",
    "--parent",
    "parent_id",
    type=str,
    default="Parent",
    show_default=True,
    help="*Parent* tag in GFF3 attribute column",
)
@click.option(
    "-g",
    "--gene_id",
    "gene_id",
    type=str,
    default="gene_id",
    show_default=True,
    help="*Gene id* tag in GFF3 attribute column",
)
@click.option(
    "-n",
    "--gene_name",
    "gene_name",
    type=str,
    default="gene_name",
    show_default=True,
    help="*Gene name* tag in GFF3 attribute column",
)
@click.option(
    "-t",
    "--gene_type",
    "gene_type",
    type=str,
    default="gene_type",
    show_default=True,
    help="*Gene type* tag in GFF3 attribute column",
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
    help="*Gene* like features to parse from GFF3 (based on GFF3 3rd column). Multiple values can be passed as -x tRNA -x rRNA...",
)
@tabix_option
@click.option(
    "--split-intron",
    is_flag=True,
    default=False,
    show_default=True,
    help="If an intron overlaps exon of another genes, split this introns into separate bed entries after removing exon overlap",
)
@verbose_option
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
    verbose: str,
) -> None:
    """
    Parse gff3 file and extract features to bed format.

    See `shoji createSlidingWindows -h` for output file use

    Basic usage: `shoji annotate -a <annotation.gff3> -o <out.bed>`

    Note:

    The default values used for --id, --parent, --gene_id, --gene_name, --gene_type are Gencode GFF3 attribute tags.
    """
    setup_logger(verbose)
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


# createSlidingWindows subcommand
@run.command("createSlidingWindows", context_settings=CONTEXT_SETTINGS)
@click.option(
    "-a",
    "--annotation",
    "annotation",
    type=click.Path(exists=True),
    required=True,
    help="Input annotation file [see `shoji annotate -h`] to create sliding windows (supports .gz files, tabix indexed files)",
)
@click.option(
    "-o",
    "--out",
    "out",
    type=click.Path(exists=False),
    required=True,
    help="Output sliding windows bed file name (supports .gz compression, tabix indexing)",
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
@tabix_option
@cores_option
@verbose_option
def create_sliding_windows(
    annotation: str,
    out: str,
    size: int,
    step: int,
    tabix: bool,
    cores: int,
    verbose: str,
) -> None:
    """
    Create sliding windows from flattened annotation.

    See `shoji annotate -h` for annotation file creation

    Basic usage: `shoji creteSlidingWindows -a <annotation.bed> -o <out_w100_s20.bed> -w 100 -s 20 -c 6`
    """
    setup_logger(verbose)
    with SlidingWindows(annotation=annotation, out=out, cores=cores) as sw:
        sw.generate_sliding_windows(step=step, size=size, use_tabix=tabix)


# mapToId subcommand
@run.command("mapToId", context_settings=CONTEXT_SETTINGS)
@click.option(
    "-a",
    "--annotation",
    "annotation",
    type=click.Path(exists=True),
    required=True,
    help="flattened annotation file from `shoji annotation -h` or sliding window file from `shoji createSlidingWindows -h`",
)
@click.option(
    "-o",
    "--out",
    "out",
    type=click.Path(exists=False),
    required=True,
    help="Output file, supports .gz compression. Region/window annotation mapped to a unique id",
)
@verbose_option
def map_to_id(annotation: str, out: str, verbose: str) -> None:
    """
    map entries in *name* column to unique ids and write in tab separated format
    """
    setup_logger(verbose)
    mid = MapToId(annotation=annotation, out=out)
    mid.map_to_id()


# extract subcommand
@run.command("extract", context_settings=CONTEXT_SETTINGS)
@click.option(
    "-b",
    "--bam",
    "bam",
    type=click.Path(exists=True),
    required=True,
    help="Alignment bam file. Must be co-ordinate sorted and indexed",
)
@click.option(
    "-o",
    "--out",
    "out",
    type=click.Path(exists=False),
    required=True,
    help="Output crosslink sites in BED6 format (supports .gz file, tabix indexing)",
)
@click.option(
    "-e",
    "--mate",
    "mate",
    type=click.Choice(["1", "2"]),
    required=True,
    multiple=False,
    help="for paired end sequencing, select the read/mate to extract the crosslink sites. For single end data, the choice is always 1",
)
@click.option(
    "-s",
    "--site",
    "site",
    type=click.Choice(["s", "m", "e", "i", "d"]),
    required=True,
    multiple=False,
    help="Crosslink site choices, s : start, m : middle, e : end, i : insertion, d : deletion",
)
@click.option(
    "-g",
    "--offset",
    "offset",
    type=int,
    default=0,
    show_default=True,
    help="Number of nucleotides to offset for crosslink sites",
)
@click.option(
    "-q",
    "--qual",
    "min_qual",
    type=int,
    default=0,
    show_default=True,
    help="Minimum alignment quality score",
)
@click.option(
    "-m",
    "--min_len",
    "min_len",
    type=int,
    default=0,
    show_default=True,
    help="Minimum read length",
)
@click.option(
    "-x",
    "--max_len",
    "max_len",
    type=int,
    default=1000,
    show_default=True,
    help="Maximum read length",
)
@click.option(
    "-l",
    "--max_interval_len",
    "max_interval_len",
    type=int,
    default=10000,
    show_default=True,
    help="Maximum read interval length",
)
@click.option(
    "--primary",
    is_flag=True,
    default=False,
    show_default=True,
    help="Flag to use only the  primary alignment positions",
)
@click.option(
    "--ignore_PCR_duplicates",
    "ignore_PCR",
    is_flag=True,
    default=False,
    show_default=True,
    help="Flag to ignore PCR duplicate reads (works only if bam file has PCR duplicate flag set using tools such as samtools markdup)",
)
@tabix_option
@tmp_option
@cores_option
@verbose_option
def extract(
    bam: str,
    out: str,
    mate: str,
    site: str,
    offset: int,
    min_qual: int,
    min_len: int,
    max_len: int,
    max_interval_len: int,
    primary: bool,
    ignore_PCR: bool,
    tabix: bool,
    cores: int,
    tmp_dir: str,
    verbose: str,
) -> None:
    """
    Extract crosslink sites from bam file.

    Crosslinks are extracted based on the read/mate position and the crosslink site choice.

    Crosslinks sites can be:
    - mapped to the start, middle, end of the reads OR
    - either insertion or deletion events in the reads

    Basic usage for eCLIP data:
    `shoji extract -b <bam> -o <out.bed> -e 1 -s s -g -1 -c 6`
    """
    setup_logger(verbose)
    with BamParser(
        bam=bam, out=out, use_tabix=tabix, cores=cores, tmp_dir=tmp_dir
    ) as bp:
        bp.extract(
            site=site,
            mate=int(mate),
            offset=offset,
            min_qual=min_qual,
            min_len=min_len,
            max_len=max_len,
            max_interval_len=max_interval_len,
            primary=primary,
            ignore_PCR=ignore_PCR,
        )


# count subcommand
@run.command("count", context_settings=CONTEXT_SETTINGS)
@click.option(
    "-a",
    "--annotation",
    "annotation",
    type=click.Path(exists=True),
    required=True,
    help="flattened annotation file from `shoji annotation -h` or sliding window file from `shoji createSlidingWindows -h`",
)
@click.option(
    "-i",
    "--input",
    "bed",
    type=click.Path(exists=True),
    required=True,
    help="Extracted crosslink sites in BED format. See `shoji extract -h` for more details.",
)
@click.option(
    "-o",
    "--out",
    "out",
    type=click.Path(exists=False),
    required=True,
    help="Output file, crosslinksite counts per window. Note: This function outputs results only in Apache Parquet format",
)
@click.option(
    "-n",
    "--name",
    "sample_name",
    type=click.STRING,
    default=None,
    help="Sample name to use as a column in the output file. If not provided, the sample name will be inferred from the input file.",
)
@cores_option
@tmp_option
@verbose_option
def count(
    annotation: str,
    bed: str,
    out: str,
    sample_name: str,
    cores: int,
    tmp_dir: str,
    verbose: str,
) -> None:
    """
    Count number of crosslink sites in a window.

    Given an annotation file (see `shoji annotation -h` OR `shoji createSlidingWindows -h`) and
    a bed file (see `shoji extract -h`) with crosslink sites, count the number of crosslink sites in each window.

    Note: This command outputs results only in Apache Parquet format.

    Basic usage for eCLIP data:
    `shoji count -a <annotation> -i <counts> -o <out> -c 6`
    """
    setup_logger(verbose)
    with Count(
        annotation=annotation,
        bed=bed,
        out=out,
        cores=cores,
        tmp_dir=tmp_dir,
        sample_name=sample_name,
    ) as cp:
        cp.count()


if __name__ == "__main__":
    run()
