import pyarrow as pa

# schema for count file
count_schema: pa.Schema = pa.schema(
    [
        pa.field("chrom", pa.string(), nullable=False),
        pa.field("begin", pa.uint32(), nullable=False),
        pa.field("end", pa.uint32(), nullable=False),
        pa.field("uniq_id", pa.string(), nullable=False),
        pa.field("gene_id", pa.string(), nullable=False),
        pa.field("gene_name", pa.string(), nullable=False),
        pa.field("gene_type", pa.string(), nullable=False),
        pa.field("feature", pa.string(), nullable=False),
        # this is string to account for split introns
        pa.field("nr_of_region", pa.string(), nullable=False),
        # this is string since nr_of_region is string and to keep the schema consistent
        pa.field("total_region", pa.string(), nullable=False),
        pa.field("window_number", pa.uint16(), nullable=True),
        pa.field("strand", pa.string(), nullable=False),
        pa.field("sample", pa.string(), nullable=False),
        # map of chromosome positions to crosslink counts
        pa.field(
            "pos_counts", pa.map_(pa.uint32(), pa.uint32()), nullable=False
        ),
    ],
    metadata={
        "chrom": "chromosome name",
        "begin": "window start position",
        "end": "window end position",
        "uniq_id": "unique id of this window",
        "gene_id": "gene id",
        "gene_name": "gene name",
        "gene_type": "gene type, eg: protein coding, lncRNA,...",
        "feature": "gene feature, intron or exon",
        "nr_or_region": "number of the current region",
        "total_region": "total number of regions",
        "window_number": "window number, is nullable",
        "strand": "strand info",
        "sample": "sample name",
        "pos_counts": "Map of chromosome positions to crosslink counts",
    },
)

# basic annotation schema for regions
# this is used for the annotation file
region_annotation_schema: pa.Schema = pa.schema(
    [
        pa.field("unique_id", pa.string(), nullable=False),
        pa.field("chromosome", pa.string(), nullable=False),
        pa.field("begin", pa.uint32(), nullable=False),
        pa.field("end", pa.uint32(), nullable=False),
        pa.field("strand", pa.string(), nullable=False),
        pa.field("gene_id", pa.string(), nullable=False),
        pa.field("gene_name", pa.string(), nullable=False),
        pa.field("gene_type", pa.string(), nullable=False),
        pa.field("gene_region", pa.string(), nullable=False),
        pa.field("Nr_of_region", pa.string(), nullable=False),
        pa.field("Total_nr_of_region", pa.string(), nullable=False),
    ],
    metadata={
        "unique_id": "unique id of this window",
        "chromosome": "chromosome name",
        "begin": "window start position",
        "end": "window end position",
        "strand": "strand info",
        "gene_id": "gene id",
        "gene_name": "gene name",
        "gene_type": "gene type, eg: protein coding, lncRNA,...",
        "gene_region": "region type, intron or exon",
        "Nr_of_region": "number of the current region",
        "Total_nr_of_region": "total number of regions",
        "window_number": "Number of this window in the region (optional)",
    },
)


def get_annotation_schema(windowed: bool = False) -> pa.Schema:
    """
    Get the schema for the annotation file.
    """
    if windowed:
        window_annotation_schema = region_annotation_schema.append(
            (pa.field("window_number", pa.uint16(), nullable=False))
        )
        return window_annotation_schema
    else:
        return region_annotation_schema
