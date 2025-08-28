   
![GitHub](https://img.shields.io/github/license/EMBL-Hentze-group/Shoji)

# Shoji
 Shoji is a a flexible command-line toolset for the analysis of iCLIP and eCLIP sequencing data. It is designed as a replacement for [htseq-clip](https://htseq-clip.readthedocs.io/en/latest/), providing streamlined workflows for annotation parsing, crosslink site extraction, counting, and matrix generation.

- **Documentation** : [readthedocs](https://shoji.readthedocs.io/en/latest/)
- **Submit an issue**: [github](https://github.com/EMBL-Hentze-group/Shoji/issues)

## Installation

```
pip install shoji
```

## Features
- Annotation Parsing: Extract and flatten features from GFF3 files to BED format.
- Sliding Window Generation: Create sliding windows over genomic annotations for downstream analysis.
- Crosslink Extraction: Extract crosslink sites from BAM files with flexible options for site type, mate, and filtering.
- Counting: Count crosslink sites per window and output results in Apache Parquet format.
- Matrix Creation: Aggregate counts across samples into R-friendly matrices (CSV or Parquet).
- Tabix Conversion: Convert BED files to bgzipped, tabix-indexed format for efficient querying.

## Major differences to htseq-clip  

- No `--splitExons` flag, Shoji cannot split exons into components  
- New `--split-intron` flag. If an intron overlaps exon from another gene, using this tag will split the intron into non overlapping chunks
- Piping output disabled. Output file names MUST be specified
- `count` function output is only available in Apache parquet format
- `createMatrix` by default do not write duplicate windows. If adjacent overlapping windows have same crosslink counts across all samples, this function now writes only the most 5' (relative to strand) window to output file.  

Developed at [Hentze lab, EMBL Heidelberg](https://www.embl.org/groups/hentze/)

Name inspired from: [Shoji](https://en.wikipedia.org/wiki/Shoji)

