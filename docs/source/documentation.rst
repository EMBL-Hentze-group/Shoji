.. raw:: html

    <style> .c1 {color:#3282b8; font-weight:bold} </style>

.. role:: c1

documentation
=============

After the successful installation, use

.. code-block:: sh

    $ shoji --help
to see the list of available commands and their options.

In Shoji, main functions are organized into four categories:

.. list-table::

    * - Annotation
      - | ``annotation``
        | ``createSlidingWindows``
      - | Parse gff3 file and extract features to bed format
        | Create sliding windows from flattened annotation.
    * - Extraction
      - ``extract``
      - Extract crosslink sites from alignment file.
    * - Counting
      - | ``count``
        | ``createMatrix``
      - | Count number of crosslink sites in a window.
        | create R friendly output matrices
    * - Helpers
      - ``toTabix``
      - Convert BED file to bgzipped, tabix indexed BED file

.. _AnnotationOverview:

Annotation
***********

.. _annotation:

:c1:`annotation`
---------------
Parse gff3 file and extract features to bed format.

**Arguments**

* ``--gff3/-a``   GFF3 file to parse (supports .gz files)
* ``--out/-o``   Output bed file name (supports .gz compression, tabix indexing)
* ``--id/-i``   ID tag in GFF3 attribute column
* ``--parent/-p``   Parent tag in GFF3 attribute column
* ``--gene_id/-g``   Gene id tag in GFF3 attribute column
* ``--gene_name/-n``   Gene name tag in GFF3 attribute column
* ``--gene_type/-t``   Gene type tag in GFF3 attribute column
* ``--feature/-f``   Gene feature to extract from GFF3 file (from GFF3 3rd column)
* ``--gene_like_features/-x``   Gene like features to parse from GFF3 (based on GFF3 3rd column). Multiple values can be passed as -x tRNA -x rRNA...
* ``--tabix`` If the output suffix is .gz, use this flag to index the output bed file using tabix. It is recommended to use .gz suffix and this flag for output
* ``--split-intron`` If an intron overlaps exon of another genes, split this introns into separate bed entries after removing exon overlap

.. Note::
    All default values are based on Gencode GFF3 file format.


**Usage:**

.. code-block:: sh

    $ shoji annotation --help

.. _createSlidingWindows:

:c1:`createSlidingWindows`
--------------------------

Create sliding windows from flattened annotation

**Arguments**

* ``--annotation/-a``   Input annotation file (see ``shoji annotate -h``) to create sliding windows (supports .gz files, tabix indexed files)
* ``--size/-w``   Size of the sliding window (in bp)
* ``--step/-s``   Step/slide (in bp), from beginning of the previous window to the beginning of the current window
* ``--tabix`` If the output suffix is .gz, use this flag to index the output bed file using tabix. It is recommended to use .gz suffix and this flag for output files
* ``--cpus/-c``   Number of cores to use for parallel processing

**Usage:**

.. code-block:: sh

    $ shoji createSlidingWindows --help

Extraction
***********

.. _extract:

:c1:`extract`
---------------
Extract crosslink sites from bam file.

**Arguments**

* ``--bam/-b``   Alignment bam file. Must be co-ordinate sorted and indexed
* ``--out/-o``   Output crosslink sites in BED6 format (supports .gz file, tabix indexing)
* ``--mate/-e``   for paired end sequencing, select the read/mate to extract the crosslink sites. For single end data, the choice is always 1
* ``--site/-s``   Crosslink site choices, s : start, m : middle, e : end, i : insertion, d : deletion
* ``--offset/-g``   Number of nucleotides to offset for crosslink sites
* ``--qual/-q``   Minimum alignment quality score
* ``--min_len/-m``   Minimum read length
* ``--max_len/-x``   Maximum read length
* ``--min_aln_len/-a``   Minimum aligned read length
* ``--aln_frac/-f``   Minimum fraction of aligned read length to total read length for crosslink site extraction. If set to 0, all reads are considered 
* ``--mismatch_frac/-y``   Maximum fraction of mismatches allowed in the read, as a fraction of aligned length. If set to 1.0, all reads are considered
* ``--max_interval_len/-l``   Maximum read interval length
* ``--primary`` Flag to use only the  primary alignment position
* ``--ignore_PCR_duplicates`` Flag to ignore PCR duplicate reads (works only if bam file has PCR duplicate flag set using tools such as samtools markdup)
* ``--tabix`` If the output suffix is .gz, use this flag to index the output bed file using tabix
* ``--tmp/-t``   Temp. directory to save intermediate outputs. If not provided, creates and uses a temporary directory in * ``--out`` parent directory
* ``--cpus/-c``   Number of cores to use for parallel processing

**Usage:**

.. code-block:: sh

    $ shoji extract --help

Counting
***********

.. _count:

:c1:`count`
---------------
Count number of crosslink sites in a window.

**Arguments**

* ``--annotation/-a``   flattened annotation file from shoji annotation -h or sliding window file from ``shoji createSlidingWindows -h``
* ``--input/-i``   Extracted crosslink sites in BED format. See ``shoji extract -h`` for more details
* ``--out/-o``  Output file, crosslinksite counts per window. Note: This function outputs results only in Apache Parquet format
* ``--name/-n``  Sample name to use as a column in the output file. If not provided, the sample name will be inferred from the input file
* ``--cpus/-c``  Number of cores to use for parallel processing          
* ``--tmp/-t``   Temp. directory to save intermediate outputs. If not provided, creates and uses a temporary directory in  ``--out`` parent directory

**Usage:**

.. code-block:: sh

    $ shoji count --help

.. _createMatrix:

:c1:`createMatrix`
------------------
create R friendly output matrices.

**Arguments**

* ``--input_dir/-i``   Input directory containing the output of shoji count, see ``shoji count -h`` for details
* ``--prefix/-p``   Prefix to filter count files in in_dir
* ``--suffix/-s``   Suffix to filter count files in in_dir. Either ``--prefix`` or ``--suffix`` MUST be provided
* ``--format/-f ``  Output formats: ``csv`` or ``parquet``. Default: csv. Csv format also supports gzipped output
* ``--annotation/-a``   Output filename for trimmed annotations
* ``--output/-o``   Output filename for aggregated crosslink count per window matrix
* ``--max/-m``   Optional output. Output filename for max. crosslink site per window matrix
* ``--allow_duplicates``       Default behavior: If adjacent overlapping windows have same crosslink counts across all samples, write only the most 5' window to output file. Use this flag disable this feature and to write all windows
* ``--cpus/-c``  Number of cores to use for parallel processing 

**Usage:**

.. code-block:: sh

    $ shoji createMatrix --help


Helpers
***********

.. _toTabix:

:c1:`toTabix`
---------------

Convert BED file to bgzipped, tabix indexed BED file  

**Arguments:**

*  ``--bed/-b``   Input BED file (bed6 format, supports .gz files)
*  ``--output/-o``  Output filename for bgzipped tabix indexed bed file
*  ``--cpus/-c``   Number of cores to use for parallel processing
*  ``--tmp/-t``   Temp. directory to save intermediate outputs. If not provided, creates and uses a temporary directory in ``--out`` parent directory 

**Usage:**

.. code-block:: sh

    $ shoji toTabix --help