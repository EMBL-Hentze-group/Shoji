.. raw:: html

    <style> .c1 {color:#3282b8; font-weight:bold} </style>

.. role:: c1


overview
=========

**Shoji**

This package is designed primarily to do the following operations:

:ref:`Annotation <AnnotationOverview>`
**********************************************

A suite of functions to process and flatten genome annotation file. 

.. _annotation:
:c1:`annotation`
----------------

:ref:`annotation function <annotation>` takes as input a GFF formatted genome annotation file and converts the annotations from GFF format to bed format.
For an example, this function converts the following GFF annotation

.. _GFFTable:

.. list-table::
   
  * - chr1
    - HAVANA
    - exon
    - 17233
    - 17368
    - .
    - \-
    - .
    - ID=exon:ENST00000488147.1:6;Parent=ENST00000488147.1;gene_id=ENSG00000227232.5;transcript_id=ENST00000488147.1;gene_type=unprocessed_pseudogene;gene_name=WASH7P;transcript_type=unprocessed_pseudogene;transcript_name=WASH7P-201;exon_number=6;exon_id=ENSE00003502542.1;level=2;transcript_support_level=NA;hgnc_id=HGNC:38034;ont=PGO:0000005;tag=basic,Ensembl_canonical;havana_gene=OTTHUMG00000000958.1;havana_transcript=OTTHUMT00000002839.1


and converts this entry into the following BED6 format

.. _BEDTable:

.. list-table::
    :header-rows: 1
    
    * - chromosome
      - start
      - end
      - name
      - score
      - strand
    * - chr1
      - 17232
      - 17368
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@exon@6/11@ENSG00000227232.5:exon0006
      - 0
      - \-
Various attributes in the name column in this BED entry is seperated by ``@`` and the
order is given below

.. _AttibTable:

.. list-table::
    :widths: 3,10
    :header-rows: 1
    

    * - atrribute
      - attribute description 
    * - ENSG00000227232.5
      - gene id
    * - WASH7P
      - gene name
    * - unprocessed_pseudogene
      - gene type
    * - exon
      - gene feature (exon, intron, CDS,...)
    * - 6/11
      - 6th exon out of a total of 11 exons of this gene
    * - ENSG00000227232.5:exon0006
      - unique id, merging gene id feature and feature number

In `Gencode V42`_ genome build, the above entry is the 6th exon of the gene WASH7P, preceded by an intron from ``chr1:17368-17605(-)``, 
but chromosomal location ``chr1:17368-17436(-)`` within this intron contains a microRNA MIR6859-1. 
In Shoji and `htseq-clip`_, the default behavior is to process the intron and microRNA with the intron as separate features:

.. list-table::
    :header-rows: 1
    
    * - chromosome
      - start
      - end
      - name
      - score
      - strand
    * - chr1
      - 17232
      - 17368
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@exon@6/11@ENSG00000227232.5:exon0006
      - 0
      - \-
    * - chr1
      - 17368
      - 17436
      - ENSG00000278267.1@MIR6859-1@miRNA@exon@1/1@ENSG00000278267.1:exon0001
      - 0
      - \-
    * - chr1
      - 17368
      - 17605
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@intron@5/10@ENSG00000227232.5:intron0005
      - 0
      - \-
    
But with ``--split-intron`` flag, Shoji will split the intron into two non-overlapping chunks, and removes the chunk that overlaps with the microRNA. The output will be:

.. list-table::
    :header-rows: 1
    
    * - chromosome
      - start
      - end
      - name
      - score
      - strand
    * - chr1
      - 17232
      - 17368
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@exon@6/11@ENSG00000227232.5:exon0006
      - 0
      - \-
    * - chr1
      - 17368
      - 17436
      - ENSG00000278267.1@MIR6859-1@miRNA@exon@1/1@ENSG00000278267.1:exon0001
      - 0
      - \-
    * - chr1
      - 17436
      - 17605
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@intron@5-1/10@ENSG00000227232.5:intron0005-1
      - 0
      - \-

Notice that the intron now spand ``chr1:17436-17605(-)``  instead of spanning over the microRNA, ie  ``chr1:17368-17605(-)``.

.. Note:: The output supports `tabix`_ indexing using the flag ``--tabix`` and  is still compatible with `htseq-clip`_ commands

.. _createSlidingWindows:
:c1:`createSlidingWindows`
--------------------------

:ref:`createSlidingWindows function <createSlidingWindows>` takes as input a flattened annotation BED file
created by the annotation function and splits each individual BED entries into overlapping windows. 
``--size`` parameter controls the size of each window and ``--step`` controls the overlap 
of each neighboring windows from the same feature

Continuing with the example entry above, the first 5 sliding windows generated from the
:ref:`BED6 flattened entry <BEDTable>` are given below:

.. _SWTable:

.. list-table::
    :header-rows: 1
        
    * - chromosome
      - start
      - end
      - name
      - score
      - strand
    * - chr1
      - 17232
      - 17282
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@exon@6/11@ENSG00000227232.5:exon0006W00019@19
      - 0
      - \-
    * - chr1
      - 17237
      - 17287
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@exon@6/11@ENSG00000227232.5:exon0006W00018@18
      - 0
      - \-
    * - chr1
      - 17242
      - 17292
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@exon@6/11@ENSG00000227232.5:exon0006W00017@17
      - 0
      - \-
    * - chr1
      - 17247
      - 17297
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@exon@6/11@ENSG00000227232.5:exon0006W00016@16
      - 0
      - \-
    * - chr1
      - 17252
      - 17302
      - ENSG00000227232.5@WASH7P@unprocessed_pseudogene@exon@6/11@ENSG00000227232.5:exon0006W00015@15
      - 0
      - \-

Each sliding window listed here is 50bp long, as default value for ``--size`` argument is ``50``  and the difference between
start positions of each is 20bp, as the default value for ``--step`` argument is ``20`` 

Following the convention in :ref:`flattened annotation <BEDTable>` the attributes in sliding windows name column are also seperated by ``@`` 
and the first 5 attributes in the name column here are exactly the same as that of :ref:`flattened annotation name column <AttibTable>`
An example is given below

.. _SWAttibTable:

.. list-table::
    :header-rows: 1

    * - atrribute
      - attribute description
      - Found in :ref:`flattend name attribute <AttibTable>`
    * - ENSG00000227232.5
      - gene id
      - Yes
    * - WASH7P
      - gene name
      - Yes
    * - unprocessed_pseudogene
      - gene type
      - Yes
    * - exon
      - gene feature (exon, intron, CDS,...)
      - Yes
    * - 6/11
      - 6th exon out of a total 11 exons in this gene
      - Yes
    * - ENSG00000227232.5:exon0006W00019
      - unique id, merging gene id feature, feature number and window number (W : window)
      - No
    * - 19
      - 19th window of this feature 
      - No
 
.. Note:: There will be zero overlap between neighboring windows from two separate gene features

:ref:`Extraction <ExtractionOverview>`
***************************************************

Extract and process crosslink sites from alignment file.

:c1:`extract`
--------------

:ref:`extract function <extract>` takes as input an alignment file (.bam) and extracts and 
writes either start, insertion, deletion, middle or end site into a BED6 formatted file.
The argument ``--site``  determines crosslink site choice.

Given below is an example paired end sequence and start, middle and end positions extracted from the second mate of this fragment

.. _AlignTable1:

.. list-table::

  * - TTATTACAGC\:K00180\:131\:H7J3YBBXX\:3:2123:15057:19918
    - 99
    - chr1
    - 17252
    - 255
    - 33M
    - \=
    - 17285
    - 41
    - TTTTAAAGGCTGAGTCCTCTGAGAATTTATTAC
    - JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
    - NH:i:1
    - HI:i:1
    - AS:i:60
    - nM:i:5
    - NM:i:4
    - MD:Z:0C0A0G0G29
    - jM:B:c,-1
    - jI:B:i,-1
    - RG:Z:foo
  * - TTATTACAGC\:K00180\:131\:H7J3YBBXX\:3:2123:15057:19918
    - 147
    - chr1
    - 17285
    - 255
    - 38M
    - \=
    - 17252
    - \-41
    - TAAAGGCTGAGTCCTCTGAGAATTTATTACTACGGATC
    - JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
    - NH:i:1
    - HI:i:1
    - AS:i:60
    - nM:i:5
    - NM:i:1
    - MD:Z:0G37
    - jM:B:c,-1
    - jI:B:i,-1
    - RG:Z:foo


**start site**

.. _StartTable:

.. list-table::

  * - chromosome
    - start
    - end
    - name
    - score
    - strand
  * - chr1
    - 17322
    - 17323
    - TTATTACAGC\:K00180\:131\:H7J3YBBXX\:3:2123:15057:19918|38
    - 1
    - \-

**middle site**

.. _MiddleTable:

.. list-table::

  * - chromosome
    - start
    - end
    - name
    - score
    - strand
  * - chr1
    - 17304
    - 17305
    - TTATTACAGC\:K00180\:131\:H7J3YBBXX\:3:2123:15057:19918|38
    - 1
    - \-

**end site**

.. _EndTable:

.. list-table::

  * - chromosome
    - start
    - end
    - name
    - score
    - strand
  * - chr1
    - 17285
    - 17286
    - TTATTACAGC\:K00180\:131\:H7J3YBBXX\:3:2123:15057:19918|38
    - 1
    - \-

.. Note:: In a paired end alignment file, argument ``--mate`` is used to choose the read/mate from which crosslink sites are extracted. The sequencing protocol used to generate the file determines whether the crosslink site is located on the first mate or the second mate. Please consult your sequencing protocol to decide which mate to use.

:ref:`Counting <CountOverview>`
**********************************************

Calculate the number of extracted crosslink sites per given gene annotation feature.

.. _countoverview:
:c1:`count`
------------------

:ref:`count function <count>` takes as input either a flattened annotation file generated by annotation function or a sliding windows
file generated by createSlidingWindows function and a crosslink sites file generated by extract function and for each entry/window in the
annotation/sliding windows file count the number of crosslink sites in the region.

In Shoji, the output is written in `Apache parquet`_ file format, with the following schema:

.. _Count-Table-Schema:

.. list-table::

  * - field (column) name
    - type
    - nullable
    - description
  * - chrom
    - string
    - False
    - chromosome name
  * - begin
    - uint32
    - False
    - window start position
  * - end
    - uint32
    - False
    - window end position
  * - uniq_id
    - string
    - False
    - unique id of this window
  * - gene_id
    - string
    - False
    - gene id
  * - gene_name
    - string
    - False
    - gene name
  * - gene_type
    - string
    - False
    - gene type, eg: protein coding, lncRNA,...
  * - feature
    - string
    - False
    - gene feature, intron or exon
  * - nr_of_region
    - string
    - False
    - number of the current region
  * - total_region
    - string
    - False
    - total number of regions
  * - window_number
    - uint16
    - True
    - window number
  * - strand
    - string
    - False
    - strand info
  * - sample
    - string
    - False
    - sample name
  * - pos_counts
    - map_(uint32, uint32)
    - False
    - Map of chromosome positions to crosslink counts

.. _createMatrix:
:c1:`createMatrix`
------------------

This function creates R friendly output matrices. 
It was seen that in some CLIP experiments the adjacent overlapping windows in the sliding windows file are duplicates of each other. 
In such cases, downstream analysis of these duplicate windows are not at all useful. As a workaround, while aggregating crosslink counts 
across all samples, this functon also filters out overlapping windows with duplicate crosslink counts across all samples, and keeps only the 5' most window (in gene direction).
This behavior can be overridden by using ``--allow_duplicates`` flag.


The following output files are generated:

* count matrix: aggregated crosslink sites per window, per sample
* annotation table : flattened annotation file with unique ids, gene names, gene features, gene ids,  etc.
* max count matrix (optional): maximum crosslink sites per window, per sample

Annotation table follows the following schema:

.. _AnnotationTable:

.. list-table::
    :header-rows: 1
  
    * - Column
      - Description
    * - unique_id
      - Unique identifier of this window
    * - chromosome
      - Chromosome name
    * - begin
      - window start position
    * - end
      - window end position
    * - strand
      - Strand info
    * - gene_id
      - Gene identifier
    * - gene_name
      - Gene name
    * - gene_type
      - Gene type, eg: protein coding, lncRNA,...
    * - gene_region
      - Gene region type, intron or exon
    * - Nr_of_region
      - number of the current region
    * - Total_nr_of_region
      - total number of regions
    * - window_number
      - (optional) Number of this window in the region



  
Further analysis
*****************

Further analysis and processing of crosslink windows is done using R/Bioconductor package `DEWSeq`_. Please refer to the
user manual of this package for requirements, installation and help. 


.. _`DEWSeq`: https://bioconductor.org/packages/release/bioc/html/DEWSeq.html
.. _`Gencode V42`: https://www.gencodegenes.org/human/release_42.html
.. _`htseq-clip`: https://htseq-clip.readthedocs.io/en/latest/overview.html
.. _`tabix`: https://www.htslib.org/doc/tabix.html
.. _`Apache parquet`: https://parquet.apache.org/
