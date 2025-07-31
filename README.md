**BAM Realignment Pipeline**
Overview

This Nextflow pipeline implements the BAM realignment strategy described in the Human Pangenome Reference Consortium (HPRC) draft paper. It performs comprehensive realignment of sequencing reads to improve alignment accuracy, particularly around indels and complex genomic regions.
What the Pipeline Does

The pipeline takes a sorted BAM file and performs the following steps:

    Chromosome Splitting: Splits the input BAM file by chromosome for parallel processing
    Read Group Addition: Adds required read group information to BAM headers
    Left-alignment: Uses FreeBayes bamleftalign to left-align indels for consistency
    Target Identification: Identifies regions requiring realignment using GATK RealignerTargetCreator
    Region Expansion: Expands target regions by 160bp using bedtools slop
    ABRA Realignment: Performs assembly-based realignment using ABRA2 on target regions
    Merging: Combines per-chromosome realigned BAMs into a final output file

**Key Features**

    Parallel Processing: Processes all chromosomes simultaneously for efficiency
    HPRC Compliance: Follows the exact methodology from the Human Pangenome draft paper
    Robust Error Handling: Includes fallback mechanisms for problematic regions
    Comprehensive Output: Saves intermediate results from each processing step

**Tools Used**

    SAMtools (v1.3.1): BAM splitting, indexing, and merging
    FreeBayes (v1.2.0): Left-alignment of indels via bamleftalign
    GATK (v3.8): Realignment target identification
    bedtools (v2.30.0): Region expansion
    ABRA2 (v2.23): Assembly-based realignment
    Picard (v2.23.0): Read group management

**Input Requirements**

    Sorted and indexed BAM file
    Reference genome FASTA with index (.fai) and dictionary (.dict)
    Genome chromosome sizes file

**Output**

    Final realigned BAM file with improved alignment accuracy
    Per-chromosome intermediate files for quality control
    Comprehensive processing reports and execution metrics

Use Cases

This pipeline is ideal for:

    Preparing BAM files for variant calling
    Improving alignment quality around indels
    Processing whole-genome sequencing data
    Following HPRC best practices for pangenome analysis
