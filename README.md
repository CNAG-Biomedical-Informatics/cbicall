<div align="center">
  <a href="https://github.com/mrueda/cbicall">
    <img src="https://raw.githubusercontent.com/mrueda/cbicall/main/docs/img/cbicall-logo.png"
         width="350" alt="CBICall">
  </a>
  <p><em>CNAG Biomedical Informatics framework for variant Calling</em></p>
</div>

![version](https://img.shields.io/badge/version-0.0.1-28a745)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)


# Table of contents
- Description
  - [Name](#name)
  - [Synopsis](#synopsis)
  - [Summary](#summary)
- [Installation](#installation)
  - [Non-Containerized](non-containerized/README.md)
  - [Containerized](docker/README.md)
- [How to run cbicall](#how-to-run-cbicall)
- [Citation](#citation)
  - [Author](#author)
- [License](#copyright-and-license)

# NAME

CBICall: CNAG Biomedical Informatics Framework for Variant Calling on Illumina DNA-seq (germline) NGS Data.

# SYNOPSIS

    cbicall -p <parameters_file.yaml> -t <n_threads> [options]

    Arguments:
      -p|param          Parameters input file (YAML)
      -t|threads        Number of CPUs/Cores/Threads

    Options:
      -debug            Debugging level (from 1 to 5; 5 is maximum verbosity)
      -h|help           Brief help message
      -man              Full documentation
      -v                Show version information
      -verbose          Enable verbose output
      -nc|no-color      Do not print colors to STDOUT
      -ne|no-emoji      Do not print emojis to STDOUT

# SUMMARY

CBICall (**C**NAG **B**iomedical **I**nformatics framework for variant **Call**ing) is a computational framework designed for variant calling analysis using Illumina Next-Generation Sequencing (NGS) data.

# HOW TO RUN CBICALL

CBICall execution requires:

- Input Files

    A folder containing Paired-End FASTQ files (e.g., `MA00001_exome/MA0000101P_ex/*{R1,R2}*fastq.gz`).

    You have a `examples/input/` directory with input data that you can use for testing.

- Parameters File

    A YAML-formatted parameters file controlling pipeline execution.

Below are the parameters that can be customized, along with their default values. Parameters must be separated from their values by whitespace or tabs.

## Essential Parameters

    mode:            single  
    pipeline:        wes          
    sample:          undef        
    sample_map:      undef
    workflow_engine:   bash
    gatk_version:      gatk3.5
    cleanup_bam:       false

## Optional Parameters (Currently Unused)

    organism:        Homo Sapiens        
    technology:      Illumina HiSeq      

CBICall will create a dedicated project directory (`cbicall_*`) to store analysis outputs. This design allows multiple independent runs concurrently without modifying original input files.

Below is a detailed description of key parameters:

- **cleanup\_bam**

    Set it to `true` to delete `01_bam/*.{bam,bai}`.

- **gatk\_version**

    Supported values: `gatk3.5` or `gatk4.6`.

- **mode**

    Two modes are supported: `single` (default, for individual samples) and `cohort` (for family/cohort-based analyses).

- **pipeline**

    Specifies the analysis pipeline. Currently available options: `wes` (whole-exome sequencing) and `mit` (mitochondrial DNA analysis). Note: to run `cohort` analysis, first complete a `single` analysis for each sample.

- **projectdir**

    The prefix for dir name (e.g., 'cancer\_sample\_001'). Note that it can also contain a path (e.g., foo/cancer\_sample\_001).

    **Note:** Such directory will be always created below the **sample** directory. The script will automatically add an unique identifier to each job.

- **sample**

    Path (relative or absolute) to the directory containing FASTQ files for analysis. See the `examples` directory for detailed guidance.

    Example:

    examples/input/CNAG999\_exome/CNAG99901P\_ex

- **sample\_map** (cohort-mode only)

    Path (relative or absolute) to the file containing the sample ids and teh paths for the GVCF files

    See example [here](https://github.com/CNAG-Biomedical-Informatics/cbicall/tree/main/examples/input/sample_map.tsv)

- **workflow\_engine**

    Supported workflow engines: `bash` or `snakemake`.

## Example Commands

    $ bin/cbicall -p param_file.yaml -t 8
    $ bin/cbicall -p param_file.yaml -t 4 -verbose
    $ bin/cbicall -p param_file.yaml -t 16 > log 2>&1
    $ $path_to_cbicall/bin/cbicall -p param_file.yaml -t 8 -debug 5

Note: For Trio analyses, unique (de novo) variant rates for probands typically should be ~1%, and ~10% for parents. Significant deviations may indicate issues.

## ADDENDUM: Nomenclature Guidelines

All parts must follow a strict character count, and everything after the underscore is mandatory.

## Directory Naming

- Format: `[ProjectCode]_[SampleType]`

    - - `ProjectCode`: Exactly 7 characters \[a-zA-Z0-9\] (e.g., `MA99999`)
    - - `SampleType`: Must be `exome` (5 characters)

    Example:

        MA99999_exome

    Total: 7 + 1 + 5 = 13 characters.

## Subdirectory Naming

- Format: `[ProjectCode][SampleID][Role]_[SampleTypeShort]`

    - - `ProjectCode`: 7 characters (\[a-zA-Z0-9\] e.g., `MA99999`)
    - - `SampleID`: 2 characters (e.g., `01`)
    - - `Role`: 1 character (e.g., `P` for Proband, `F` for Father, `M` for Mother)
    - - `SampleTypeShort`: Must be `ex` (2 characters)

    Example:

        MA9999901F_ex

    Total: 7 + 2 + 1 + 1 + 2 = 13 characters (excluding any file extension).

## FASTQ Naming Convention

This convention is adapted from the following document:

[Document](https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)

We have added a custom suffix to indicate sequencing type:

    - _ex for exome sequencing
    - _wg for whole genome sequencing

Example:

    MA0004701P_ex_S5_L001_R1_001.fastq.gz

In summary, you need to have something like this:

    MA00001_exome/MA0000101P_ex/MA0000101P_ex_S1_L001_R?_001.fastq.gz

# SYSTEM REQUIREMENTS

CBICall is optimized for multi-core Linux desktop, workstation, or server environments. Snakemake-based workflows can also be adapted for HPC clusters.

Recommended specifications:

    * Works in amd64 and arm64 archs (M-based Macs).
    * Ideally a Debian-based distribution (Ubuntu or Mint), but any other (e.g., CentOS, OpenSUSE) should do as well (untested).
    * >= 8 GB RAM.
    * >= 4 CPU cores (Intel i7 or Xeon preferred).
    * >= 250 GB HDD space.
    * Perl >= 5.36 and required CPAN modules (install via C<cpanm --notest --installdeps .>).
    * Java 8 (install via C<sudo apt install openjdk-8-jdk>).
    * Snakemake (install via C<pip3 install -r requirements.txt>).

Perl scripts in CBICall use minimal RAM (~2% of a 16 GB system). Genome mapping with BWA benefits from higher memory but lacks built-in RAM limits. Its usage depends on thread count and reference size. To constrain BWA's memory, external tools like shell `ulimit` are required. In contrast, GATK and Picard default to 8 GB RAM, adjustable via the configuration file.

Parallel execution is supported but does not scale linearly. Optimal performance is achieved using ~ 4 threads per task. For example, with 12 cores, running 3 tasks in parallel with 4 cores each is typically more efficient than one task with all 12 cores. See example in figure below:

![Time](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/docs/img/run-time.png)

Unit/integration tests are conducted manually by verifying CSV and VCF outputs against established test datasets.

# SUPPORTED PIPELINES

The following table shows valid pipeline and mode combinations for each GATK version:

| GATK Version | wes\_single | wes\_cohort | wgs\_single | wgs\_cohort | mit\_single | mit\_cohort |
|--------------|------------|------------|------------|------------|------------| -----------|
| gatk3.5      | +          | +          | -          | -          |   +        | +          | 
| gatk4.6      | +          | +          | +          | +          |   -        | -          | 

Date: Oct-2025

## Capture Kits

\* For GATK version 3.5: Exome capture is based on Agilent SureSelect.

\* For GATK version 4.6: Exome and WGS reference is based on the GATK bundle (b37).

# COMMON ERRORS AND TROUBLESHOOTING

- GATK|Picard Errors (wes\_single.sh or wes\_cohort.sh)
    - Error: `NaN LOD value assigned` in recalibration steps.

        Occurs due to insufficient INDEL variants (typically fewer than 8000) for negative model training. The default threshold is 8000.

        Solution: Increase minimum INDEL count (e.g., to >8000) in relevant pipeline scripts. Only rerun failed samples.

    - Error: `there aren't enough columns for line ... dbsnp_137.hg19.vcf`

        Solution: Remove the problematic line from the VCF file and document changes in a README file.

    - Error: `Error parsing text SAM file. Not enough fields; File /dev/stdin; Line 105120626...`

        Certain SRA-derived or dbGaP datasets can contain duplicate reads.

        When piping BWA output into `AddOrReplaceReadGroups`, you may need to remove secondary (`0x100`) and supplementary (`0x800`) alignments.

        This can prevent duplicate-read collisions in downstream Picard/GATK steps.

        **Solution**: Uncomment the following line in `wes_single.sh`:

        `| $SAM view -bSh -F 0x900 -`
- **MTOOLBOX Errors**

    \- Failure related to unsupported N\_CIGAR:  
      Add flag `--filter_reads_with_N_cigar` in Mtoolbox.sh (line ~386).

    \- Samples with coverage below ~10x yield unreliable heteroplasmy fractions (HF). Extremely low coverage (<10x) can render HFs meaningless, despite generally consistent results across widely varying coverage levels.

# CITATION

To be determined.

# AUTHOR

Written by Manuel Rueda (mrueda). GitHub repository: [https://github.com/CNAG-Biomedical-Informatics/cbicall](https://github.com/CNAG-Biomedical-Informatics/cbicall).

# COPYRIGHT AND LICENSE

Please see the included LICENSE file for distribution and usage terms.
