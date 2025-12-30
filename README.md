<div align="center">
  <a href="https://github.com/mrueda/cbicall">
    <img src="https://raw.githubusercontent.com/mrueda/cbicall/main/docs/img/cbicall-logo.png"
         width="350" alt="CBICall">
  </a>
  <p><em>CNAG Biomedical Informatics framework for variant Calling</em></p>
</div>


[![Build](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/build-and-test.yml)
![version](https://img.shields.io/badge/version-0.0.1-28a745)
[![Coverage Status](https://coveralls.io/repos/github/CNAG-Biomedical-Informatics/cbicall/badge.svg?branch=main)](https://coveralls.io/github/CNAG-Biomedical-Informatics/cbicall?branch=main)
[![Docker Build](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/docker-build-multi-arch.yml/badge.svg?branch=main)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/docker-build-multi-arch.yml)
[![Docker Image Size](https://img.shields.io/docker/image-size/manuelrueda/cbicall/latest?logo=docker&label=image%20size)](https://hub.docker.com/r/manuelrueda/cbicall/)
[![Docker Pulls](https://img.shields.io/docker/pulls/manuelrueda/cbicall.svg?logo=docker&label=pulls)](https://hub.docker.com/r/manuelrueda/cbicall/)
[![Documentation Status](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/documentation.yml/badge.svg)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/documentation.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

**📘 Documentation:** <a href="https://cnag-biomedical-informatics.github.io/cbicall" target="_blank">https://cnag-biomedical-informatics.github.io/cbicall</a>

**🐳 Docker Hub Image:** <a href="https://hub.docker.com/r/manuelrueda/cbicall/tags" target="_blank">https://hub.docker.com/r/manuelrueda/cbicall/tags</a>



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


# Name

CBICall: CNAG Biomedical Informatics Framework for Variant Calling on Illumina DNA-seq (germline) NGS Data.

# Synopsis

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

# Summary

CBICall (**C**NAG **B**iomedical **I**nformatics framework for variant **Call**ing) is a computational framework designed for variant calling analysis using Illumina Next-Generation Sequencing (NGS) data.

# How to run CBICall

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
    gatk_version:      gatk-4.6
    genome:            b37
    cleanup_bam:       false

## Optional Parameters (Currently Unused)

    organism:        Homo Sapiens        
    technology:      Illumina HiSeq      

CBICall will create a dedicated project directory (`cbicall_*`) to store analysis outputs. This design allows multiple independent runs concurrently without modifying original input files.

Below is a detailed description of key parameters:

- **cleanup\_bam**

    Set it to `true` to delete `01_bam/*.{bam,bai}`.

- **gatk\_version**

    Supported values: `gatk-3.5` or `gatk-4.6`.

- **genome**

    Supported values: `b37` (default), `hg38` or `rsrs` (mtDNA).

- **mode**

    Two modes are supported: `single` (default, for individual samples) and `cohort` (for family/cohort-based analyses).

- **pipeline**

    Specifies the analysis pipeline. Currently available options: `wes` (whole-exome sequencing), `wgs` (whole-genome sequencing) and `mit` (mitochondrial DNA analysis). Note: to run `cohort` analysis, first complete a `single` analysis for each sample.

- **projectdir**

    The prefix for dir name (e.g., 'cancer\_sample\_001'). Note that it can also contain a path (e.g., foo/cancer\_sample\_001).

    **Note:** Such directory will be always created below the **sample** directory. The script will automatically add an unique identifier to each job.

- **sample**

    Path (relative or absolute) to the directory containing FASTQ files for analysis. See the `examples` directory for detailed guidance.

    Example:

    examples/input/CNAG999\_exome/CNAG99901P\_ex

- **sample\_map** (cohort-mode only)

    Path (relative or absolute) to the file containing the sample ids and the paths for the GVCF files

    See example [here](https://github.com/CNAG-Biomedical-Informatics/cbicall/tree/main/examples/input/sample_map.tsv)

- **workflow\_engine**

    Supported workflow engines: `bash` or `snakemake`.

## Example Commands

    $ bin/cbicall -p param_file.yaml -t 8
    $ bin/cbicall -p param_file.yaml -t 4 -verbose
    $ bin/cbicall -p param_file.yaml -t 16 > log 2>&1
    $ $path_to_cbicall/bin/cbicall -p param_file.yaml -t 8 -debug 5
    $ nohup bin/cbicall -p param_file.yaml -t 4 &

# Recommended specifications:

CBICall is optimized for multi-core Linux desktop, workstation, or server environments.

    * Works in amd64 and arm64 archs (M-based Macs).
    * Ideally a Debian-based distribution (Ubuntu or Mint), but any other (e.g., CentOS, OpenSUSE) should do as well (untested).
    * >= 8 GB RAM.
    * >= 4 CPU cores (Intel i7 or Xeon preferred).
    * >= 250 GB HDD space.
    * Python >= 3.8
    * Java 8 (install via C<sudo apt install openjdk-8-jdk>).
    * Snakemake (install via C<pip3 install -r requirements.txt>).

# Citation

CBICall: a configuration-driven framework for variant calling in large DNA-seq cohorts. _Manuscript In preparation._

# Author

Written by Manuel Rueda (mrueda). GitHub repository: [https://github.com/CNAG-Biomedical-Informatics/cbicall](https://github.com/CNAG-Biomedical-Informatics/cbicall).

# Copyright and license

Please see the included LICENSE file for distribution and usage terms.

