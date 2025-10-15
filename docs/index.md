<div align="center">
<strong>Welcome to the documentation for CBICall</strong>
</div>

<p align="center">
  <a href=""><img src="img/cbicall-logo.png" width="300" alt="CBICall"></a>
</p>

**CBICall** (**C**NAG **B**iomedical **I**nformatics framework for variant **Call**ing) is a lightweight, reproducible framework for germline variant calling developed at CNAG. Built for Illumina WES/WGS data, CBICall wraps established best-practices (BWA → GATK → VQSR / hard-filters) into easy-to-run **Bash** and **Snakemake** workflows so labs can produce high-quality single-sample and cohort VCFs with minimal fuss. 🧬

## Why CBICall?
- Implements GATK best practices (GATK-4.6 and legacy GATK3.5) tuned for real project needs.  
- Supports both single-sample and cohort pipelines (WES / WGS) with GenomicsDBImport and optional per-chromosome sharding.  
- Handles mitochondrial DNA (mtDNA) analysis through MToolBox integration for heteroplasmy-aware calling and mtDNA-specific processing.  
- Simple YAML configuration and sensible defaults to get you running quickly.  
- Transparent, auditable logs and outputs for QC and downstream analysis.

## Key features
- Per-sample preprocessing: alignment, read groups, merging, duplicate marking, BQSR.  
- Per-sample GVCF generation and scalable joint genotyping (GenomicsDBImport → GenotypeGVCFs).  
- Variant quality control: VQSR when cohort size permits, with reproducible hard-filter fallbacks.  
- mtDNA support: MToolBox-based workflows for mitochondrial assembly, heteroplasmy estimation and annotation.  
- Small-footprint workflows: Bash and Snakemake variants, plus optional containerized deployment.  
- Handy utility scripts: coverage stats, sex determination, basic cohort QC.
