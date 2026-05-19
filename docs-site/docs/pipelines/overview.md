# Included Pipelines

CBIcall ships curated variant-calling pipelines for Illumina DNA sequencing.
These are the analyses available out of the box.

| Included pipeline | Modes | Main output |
| --- | --- | --- |
| WES | Single-sample and cohort joint genotyping | Germline VCF/gVCF outputs for exome data |
| WGS | Single-sample and cohort joint genotyping | Germline VCF/gVCF outputs for genome data |
| mtDNA | Single-sample and cohort/family analysis | mtDNA VCF, prioritized variant table, and browser report |

The selected pipeline is configured with `pipeline` and `mode`:

```yaml
pipeline: wes
mode: single
```

The workflow backend is configured separately with `workflow_engine`. See
[Native Backends](../backends/native) for Bash, Snakemake, Nextflow, and
external nf-core execution.

## Pipeline Guides

| Goal | Page |
| --- | --- |
| Process one WES/WGS sample from FASTQ | [WES/WGS Single-Sample](wes-wgs-single) |
| Joint-genotype a WES/WGS cohort from gVCFs | [WES/WGS Cohort](wes-wgs-cohort) |
| Run mitochondrial variant calling | [mtDNA](mtdna) |

## Backend Availability

Backend coverage differs by pipeline family:

| Pipeline | Bash | Snakemake | Nextflow |
| --- | --- | --- | --- |
| WES | Yes | Yes | Yes |
| WGS | Yes | Yes | Yes |
| mtDNA | Yes | No | No |

External nf-core workflows are registered through the Nextflow backend, but they
are not part of the included WES/WGS/mtDNA pipeline set.
