# WES/WGS Cohort Joint-Genotyping Pipeline

A user-oriented guide for multi-sample joint genotyping using GenomicsDB, GenotypeGVCFs, and VQSR.


**Sources:** [Bash](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-4.6/wes_cohort.sh), [Snakemake](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/snakemake/gatk-4.6/wes_cohort.smk), [Nextflow](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/nextflow/gatk-4.6/wes_cohort.nf), [Cromwell](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/cromwell/gatk-4.6/wes_cohort.wdl)

---

## Diagram: Cohort Joint-Genotyping Workflow

![WES/WGS cohort joint-genotyping workflow](/img/diagram-wes-wgs-cohort.svg)

---

## Purpose

This pipeline combines many per-sample gVCFs and performs **cohort-level joint genotyping**, producing a single VCF for all samples with consistent genotype calls and variant filtering.

Use this when:

- You want consistent genotypes across a family or population.
- You plan to build VQSR models on cohort-level variant distributions.
- You are preparing a joint VCF for association or segregation analyses.

---

## Inputs

- **Sample map file** (`sample_map` in the CBIcall YAML):
  TSV used by `GenomicsDBImport` (sample name → gVCF path).
- **Per-sample gVCFs** from the single-sample pipeline.
- **Reference genome** (`REF` from `env.sh`).
- **VQSR resources** (SNP and INDEL training sets).
- Optional **interval list** for WES mode.

---

## Workflow

### 1. GenomicsDBImport

- Imports all gVCFs into a **GenomicsDB workspace**.
- Handles both WES (interval-limited) and WGS (whole genome) modes.
- Output: on-disk database accessed as `gendb://<workspace>`.

### 2. Joint Genotyping (GenotypeGVCFs)

- Runs `GenotypeGVCFs` on `gendb://<workspace>` and the reference.
- Produces the cohort-level VCF:

  - `cohort.gv.raw.vcf.gz`

### 3. Count Variants and Decide on VQSR

- Counts SNPs and INDELs in the raw cohort VCF.
- Compares counts to configurable thresholds:
  - `MIN_SNP_FOR_VQSR` (default 1000)
  - `MIN_INDEL_FOR_VQSR` (default 8000)
- Determines whether to build SNP and/or INDEL VQSR models.

### 4. Build VQSR SNP Model

- If enough SNPs:
  - Run `VariantRecalibrator` in SNP mode.
  - Uses training resources and multiple annotations:
    - QD, FS, MQ, MQRankSum, ReadPosRankSum.
  - Outputs:
    - `cohort.snp.recal.vcf.gz`
    - `cohort.snp.tranches.txt`.

### 5. Build VQSR INDEL Model

- If enough INDELs:
  - Run `VariantRecalibrator` in INDEL mode.
  - Uses annotations:
    - QD, FS, ReadPosRankSum.
  - Outputs:
    - `cohort.indel.recal.vcf.gz`
    - `cohort.indel.tranches.txt`.

### 6. Apply VQSR

- If SNP model exists:
  - Apply SNP VQSR → `cohort.post_snp.vcf.gz`.
- If INDEL model exists:
  - Apply INDEL VQSR → `cohort.vqsr.vcf.gz`.

The best available VCF (VQSR-filtered or raw) is used as input to the next step.

### 7. Hard Filtering and QC VCF

- Run `VariantFiltration` with the GATK 4.6 hard filters below.
- Output: `cohort.gv.QC.vcf.gz`.

| Filter name | Expression |
| --- | --- |
| `LowQUAL` | `QUAL < 30.0` |
| `QD2` | `QD < 2.0` |
| `FS60` | `FS > 60.0` |
| `MQ40` | `MQ < 40.0` |
| `MQRS-12.5` | `MQRankSum < -12.5` |
| `RPRS-8` | `ReadPosRankSum < -8.0` |
| `QD2_indel` | `QD < 2.0` |
| `FS200` | `FS > 200.0` |
| `RPRS-20` | `ReadPosRankSum < -20.0` |

This QC VCF is the primary cohort workflow output for downstream tools and
project-level review.

---

## Output Files

| File                                | Description                                        |
|-------------------------------------|----------------------------------------------------|
| `cohort.gv.raw.vcf.gz`              | Raw cohort joint-genotyped VCF                    |
| `cohort.snp.recal.vcf.gz`           | SNP VQSR model VCF                                |
| `cohort.snp.tranches.txt`           | SNP VQSR tranches and diagnostics                 |
| `cohort.indel.recal.vcf.gz`         | INDEL VQSR model VCF                              |
| `cohort.indel.tranches.txt`         | INDEL VQSR tranches and diagnostics               |
| `cohort.post_snp.vcf.gz`            | VCF after applying SNP VQSR                       |
| `cohort.vqsr.vcf.gz`                | VCF after applying SNP and INDEL VQSR             |
| `cohort.gv.QC.vcf.gz`               | Final hard-filtered cohort QC VCF (recommended)   |
| `logs/cohort_joint_genotyping.log`  | Main log file for the cohort pipeline             |

---

## When to Use This Pipeline

- After you have gVCFs from the **single-sample** pipeline.
- When you need a single VCF for all samples for:
  - Family-based segregation analysis.
  - Case/control or population cohorts.
  - Downstream tools that expect joint genotypes.
