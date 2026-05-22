import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# WES/WGS Single-Sample Pipeline
A user-focused guide to processing whole-exome (WES) and whole-genome (WGS) data using GATK Best Practices.


<Tabs groupId="workflow-mode">
<TabItem value="bash" label="Bash" default>

**Source:** [View source](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-4.6/wes_single.sh)

</TabItem>
<TabItem value="snakemake" label="Snakemake">

**Source:** [View source](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/snakemake/gatk-4.6/wes_single.smk)

</TabItem>
<TabItem value="nextflow" label="Nextflow">

**Source:** [View source](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/nextflow/gatk-4.6/wes_single.nf)

</TabItem>
</Tabs>
---

## Diagram: Single-Sample WES/WGS Workflow

![WES/WGS single-sample workflow](/img/diagram-wes-wgs-single.svg)

---

## Purpose

This pipeline processes **one sample** at a time and produces a filtered VCF
and gVCF organized for downstream tools and project QC. It automatically adapts
to:

- **WES**: restricted to an exome interval list
- **WGS**: whole genome (no interval restriction)

---

## What the Pipeline Does

### 1. Alignment & Read Groups

- Align paired-end FASTQ files using **BWA-MEM**.
- Add read groups (sample, library, lane, platform) required by GATK.
- Output: lane-level BAMs with correct RG tags.

### 2. Lane Merging

- Merge all lane BAMs for the same sample into a single BAM.
- Ensures duplicate marking and BQSR operate on the full dataset.

### 3. Duplicate Marking

- Use GATK `MarkDuplicates` on the merged BAM.
- Flags PCR/optical duplicates to prevent them from inflating support for artefactual variants.

### 4. Base Quality Score Recalibration (BQSR)

- Two-step process: `BaseRecalibrator` then `ApplyBQSR`.
- Uses known variant databases (dbSNP, Mills, 1000G indels) to model and correct systematic base-quality errors.
- Output: recalibrated BAM used for variant calling.

### 5. Variant Calling (HaplotypeCaller, gVCF)

- Run GATK `HaplotypeCaller` in **GVCF mode** (`-ERC GVCF`).
- WES: uses exome intervals; WGS: full genome.
- Output: `<id>.hc.g.vcf.gz` (per-sample gVCF).

### 6. GenotypeGVCFs (Raw VCF)

- Run GATK `GenotypeGVCFs` on the sample gVCF.
- Output: `<id>.hc.raw.vcf.gz` (raw VCF with SNPs and indels).

### 7. Variant Quality Score Recalibration (VQSR)

- If there are enough variants (SNPs and indels), build VQSR models:
  - `VariantRecalibrator` for SNPs and indels separately.
  - Uses multiple annotations (QD, MQ, FS, MQRankSum, ReadPosRankSum).
- Output: recalibration VCFs and tranche files.

### 8. Apply VQSR or Hard Filters

- If models exist:
  - Apply SNP VQSR.
  - Then apply INDEL VQSR.
  - Output: `<id>.hc.vqsr.vcf.gz`.
- If not:
  - Skip directly to hard filters on the raw VCF or post-SNP VQSR VCF.

### 9. Generate Final QC VCF

- Run `VariantFiltration` with the GATK 4.6 hard filters below.
- Output: `<id>.hc.QC.vcf.gz` (final QC VCF).

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

### 10. Coverage & Sex Determination

- Extract chromosome 1 reads from raw and recalibrated BAMs.
- Compute coverage statistics.
- Infer sample sex from final VCF using a dedicated script.
- Outputs:
  - `03_stats/<id>.coverage.txt`
  - `03_stats/<id>.sex.txt`

---

## Output Files

| File                                   | Meaning                                             |
|----------------------------------------|-----------------------------------------------------|
| `02_varcall/<id>.hc.g.vcf.gz`          | Per-sample gVCF (HaplotypeCaller)                   |
| `02_varcall/<id>.hc.raw.vcf.gz`        | Raw VCF after GenotypeGVCFs                         |
| `02_varcall/<id>.hc.vqsr.vcf.gz`       | VCF after VQSR (if VQSR was applied)                |
| `02_varcall/<id>.hc.QC.vcf.gz`         | Final QC-filtered VCF (recommended)                 |
| `03_stats/<id>.coverage.txt`           | Coverage metrics                                    |
| `03_stats/<id>.sex.txt`                | Sex determination result                            |
| `logs/<id>.log`                        | Main pipeline log                                   |

---

## When to Use This Pipeline

- Standard research **WES** or **WGS** processing.
- Generating gVCFs for **cohort joint genotyping**.
- Producing filtered single-sample VCFs for downstream review or interpretation.
