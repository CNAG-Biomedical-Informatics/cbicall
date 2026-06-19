import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

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

:::note[WES interval source]
The interval resource depends on the software stack. Legacy `gatk-3.5` Bash WES
uses Agilent SureSelect hg19 BED files, while current `gatk-4.6` native WES uses
the GATK bundle / Broad b37 exome interval list. See the
[FAQ](../help/faq#wes--wgs) for details.
:::

---

## Execution Modes

<Tabs groupId="cohort-execution-mode">
<TabItem value="standard" label="Standard run" default>

Standard cohort mode runs as one job from GenomicsDB import through final
filtering.

### 1. GenomicsDBImport

- Imports all gVCFs into a **GenomicsDB workspace**.
- Handles both WES (interval-limited) and WGS (whole genome) modes.
- Output: on-disk database under `01_genomicsdb/`, accessed as
  `gendb://<workspace>`.

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
| `QD2` | `QD < 2.0`, when `QD` is present |
| `FS60` | `FS > 60.0` |
| `MQ40` | `MQ < 40.0` |
| `MQRS-12.5` | `MQRankSum < -12.5`, when `MQRankSum` is present |
| `RPRS-8` | `ReadPosRankSum < -8.0`, when `ReadPosRankSum` is present |
| `QD2_indel` | `QD < 2.0`, when `QD` is present |
| `FS200` | `FS > 200.0` |
| `RPRS-20` | `ReadPosRankSum < -20.0`, when `ReadPosRankSum` is present |

This QC VCF is the primary cohort workflow output for downstream tools and
project-level review.

</TabItem>
<TabItem value="staged" label="Staged run">

Staged execution splits a large native GATK 4.6 cohort run into shard jobs and
one final filtering job.

For large WES/WGS cohorts, a single GenomicsDB import and joint-genotyping job
can be slow or memory-heavy. Native GATK 4.6 cohort workflows can split the
cohort run into two explicit stages:

| Stage | What it does | Main output |
| --- | --- | --- |
| `all` | Default behavior: import, genotype, VQSR/hard-filter in one run. | `02_varcall/<basename>.gv.QC.vcf.gz` |
| `shard` | Import and genotype one interval shard, then stop before VQSR/filtering. | `02_varcall/<basename>.gv.raw.vcf.gz` |
| `finalize` | Start from a user-gathered raw cohort VCF and run global VQSR/hard-filtering. | `02_varcall/<basename>.gv.QC.vcf.gz` |

This is intentionally a manual checkpoint. Run one `shard` job per contig or
interval shard, confirm that all shard jobs finished, concatenate the raw VCFs,
then run `finalize` once on the gathered raw VCF. Final filtering stays global,
which is easier to reason about than applying VQSR independently to each shard.

:::info[Backend support]
Staged cohort execution is supported for CBIcall-native GATK 4.6 cohort runs
with `workflow_backend: bash`, `snakemake`, `nextflow`, or `cromwell`.
:::

### Shard runs

Each shard run still uses the full `sample_map`; the shard only restricts the
genomic intervals sent to `GenomicsDBImport` and `GenotypeGVCFs`.

```yaml
mode:            cohort
pipeline:        wgs
workflow_backend: bash
software_stack:    gatk-4.6
genome:          hg38
sample_map:      ./sample_map.tsv
cohort_stage:    shard
interval_shard:  chr1
output_basename: cohort.chr1
```

For WGS, `interval_shard` must match a contig in the reference dictionary, such
as `chr1` for hg38. For WES, CBIcall filters the configured WES interval list to
records whose contig column matches `interval_shard`.

:::note[GenomicsDB workspace names]
Do not set a GenomicsDB workspace in the parameters YAML. CBIcall creates a unique
workspace for each run under `01_genomicsdb/cohort.genomicsdb.<run-id>`, so
parallel shard jobs do not collide. Use `output_basename` to control the VCF
stems users see, for example `cohort.chr1`.
:::

Run one YAML per shard, changing `interval_shard` and `output_basename`:

```bash
bin/cbicall run -p cohort.chr1.yaml -t 12
bin/cbicall run -p cohort.chr2.yaml -t 12
```

#### Run Chromosome Shards With GNU Parallel

On a workstation, GNU parallel is a convenient way to launch a small number of
chromosome shard jobs concurrently. On HPC, use the same YAML keys but submit
each chromosome as a scheduler job instead of running GNU parallel on the login
node.

<details>
<summary>Workstation GNU parallel example</summary>

This example launches one WES b37 shard per chromosome with two concurrent
CBIcall runs (`-j 2`). It assumes one sample map per chromosome, named
`sample_map.chr1.txt` ... `sample_map.chr22.txt`, `sample_map.chrX.txt`, and
`sample_map.chrY.txt`. Each sample map should contain absolute gVCF paths.

```bash
#!/usr/bin/env bash
set -euo pipefail

export CBICALL=/path/to/cbicall/bin/cbicall
export ROOT=/path/to/project

parallel --halt soon,fail=1 --joblog cbicall.shards.joblog -j 2 \
  --env CBICALL --env ROOT '
  set -euo pipefail
  chr={1}
  yaml="cbicall.chr${chr}.yaml"

  printf "%s\n" \
    "mode: cohort" \
    "pipeline: wes" \
    "workflow_backend: bash" \
    "software_stack: gatk-4.6" \
    "genome: b37" \
    "sample_map: ${ROOT}/sample_map.chr${chr}.txt" \
    "cohort_stage: shard" \
    "interval_shard: ${chr}" \
    "output_basename: cohort.chr${chr}" \
    > "$yaml"

  "$CBICALL" -p "$yaml" -t 4
' ::: $(seq 22 -1 1) X Y
```

</details>

For WES b37, `interval_shard` uses bare contig labels (`1`, `2`, ..., `22`,
`X`, `Y`) because the bundled b37 exome interval list uses bare contig names.
For WGS hg38, use reference-dictionary labels such as `chr1`.

After all shard jobs finish, concatenate the raw shard VCFs in genomic order:

```bash
set -euo pipefail

for chr in $(seq 1 22) X Y; do
  ls -1 cbicall_*/02_varcall/cohort.chr${chr}.gv.raw.vcf.gz
done > raw_vcfs.list

bcftools concat -f raw_vcfs.list -Oz -o cohort.gathered.gv.raw.vcf.gz
bcftools index -t cohort.gathered.gv.raw.vcf.gz
```

### Finalize run

The finalize run does not need `sample_map`. It starts from the gathered raw VCF
and produces the same kind of filtered cohort output as the default `all` run.

```yaml
mode:            cohort
pipeline:        wgs
workflow_backend: bash
software_stack:    gatk-4.6
genome:          hg38
cohort_stage:    finalize
input_vcf:       ./cohort.gathered.gv.raw.vcf.gz
output_basename: cohort
```

```bash
bin/cbicall run -p cohort.finalize.yaml -t 12
```

</TabItem>
</Tabs>

---

## Output Files

In these filenames, `gv` means `GenotypeGVCFs`: the raw VCF is the direct
joint-genotyped output from that GATK step, and the QC VCF is the filtered
version used downstream.

The default basename is `cohort`. If `output_basename` is set, replace `cohort`
below with that value.

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
