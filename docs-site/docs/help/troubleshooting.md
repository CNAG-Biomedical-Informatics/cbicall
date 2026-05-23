# Troubleshooting

Use the error text from the terminal or workflow log to find the matching section. Most failures fall into three groups: missing external data, GATK/Picard input problems, or mtDNA-specific MToolBox issues.

:::tip[Identify the execution path]
Before debugging, confirm whether the run is native CBIcall or external nf-core.

| Run type | YAML signal | First place to look |
| --- | --- | --- |
| Native CBIcall | `workflow_provider: cbicall` or omitted | CBIcall run log, `log.json`, and backend logs such as `logs/*.log` |
| External nf-core | `workflow_provider: nf-core` | Nextflow launcher log, nf-core output directory, `pipeline_info/`, and task work directories |

For all runs, check `log.json` to confirm the resolved `input_dir`, `sample_map`, `genome`, workflow, resource, and run directory.
:::

## Installation and External Data

<details>
<summary>External data or tool path not found</summary>

**Symptom**

```text
/usr/bin/bash: line 9: /media/mrueda/2TBS/NGSutils/gatk/gatk-4.6.2.0/gatk: No such file or directory
```

**Likely cause**

`DATADIR` does not point to the directory where databases and external tools are installed or mounted.

**Fix**

Update the data directory in the workflow configuration:

```text
workflows/bash/gatk-4.6/env.sh
workflows/snakemake/gatk-4.6/config.yaml
workflows/nextflow/gatk-4.6/config.yaml
workflows/cromwell/gatk-4.6/config.yaml
```

For containers, make sure the host data directory is bind-mounted at the same path used by the workflow configuration. For external nf-core workflows, missing tools or references are usually controlled by the selected nf-core profile, container/cache setup, or nf-core parameters rather than the CBIcall resource bundle.

</details>

<details>
<summary>Relative input paths resolve somewhere unexpected</summary>

**Symptom**

CBIcall cannot find FASTQ files, BAMs, or `sample_map.tsv`, even though the path looks correct from your current shell.

**Likely cause**

Relative `input_dir` and `sample_map` paths are resolved from the YAML file location.

**Fix**

Use absolute paths, or keep the YAML file next to the relative paths it references. Confirm the resolved paths in `log.json`.

</details>

## GATK and Picard

<details>
<summary>NaN LOD value during recalibration</summary>

**Symptom**

```text
NaN LOD value assigned
```

**Likely cause**

There are too few variants to train a reliable VQSR model, often too few INDELs.

**Fix**

Use the existing thresholds that skip VQSR when the variant count is too small, or increase the minimum threshold before rerunning. The final `*.QC.vcf.gz` is still produced by hard filtering when VQSR is skipped.

</details>

<details>
<summary>Not enough columns in dbSNP line</summary>

**Symptom**

```text
there aren't enough columns for line ... dbsnp_137.hg19.vcf
```

**Likely cause**

The dbSNP VCF contains malformed or truncated records.

**Fix**

Inspect the reported line in the dbSNP VCF, replace the database file if possible, or correct the malformed record locally and document the change.

</details>

<details>
<summary>Error parsing text SAM file</summary>

**Symptom**

```text
Error parsing text SAM file. Not enough fields; File /dev/stdin; Line ...
```

**Likely cause**

Secondary or supplementary alignments can introduce records that Picard/GATK rejects when the alignment stream is passed directly into read-group assignment.

**Fix**

Filter secondary and supplementary alignments before adding read groups:

```bash
bwa mem -M -t "$THREADS" "$REFGZ" "$R1" "$R2" \
  | samtools view -bSh -F 0x900 - \
  | gatk AddOrReplaceReadGroups ...
```

</details>

## mtDNA and MToolBox

<details>
<summary>MToolBox fails on ARM / aarch64</summary>

**Symptom**

```text
mit_single cannot be performed with: aarch64
```

or:

```text
mit_cohort cannot be performed with: aarch64
```

**Likely cause**

The bundled MToolBox workflow is x86_64-only.

**Fix**

Run mtDNA workflows on an x86_64 Linux host. WES/WGS GATK 4.6 workflows can still run on supported ARM systems.

</details>

<details>
<summary>No usable BAM found for mtDNA</summary>

**Symptom**

```text
ERROR: Could not find BAM for ID ...
```

or:

```text
ERROR: No usable sample BAMs found. Nothing to do.
```

**Likely cause**

The mtDNA workflow expects BAMs from previous WES/WGS single-sample runs in the expected project layout.

**Fix**

Run WES/WGS single-sample processing first, keep the `01_bam` outputs, and then rerun the mtDNA workflow from the sample or project directory described in the mtDNA example.

</details>

<details>
<summary>Unsupported N CIGAR operations</summary>

**Symptom**

MToolBox fails due to unsupported `N` operations in CIGAR strings.

**Likely cause**

Some reads contain skipped-region CIGAR operations that MToolBox cannot process.

**Fix**

Add this flag in the relevant MToolBox alignment or SAM-processing step:

```text
--filter_reads_with_N_cigar
```

</details>

<details>
<summary>Low coverage and unreliable heteroplasmy fractions</summary>

**Symptom**

mtDNA coverage is low, or heteroplasmy fraction estimates look unstable.

**Likely cause**

Below roughly 10x mtDNA coverage, heteroplasmy fraction estimates are unreliable.

**Fix**

Flag samples below 10x median mtDNA coverage, interpret HF values cautiously, exclude low-coverage samples from HF-based analyses when needed, and consider resequencing if mtDNA interpretation is critical.

</details>

## Variant Interpretation

<details>
<summary>Unexpected de novo rates in trios</summary>

**Symptom**

Observed de novo rates differ strongly from expectations.

**Likely cause**

Large deviations can indicate sample, data-quality, annotation, or pipeline issues.

**Reference values**

| Sample type | Typical de novo rate |
| --- | ---: |
| Proband | ~1% |
| Parent | ~10% |

**Fix**

Check sample identity, pedigree labels, coverage, variant filters, and annotation assumptions before interpreting the result biologically.

</details>

## Next Steps

- Confirm generated files in [Outputs](outputs).
- Confirm YAML choices in [Configuration Reference](configuration-reference).
- Review runtime guidance in [Performance](performance).
