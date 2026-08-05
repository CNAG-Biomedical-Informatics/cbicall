import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# End-to-end WES and WGS examples

These examples use the bundled `cbicall-core` GATK 4.6 workflows. Complete the
[installation](../installation/non-containerized), install the resource bundle,
and run `cbicall doctor` before starting.

<Tabs groupId="workflow-mode">
<TabItem value="single-sample" label="Single sample" default>

## 1. Prepare paired FASTQ files

Place the read pairs in one sample directory. Their names must share a prefix
and distinguish R1 from R2:

```text
CNAG999_exome/CNAG99901P_ex/
  CNAG99901P_ex_S1_L001_R1_001.fastq.gz
  CNAG99901P_ex_S1_L001_R2_001.fastq.gz
```

See [Naming Conventions](../help/naming-conventions) for supported file names
and multi-lane samples.

## 2. Create the parameters YAML

Create `wes_single.yaml`:

```yaml
mode:             single
pipeline:         wes
workflow_backend: bash
software_stack:   gatk-4.6
input_dir:        CNAG999_exome/CNAG99901P_ex
genome:           b37
cleanup_bam:      false
```

`workflow_backend` can be `bash`, `snakemake`, `nextflow`, or `cromwell` for
the bundled GATK 4.6 WES/WGS workflows. The selected backend must be installed
and available to CBIcall.

:::note[Running WGS]
Use a directory containing WGS FASTQ files and change `pipeline` to `wgs`.
Select `b37` or `hg38` to match the intended WGS reference resources:

```yaml
mode:             single
pipeline:         wgs
workflow_backend: bash
software_stack:   gatk-4.6
input_dir:        CNAG999_genome/CNAG99901P_wg
genome:           hg38
cleanup_bam:      false
```
:::

## 3. Run CBIcall

```bash
cbicall run -p wes_single.yaml -t 4
```

CBIcall prints the resolved workflow and run directory before launch. A
successful run ends with the workflow log and the paths to `run-report.json`
and `run-report.html`.

## 4. Inspect the outputs

The generated run directory contains:

```text
cbicall_bash_gatk-4.6_wes_single_b37_*/
  01_bam/
  02_varcall/
  03_stats/
  logs/
  log.json
  run-report.json
  run-report.html
```

The principal files are:

| File | Use |
| --- | --- |
| `02_varcall/<id>.hc.QC.vcf.gz` | Final filtered VCF. |
| `02_varcall/<id>.hc.g.vcf.gz` | Per-sample gVCF for later cohort joint genotyping. |
| `03_stats/<id>.coverage.txt` | Lightweight coverage summary. |
| `03_stats/<id>.sex.txt` | VCF-derived sex inference for QC. |
| `run-report.html` | Human-readable run and provenance report. |

</TabItem>
<TabItem value="cohort" label="Cohort">

## 1. Prepare per-sample gVCFs

Cohort mode starts from gVCFs produced by WES or WGS single-sample runs. Create
a two-column, tab-separated sample map with one absolute gVCF path per sample:

```text
SAMPLE_001	/absolute/path/SAMPLE_001.hc.g.vcf.gz
SAMPLE_002	/absolute/path/SAMPLE_002.hc.g.vcf.gz
```

Each gVCF must have its corresponding index.

## 2. Create the parameters YAML

Create `wes_cohort.yaml`:

```yaml
mode:             cohort
pipeline:         wes
workflow_backend: bash
software_stack:   gatk-4.6
genome:           b37
sample_map:       ./sample_map.tsv
```

Relative paths such as `sample_map` are resolved from the YAML file location.

## 3. Run joint genotyping

```bash
cbicall run -p wes_cohort.yaml -t 4
```

For large cohorts, use the documented `shard` and `finalize` stages instead of
constructing independent filtered chromosome VCFs. See
[WES/WGS Cohort Joint Genotyping](../pipelines/wes-wgs-cohort#execution-modes).

## 4. Inspect the outputs

```text
cbicall_bash_gatk-4.6_wes_cohort_b37_*/
  01_genomicsdb/
  02_varcall/
  logs/
  log.json
  run-report.json
  run-report.html
```

The final joint callset is `02_varcall/cohort.gv.QC.vcf.gz`. The raw
joint-genotyped VCF and conditional VQSR outputs are also retained in
`02_varcall/`.

</TabItem>
</Tabs>

See [Outputs](../help/outputs) for the complete file reference and
[Configuration Reference](../help/configuration-reference) for all accepted
YAML keys.
