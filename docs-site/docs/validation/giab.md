import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# GIAB Benchmarking

Genome in a Bottle (GIAB) benchmarking is analytical validation of the selected
variant-calling workflow, not validation of CBIcall itself. Use this page to
record truth-set comparisons for native CBIcall WES/WGS workflows.

:::important[Scope]
GIAB benchmarking evaluates the analytical behavior of the selected workflow and
resource bundle against a truth set. It complements, but does not replace,
CBIcall run comparison, which audits whether executions used the same workflow,
resources, runtime context, and final-output fingerprints.
:::

## Benchmark Plan

The GIAB validation flow is organized into two WGS benchmarks:

| Benchmark | Purpose | Samples |
| --- | --- | --- |
| Single-sample WGS on GRCh38 | Compare each native CBIcall single-sample WGS callset against its matching GIAB truth set. | Three Ashkenazim trio samples, reported independently. |
| Trio WGS joint genotyping | Evaluate the native CBIcall WGS cohort workflow on the trio after single-sample gVCF generation. | Father, mother, and son jointly genotyped as a trio. |

## Shared FASTQ Preparation

The benchmark starts from the GIAB Ashkenazim trio FASTQs. The download was
driven by a plain-text aria2 input manifest:

```text
giab_AJtrio_persondir_aria2.txt
```

Each line in that manifest points to one FASTQ URL from the GIAB Ashkenazim trio
directory. The download was run as a resumable background job:

<details className="cbicallCodeDetails">
<summary>FASTQ download command</summary>
<div>

```bash
aria2c -c \
  -x 1 \
  -s 1 \
  -j 2 \
  --retry-wait=120 \
  --max-tries=0 \
  --timeout=120 \
  --connect-timeout=30 \
  --auto-file-renaming=false \
  --allow-overwrite=false \
  -i ./giab_AJtrio_persondir_aria2.txt \
  > giab_AJtrio_persondir_aria2.log 2>&1 &
```

</div>
</details>

Before downsampling, verify that every downloaded FASTQ has its expected mate:

<details className="cbicallCodeDetails">
<summary>FASTQ mate-pair check</summary>
<div>

```bash
#!/usr/bin/env bash

find . -type f -name "*.fastq.gz" | while read -r f; do
  if [[ $f == *_R1_* ]]; then
    mate="${f/_R1_/_R2_}"
  elif [[ $f == *_R2_* ]]; then
    mate="${f/_R2_/_R1_}"
  else
    continue
  fi

  [[ -f "$mate" ]] || printf "MISSING\t%s\nEXPECTED\t%s\n\n" "$f" "$mate"
done
```

</div>
</details>

No output means all detected `R1`/`R2` files have their corresponding mate. The
source dataset is high-depth WGS; the benchmark then downsamples the trio FASTQs
before running CBIcall.

For each sample, merge and downsample the FASTQ chunks to target approximately
30x coverage. The example below shows the HG002 downsampling command:

:::note[Downsampling tools]
`seqtk` is a lightweight toolkit for FASTA/FASTQ processing; here it performs
reproducible random read sampling with seed `42`. `pigz` is a parallel gzip
implementation used to compress the downsampled FASTQ output faster than
single-threaded `gzip`.
:::

<details className="cbicallCodeDetails">
<summary>Downsample FASTQs to approximately 30x</summary>
<div>

```bash
ls -1 *_R1_*.fastq.gz | sort > r1.list
ls -1 *_R2_*.fastq.gz | sort > r2.list

cat $(<r1.list) \
  | seqtk sample -s42 - 0.10 \
  | pigz -p 16 > HG002_30x_R1_001.fastq.gz

cat $(<r2.list) \
  | seqtk sample -s42 - 0.10 \
  | pigz -p 16 > HG002_30x_R2_001.fastq.gz
```

</div>
</details>

## Benchmarks

<Tabs groupId="giab-benchmark">
<TabItem value="single-wgs" label="Single-Sample WGS" default>

This benchmark compares each **native CBIcall WGS single-sample pipeline using
GATK 4.6** on GRCh38 against its matching GIAB truth set. The comparison uses
`hap.py` with the `vcfeval` engine.

| Field | Value |
| --- | --- |
| Samples | Three Ashkenazim trio samples, reported independently |
| Sequencing data | Downsampled GIAB high-depth WGS FASTQs |
| CBIcall workflow | Native WGS single-sample pipeline, `workflow_backend: bash`, `software_stack: gatk-4.6`, `genome: hg38` |
| Truth VCF/BED | Matching GIAB GRCh38 benchmark VCF and BED per sample |
| Reference | Broad GRCh38 resource FASTA: `resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta` |
| Comparison tool | `hap.py` with `--engine vcfeval` |
| Output summary | `results_<sample>_happy.summary.csv` |

### Run CBIcall

Create one WGS single-sample parameters YAML per downsampled trio sample. This
example is for the son sample, HG002 / NA24385:

<details className="cbicallCodeDetails">
<summary>Single-sample CBIcall YAML</summary>
<div>

```yaml
mode: single
pipeline: wgs
workflow_backend: bash
software_stack: gatk-4.6
cleanup_bam: true
genome: hg38
input_dir: AshkenazimTrio_fastq/HG002_NA24385_son/
```

</div>
</details>

For the general execution pattern, see [General Usage](../usage). Then run CBIcall as usual:

<details className="cbicallCodeDetails">
<summary>Single-sample CBIcall command</summary>
<div>

```bash
bin/cbicall run -p <sample>_wgs_hg38.yaml -t 4
```

</div>
</details>

Keep the generated CBIcall audit files with the benchmark output:

```text
log.json
cbicall-execution-contract.json
run-report.json
run-report.html
<workflow>.log
```

### Prepare hap.py Environment

On the workstation used for this run, `hap.py` was loaded through a shell
bootstrap script:

<details className="cbicallCodeDetails">
<summary>hap.py module environment</summary>
<div>

```bash
source ~/load_happy.sh
```

The script loaded these modules:

```bash
module load GCC/10.2.0
module load OpenSSL/1.1
module load Python/2.7.18-GCCcore-10.2.0
module load Java
```

</div>
</details>

### Run hap.py

Run the comparison from the CBIcall final-VCF directory:

<details className="cbicallCodeDetails">
<summary>Enter the final-VCF directory</summary>
<div>

```bash
cd /path/to/<sample>_cbicall_bash_wgs_single_hg38_gatk-4.6_<run_id>/02_varcall/
```

</div>
</details>

Create the `vcfeval` reference SDF once from the same GRCh38 FASTA used for the
benchmark:

<details className="cbicallCodeDetails">
<summary>Create the RTG reference SDF</summary>
<div>

```bash
rtg format \
  -o ../../../benchmark/resources_broad_hg38_v0_Homo_sapiens_assembly38 \
  ../../../benchmark/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
```

</div>
</details>

Then run `hap.py`:

<details className="cbicallCodeDetails">
<summary>hap.py comparison command</summary>
<div>

```bash
hap.py \
  ../../../benchmark/<sample>_GRCh38_1_22_benchmark.vcf.gz \
  <sample>.hc.QC.vcf.gz \
  -f ../../../benchmark/<sample>_GRCh38_1_22_benchmark.bed \
  -r ../../../benchmark/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
  -o results_<sample>_happy \
  --engine vcfeval \
  --engine-vcfeval-template ../../../benchmark/resources_broad_hg38_v0_Homo_sapiens_assembly38 \
  --threads 8 \
  > results_<sample>_happy.log 2>&1
```

</div>
</details>

The main summary table is:

```text
results_<sample>_happy.summary.csv
```

</TabItem>
<TabItem value="trio-wgs" label="Trio Joint Genotyping">

This benchmark evaluates the **native CBIcall WGS cohort pipeline using GATK
4.6** on the three downsampled Ashkenazim trio samples after single-sample gVCF
generation.

| Field | Value |
| --- | --- |
| Samples | Father, mother, and son from the Ashkenazim trio |
| Input to cohort step | Single-sample gVCFs generated by CBIcall from the same downsampled FASTQs |
| CBIcall workflow | Native WGS cohort pipeline, `workflow_backend: bash`, `software_stack: gatk-4.6`, `genome: hg38` |
| Benchmark outputs | Joint-genotyped trio VCF and CBIcall audit files |
| Planned checks | Truth-set comparison per sample and trio-level consistency checks |

Create a cohort sample map from the three single-sample gVCFs:

<details className="cbicallCodeDetails">
<summary>Trio cohort sample map</summary>
<div>

```text
sample	gvcf
<father>	/path/to/<father>.hc.g.vcf.gz
<mother>	/path/to/<mother>.hc.g.vcf.gz
<son>	/path/to/<son>.hc.g.vcf.gz
```

</div>
</details>

Then run the WGS cohort workflow:

<details className="cbicallCodeDetails">
<summary>Trio cohort CBIcall YAML</summary>
<div>

```yaml
mode: cohort
pipeline: wgs
workflow_backend: bash
software_stack: gatk-4.6
resource: "cbicall-germline-resources-v1"
genome: hg38
sample_map: /path/to/giab_trio_sample_map.tsv
```

</div>
</details>

<details className="cbicallCodeDetails">
<summary>Trio cohort CBIcall command</summary>
<div>

```bash
bin/cbicall run -p giab_trio_wgs_hg38.yaml -t 4
```

</div>
</details>

Keep the same CBIcall audit files as for the single-sample benchmark, plus the
cohort final VCF and any downstream truth-set or trio-consistency summaries.

</TabItem>
</Tabs>
