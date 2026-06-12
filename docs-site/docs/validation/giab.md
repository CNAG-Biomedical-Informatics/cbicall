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

The GIAB validation flow is one WGS joint-genotyping benchmark evaluated as
three sample-level truth comparisons. CBIcall produces one joint-genotyped trio
VCF for HG002, HG003, and HG004; each sample is then extracted from that joint
VCF and compared against that sample's own GIAB v4.2.1 GRCh38 truth VCF and
confident BED.

:::note[No joint trio truth VCF]
For GIAB small-variant benchmarking, the standard comparison is sample-level:
query VCF versus truth VCF inside confident regions. This page does not use a
single "joint trio truth VCF". The joint-called trio VCF is evaluated by
extracting HG002, HG003, and HG004 and benchmarking the three sample VCFs
separately; Mendelian consistency is a separate trio-level QC.
:::

| Sample | Relationship | CBIcall query | GIAB truth resources |
| --- | --- | --- | --- |
| HG002 / NA24385 | Son | HG002 extracted from the CBIcall joint-genotyped trio VCF | `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` and `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |
| HG003 / NA24149 | Father | HG003 extracted from the CBIcall joint-genotyped trio VCF | `HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` and `HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |
| HG004 / NA24143 | Mother | HG004 extracted from the CBIcall joint-genotyped trio VCF | `HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` and `HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |

## Shared FASTQ Preparation

The benchmark starts from the official GIAB Ashkenazim trio FASTQs hosted by
NCBI under `ReferenceSamples/giab/data/AshkenazimTrio/`. The local benchmark
keeps one directory per trio member:

| GIAB sample | Relationship | NCBI GIAB FTP source | Local FASTQ directory |
| --- | --- | --- | --- |
| HG002 / NA24385 | Son | `HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/` | `AshkenazimTrio_fastq/HG002_NA24385_son/` |
| HG003 / NA24149 | Father | `HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003_HiSeq300x_fastq/` | `AshkenazimTrio_fastq/HG003_NA24149_father/` |
| HG004 / NA24143 | Mother | `HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq/` | `AshkenazimTrio_fastq/HG004_NA24143_mother/` |

For this run, the selected FASTQ URLs were stored in a local `aria2` input file
and downloaded as a resumable background job. Any equivalent download method is
acceptable if it records the same source URLs and preserves the sample-specific
directory structure.

<details className="cbicallCodeDetails">
<summary>Example aria2 input records</summary>
<div>

```text
https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R1_001.fastq.gz
  dir=/media/mrueda/4TBA/GIAB/AshkenazimTrio_fastq/HG002_NA24385_son
  out=140528_D00360_0018_AH8VC6ADXX__Project_RM8391_RM8392__Sample_2A1__2A1_CGATGT_L001_R1_001.fastq.gz
  checksum=md5=c2ae5e412fb211974f9a9a46a5392428

https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/Project_RM8391_RM8392/Sample_2A1/2A1_CGATGT_L001_R2_001.fastq.gz
  dir=/media/mrueda/4TBA/GIAB/AshkenazimTrio_fastq/HG002_NA24385_son
  out=140528_D00360_0018_AH8VC6ADXX__Project_RM8391_RM8392__Sample_2A1__2A1_CGATGT_L001_R2_001.fastq.gz
  checksum=md5=83826a956fc90c501645391314b2abf3
```

</div>
</details>

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

## Truth Set Preparation

Download the per-sample GIAB v4.2.1 GRCh38 truth VCFs, VCF indexes, and
confident-region BED files for HG002, HG003, and HG004:

<details className="cbicallCodeDetails">
<summary>Download GIAB truth sets</summary>
<div>

```bash
mkdir -p /media/mrueda/4TBA/GIAB/truthsets

wget -c -P /media/mrueda/4TBA/GIAB/truthsets \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

</div>
</details>

Expected files:

```text
HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

## Benchmark Workflow

This benchmark compares the **native CBIcall WGS joint-genotyping pipeline using
GATK 4.6** on GRCh38 against the matching GIAB truth sets. The joint VCF is
split into three sample VCFs, and each sample comparison uses `hap.py` with the
`vcfeval` engine.

| Field | Value |
| --- | --- |
| Samples | Three Ashkenazim trio samples extracted from one joint-genotyped VCF |
| Upstream sequencing data | Downsampled GIAB high-depth WGS FASTQs |
| CBIcall input | Three gVCFs generated from the downsampled trio FASTQs |
| CBIcall workflow | Native WGS joint-genotyping pipeline, `workflow_backend: bash`, `software_stack: gatk-4.6`, `genome: hg38` |
| Truth VCF/BED | Matching GIAB GRCh38 benchmark VCF and BED per sample |
| Reference | Broad GRCh38 resource FASTA: `resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta` |
| Comparison tool | Locally built `local/hap.py` Docker image, `hap.py` with `--engine vcfeval` |
| Output summary | `results_<sample>_happy.summary.csv` |

### Run CBIcall

Create a sample map that points to the three trio gVCFs generated from the same
downsampled FASTQs:

<details className="cbicallCodeDetails">
<summary>Trio gVCF sample map</summary>
<div>

```text
sample	gvcf
HG002	/path/to/HG002.hc.g.vcf.gz
HG003	/path/to/HG003.hc.g.vcf.gz
HG004	/path/to/HG004.hc.g.vcf.gz
```

</div>
</details>

Then run the WGS joint-genotyping workflow:

<details className="cbicallCodeDetails">
<summary>Joint-genotyping CBIcall YAML</summary>
<div>

```yaml
mode: cohort
pipeline: wgs
workflow_backend: bash
software_stack: gatk-4.6
genome: hg38
sample_map: /path/to/giab_trio_sample_map.tsv
```

</div>
</details>

<details className="cbicallCodeDetails">
<summary>Joint-genotyping CBIcall command</summary>
<div>

```bash
bin/cbicall run -p giab_trio_wgs_hg38.yaml -t 4
```

</div>
</details>

Extract one query VCF per trio member from the CBIcall joint VCF:

<details className="cbicallCodeDetails">
<summary>Extract sample VCFs from the joint VCF</summary>
<div>

```bash
bcftools view -s HG002 -Oz -o HG002.joint.vcf.gz trio.joint.vcf.gz
bcftools view -s HG003 -Oz -o HG003.joint.vcf.gz trio.joint.vcf.gz
bcftools view -s HG004 -Oz -o HG004.joint.vcf.gz trio.joint.vcf.gz

tabix -p vcf HG002.joint.vcf.gz
tabix -p vcf HG003.joint.vcf.gz
tabix -p vcf HG004.joint.vcf.gz
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

### Prepare hap.py Docker Image

The prebuilt `pkrusche/hap.py:latest` image could not be pulled because Docker
rejected its deprecated image manifest format. For this benchmark, `hap.py` was
therefore built locally from the checked-out `hap.py` source tree:

<details className="cbicallCodeDetails">
<summary>Build the local hap.py image</summary>
<div>

```bash
cd /media/mrueda/4TBA/GIAB/hap.py
docker build -t local/hap.py .
```

</div>
</details>

The local image name is `local/hap.py`. The build also ran
`bin/test_haplotypes`, which passed without errors.

<details className="cbicallCodeDetails">
<summary>Smoke-test the local hap.py image</summary>
<div>

```bash
docker run --rm local/hap.py /opt/hap.py/bin/hap.py --version
```

Expected version:

```text
Hap.py v0.3.15
```

</div>
</details>

### Run hap.py

Run each comparison from the directory that contains the extracted sample VCFs:

<details className="cbicallCodeDetails">
<summary>Enter the final-VCF directory</summary>
<div>

```bash
cd /path/to/giab_joint_benchmark/
```

</div>
</details>

The hg38 reference directory used in this benchmark is:

```text
/media/mrueda/4TBB/cbicall-data/Databases/GATK_bundle/hg38
```

The FASTA passed to `hap.py -r` is the Broad/GATK hg38 FASTA:

```text
/media/mrueda/4TBB/cbicall-data/Databases/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
```

The reference already had the required `.fai` and `.dict` companion files:

```text
resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai
resources_broad_hg38_v0_Homo_sapiens_assembly38.dict
```

For `--engine=vcfeval`, RTG requires an SDF reference created from the **exact
same FASTA** passed to `hap.py -r`. Create this once:

<details className="cbicallCodeDetails">
<summary>Create the RTG reference SDF</summary>
<div>

```bash
docker run --rm \
  -v /media/mrueda/4TBB/cbicall-data/Databases/GATK_bundle/hg38:/ref \
  local/hap.py \
  /opt/hap.py/libexec/rtg-tools-install/rtg format \
  -o /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
  /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
```

RTG reported:

```text
Detected: Human GRCh38 with UCSC naming
Number of sequences: 3366
Total residues: 3217346917
```

</div>
</details>

Then mount both the dataset and the reference directory when running `hap.py`:

<details className="cbicallCodeDetails">
<summary>hap.py comparison command</summary>
<div>

```bash
docker run --rm \
  -v /path/to/HG002:/data \
  -v /media/mrueda/4TBB/cbicall-data/Databases/GATK_bundle/hg38:/ref \
  -v /media/mrueda/4TBA/GIAB/truthsets:/truth \
  local/hap.py \
  /opt/hap.py/bin/hap.py \
  /truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /data/HG002.joint.vcf.gz \
  -f /truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
  --engine=vcfeval \
  --engine-vcfeval-template /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
  -o /data/HG002.happy
```

</div>
</details>

Expected outputs are:

```text
HG002.happy.summary.csv
HG002.happy.extended.csv
HG002.happy.metrics.json
HG002.happy.roc.*.csv.gz
HG002.happy.vcf.gz
HG002.happy.vcf.gz.tbi
```

:::tip[Engine and reference consistency]
Use `--engine=vcfeval` for GIAB/WGS benchmarking. The default `xcmp` engine is
better reserved for quick checks or reproducing older `hap.py` behavior. The
query VCF, truth VCF, confident BED, reference FASTA, and SDF must all use
matching coordinates and contig naming; here the reference is Broad/GATK hg38
with UCSC-style contigs such as `chr1`.
:::
