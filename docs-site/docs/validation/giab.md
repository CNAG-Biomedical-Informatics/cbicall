# GIAB Benchmarking

Genome in a Bottle (GIAB) benchmarking checks the analytical behavior of a
selected variant-calling workflow against independent truth sets. In CBIcall,
this is different from run-comparison validation: GIAB asks whether the calls
agree with a reference benchmark, while CBIcall reports ask whether two
executions used the same workflow, resources, runtime context, and outputs.

:::important[Scope]
This page documents a WGS benchmark for the native CBIcall GATK 4.6 pipeline on
GRCh38. It is a workflow-level analytical benchmark, not a validation of the
CBIcall Python wrapper itself.
:::

## Benchmark Design

The benchmark uses the GIAB Ashkenazim trio:

| Sample | Relationship | Query VCFs | GIAB truth resources |
| --- | --- | --- | --- |
| HG002 / NA24385 | Son | Single-sample WGS VCF and extracted joint-genotyped VCF | `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` and `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |
| HG003 / NA24149 | Father | Single-sample WGS VCF and extracted joint-genotyped VCF | `HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` and `HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |
| HG004 / NA24143 | Mother | Single-sample WGS VCF and extracted joint-genotyped VCF | `HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` and `HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |

:::info[Flow]
Download GIAB FASTQs -> check mate pairs -> downsample to approximately 30x ->
run CBIcall WGS single-sample mode for each trio member -> run CBIcall WGS
cohort mode on the three gVCFs -> split `cohort.gv.QC.vcf.gz` by sample ->
compare each query VCF with hap.py.
:::

:::warning[No joint trio truth VCF]
GIAB small-variant truth sets are sample-specific. A joint-genotyped trio VCF is
therefore benchmarked by extracting HG002, HG003, and HG004 and comparing each
sample separately against its own truth VCF and confident-region BED. Mendelian
consistency is a separate trio-level QC.
:::

## 1. Prepare FASTQs

Download the official GIAB Ashkenazim trio FASTQs from NCBI:

| Sample | NCBI GIAB FTP source | Local directory |
| --- | --- | --- |
| HG002 | `ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/` | `AshkenazimTrio_fastq/HG002_NA24385_son/` |
| HG003 | `ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/HG003_HiSeq300x_fastq/` | `AshkenazimTrio_fastq/HG003_NA24149_father/` |
| HG004 | `ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/HG004_HiSeq300x_fastq/` | `AshkenazimTrio_fastq/HG004_NA24143_mother/` |

Use an `aria2` manifest with one record per FASTQ. The important part is to keep
the source URL, output name, destination directory, and checksum together:

<details className="cbicallCodeDetails">
<summary>aria2 manifest record</summary>
<div>

```text
https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/<sample>/<project>/<file>.fastq.gz
  dir=/path/to/GIAB/AshkenazimTrio_fastq/<sample>
  out=<unique-local-file-name>.fastq.gz
  checksum=md5=<expected-md5>
```

</div>
</details>

Run the download as a resumable job:

<details className="cbicallCodeDetails">
<summary>Download FASTQs</summary>
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
  -i /path/to/giab_trio_aria2.txt \
  > /path/to/giab_trio_aria2.log 2>&1 &
```

</div>
</details>

Check that every downloaded FASTQ has its mate:

<details className="cbicallCodeDetails">
<summary>Mate-pair check</summary>
<div>

```bash
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

No output means that all detected `R1`/`R2` files have their expected mate.

:::tip[Downsampling tools]
`seqtk` performs reproducible FASTQ subsampling; `pigz` is a parallel gzip
implementation used to compress the downsampled FASTQs.
:::

Downsample each sample to approximately 30x:

<details className="cbicallCodeDetails">
<summary>Downsample one sample</summary>
<div>

```bash
SAMPLE=HG002
SAMPLING_FRACTION=0.10

ls -1 *_R1_*.fastq.gz | sort > r1.list
ls -1 *_R2_*.fastq.gz | sort > r2.list

cat $(<r1.list) \
  | seqtk sample -s42 - "$SAMPLING_FRACTION" \
  | pigz -p 16 > "${SAMPLE}_30x_R1_001.fastq.gz"

cat $(<r2.list) \
  | seqtk sample -s42 - "$SAMPLING_FRACTION" \
  | pigz -p 16 > "${SAMPLE}_30x_R2_001.fastq.gz"
```

</div>
</details>

## 2. Download Truth Sets

Download the GIAB v4.2.1 GRCh38 truth VCFs, VCF indexes, and confident-region
BED files:

<details className="cbicallCodeDetails">
<summary>Download GIAB truth resources</summary>
<div>

```bash
TRUTH_DIR=/path/to/GIAB/truthsets
mkdir -p "$TRUTH_DIR"

wget -c -P "$TRUTH_DIR" \
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

<details className="cbicallCodeDetails">
<summary>Expected truth-set files</summary>
<div>

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

</div>
</details>

## 3. Run CBIcall

Run WGS single-sample mode once per trio member:

<details className="cbicallCodeDetails">
<summary>Single-sample WGS YAML</summary>
<div>

```yaml
mode: single
pipeline: wgs
workflow_backend: bash
software_stack: gatk-4.6
genome: hg38
cleanup_bam: true
input_dir: /path/to/AshkenazimTrio_fastq/HG002_NA24385_son/
```

</div>
</details>

<details className="cbicallCodeDetails">
<summary>Run single-sample WGS</summary>
<div>

```bash
bin/cbicall run -p HG002_wgs_single_hg38.yaml -t 12
```

</div>
</details>

Each single-sample run produces:

<details className="cbicallCodeDetails">
<summary>Single-sample outputs</summary>
<div>

```text
02_varcall/HG002.hc.QC.vcf.gz
02_varcall/HG002.hc.g.vcf.gz
```

</div>
</details>

Create a sample map for joint genotyping:

<details className="cbicallCodeDetails">
<summary>Joint-genotyping sample map</summary>
<div>

```text
sample	gvcf
HG002	/path/to/HG002.hc.g.vcf.gz
HG003	/path/to/HG003.hc.g.vcf.gz
HG004	/path/to/HG004.hc.g.vcf.gz
```

</div>
</details>

Run WGS cohort mode:

<details className="cbicallCodeDetails">
<summary>WGS cohort YAML</summary>
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
<summary>Run WGS cohort mode</summary>
<div>

```bash
bin/cbicall run -p giab_trio_wgs_hg38.yaml -t 4
```

</div>
</details>

The final joint-genotyped VCF is:

<details className="cbicallCodeDetails">
<summary>Joint-genotyped output</summary>
<div>

```text
02_varcall/cohort.gv.QC.vcf.gz
```

</div>
</details>

:::success[Keep the audit files]
Store `log.json`, `cbicall-execution-contract.json`, `run-report.json`,
`run-report.html`, and the workflow log with the GIAB results. These files
record the CBIcall configuration, selected workflow, resource bundle, runtime
context, and output fingerprints.
:::

## 4. Split the Joint VCF

Extract one query VCF per trio member from `cohort.gv.QC.vcf.gz`:

<details className="cbicallCodeDetails">
<summary>Split the joint VCF by sample</summary>
<div>

```bash
bcftools view -s HG002 -Oz -o HG002.joint.QC.vcf.gz cohort.gv.QC.vcf.gz
bcftools view -s HG003 -Oz -o HG003.joint.QC.vcf.gz cohort.gv.QC.vcf.gz
bcftools view -s HG004 -Oz -o HG004.joint.QC.vcf.gz cohort.gv.QC.vcf.gz

bcftools index -t HG002.joint.QC.vcf.gz
bcftools index -t HG003.joint.QC.vcf.gz
bcftools index -t HG004.joint.QC.vcf.gz
```

</div>
</details>

## 5. Prepare hap.py

The prebuilt `pkrusche/hap.py:latest` image may fail on modern Docker versions
because it uses a deprecated image manifest. Build hap.py locally instead:

<details className="cbicallCodeDetails">
<summary>Build and test hap.py</summary>
<div>

```bash
cd /path/to/hap.py
docker build -t local/hap.py .
docker run --rm local/hap.py /opt/hap.py/bin/hap.py --version
```

</div>
</details>

Expected version:

<details className="cbicallCodeDetails">
<summary>Expected hap.py version</summary>
<div>

```text
Hap.py v0.3.15
```

</div>
</details>

For `--engine=vcfeval`, RTG also needs an SDF reference created from the exact
same FASTA passed to `hap.py -r`:

<details className="cbicallCodeDetails">
<summary>Create RTG SDF reference</summary>
<div>

```bash
REF_DIR=/path/to/cbicall-data/Databases/GATK_bundle/hg38

docker run --rm \
  -v "$REF_DIR":/ref \
  local/hap.py \
  /opt/hap.py/libexec/rtg-tools-install/rtg format \
  -o /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
  /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
```

</div>
</details>

:::caution[Reference consistency]
The query VCF, truth VCF, confident BED, FASTA, and RTG SDF must use matching
coordinates and contig naming. This benchmark uses the Broad/GATK GRCh38 FASTA
with UCSC-style contigs such as `chr1`.
:::

## 6. Run hap.py

For a single-sample CBIcall WGS run:

<details className="cbicallCodeDetails">
<summary>hap.py for single-sample WGS</summary>
<div>

```bash
docker run --rm \
  -v /path/to/cbicall-run/02_varcall:/query \
  -v /path/to/cbicall-data/Databases/GATK_bundle/hg38:/ref \
  -v /path/to/GIAB/truthsets:/truth \
  -v /path/to/GIAB/happy_results:/out \
  local/hap.py \
  /opt/hap.py/bin/hap.py \
  /truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /query/HG002.hc.QC.vcf.gz \
  -f /truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
  --engine=vcfeval \
  --engine-vcfeval-template /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
  -o /out/HG002.single_wgs.happy
```

</div>
</details>

For a sample extracted from the joint-genotyped trio VCF:

<details className="cbicallCodeDetails">
<summary>hap.py for joint-genotyped sample</summary>
<div>

```bash
docker run --rm \
  -v /path/to/extracted-joint-vcfs:/query \
  -v /path/to/cbicall-data/Databases/GATK_bundle/hg38:/ref \
  -v /path/to/GIAB/truthsets:/truth \
  -v /path/to/GIAB/happy_results:/out \
  local/hap.py \
  /opt/hap.py/bin/hap.py \
  /truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /query/HG002.joint.QC.vcf.gz \
  -f /truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
  --engine=vcfeval \
  --engine-vcfeval-template /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
  -o /out/HG002.joint_wgs.happy
```

</div>
</details>

Expected output files:

<details className="cbicallCodeDetails">
<summary>hap.py outputs</summary>
<div>

```text
HG002.single_wgs.happy.summary.csv
HG002.single_wgs.happy.extended.csv
HG002.single_wgs.happy.metrics.json
HG002.single_wgs.happy.roc.*.csv.gz
HG002.single_wgs.happy.vcf.gz
HG002.single_wgs.happy.vcf.gz.tbi
```

</div>
</details>

:::tip[Comparison engine]
Use `--engine=vcfeval` for GIAB/WGS benchmarking. The default `xcmp` engine is
better reserved for quick checks or reproducing older hap.py behavior.
:::

## Joint Genotyping Benchmark Results

The final GIAB report will summarize only the three sample-level comparisons
extracted from the joint-genotyped trio VCF. Single-sample WGS runs are used to
generate gVCFs for joint genotyping, but they are not reported as final GIAB
results.

<div className="giabResultsHero">
  <div>
    <span className="giabResultsLabel">Final reporting target</span>
    <h3>Trio joint genotyping, evaluated per sample</h3>
    <p>
      Query source: <code>cohort.gv.QC.vcf.gz</code>, split into HG002, HG003, and HG004<br />
      Truth: GIAB v4.2.1 GRCh38 per-sample truth VCFs<br />
      Confident regions: GIAB per-sample <code>*_benchmark_noinconsistent.bed</code> files
    </p>
  </div>
  <div className="giabResultsStatus">layout example</div>
</div>

<div className="giabHeatmap">
  <div className="giabHeatmapHeader giabHeatmapCorner">Sample / filter</div>
  <div className="giabHeatmapHeader">Joint-genotyped SNP</div>
  <div className="giabHeatmapHeader">Joint-genotyped INDEL</div>
  <div className="giabHeatmapSample giabHeatmapHg002">HG002 <span>PASS</span></div>
  <div className="giabHeatmapHg002 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 95.4%</span>
    <span><b>R</b> 91.3%</span>
    <span><b>P</b> 99.9%</span>
  </div>
  <div className="giabHeatmapHg002 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 96.1%</span>
    <span><b>R</b> 93.4%</span>
    <span><b>P</b> 99.0%</span>
  </div>
  <div className="giabHeatmapSample giabHeatmapHg002">HG002 <span>ALL</span></div>
  <div className="giabHeatmapHg002 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 99.1%</span>
    <span><b>R</b> 99.2%</span>
    <span><b>P</b> 99.0%</span>
  </div>
  <div className="giabHeatmapHg002 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 98.5%</span>
    <span><b>R</b> 98.3%</span>
    <span><b>P</b> 98.7%</span>
  </div>
  <div className="giabHeatmapSample giabHeatmapHg003">HG003 <span>PASS</span></div>
  <div className="giabHeatmapHg003 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 95.1%</span>
    <span><b>R</b> 91.0%</span>
    <span><b>P</b> 99.6%</span>
  </div>
  <div className="giabHeatmapHg003 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 95.8%</span>
    <span><b>R</b> 92.9%</span>
    <span><b>P</b> 98.9%</span>
  </div>
  <div className="giabHeatmapSample giabHeatmapHg003">HG003 <span>ALL</span></div>
  <div className="giabHeatmapHg003 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 99.0%</span>
    <span><b>R</b> 99.1%</span>
    <span><b>P</b> 98.9%</span>
  </div>
  <div className="giabHeatmapHg003 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 98.4%</span>
    <span><b>R</b> 98.2%</span>
    <span><b>P</b> 98.6%</span>
  </div>
  <div className="giabHeatmapSample giabHeatmapHg004">HG004 <span>PASS</span></div>
  <div className="giabHeatmapHg004 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 95.3%</span>
    <span><b>R</b> 91.4%</span>
    <span><b>P</b> 99.7%</span>
  </div>
  <div className="giabHeatmapHg004 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 96.0%</span>
    <span><b>R</b> 93.2%</span>
    <span><b>P</b> 99.0%</span>
  </div>
  <div className="giabHeatmapSample giabHeatmapHg004">HG004 <span>ALL</span></div>
  <div className="giabHeatmapHg004 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 99.1%</span>
    <span><b>R</b> 99.2%</span>
    <span><b>P</b> 99.0%</span>
  </div>
  <div className="giabHeatmapHg004 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 98.5%</span>
    <span><b>R</b> 98.3%</span>
    <span><b>P</b> 98.7%</span>
  </div>
</div>

:::note[Illustrative values]
The heatmap values above are temporary placeholders to show the final reporting
layout. They must be replaced with hap.py recall, precision, and F1 values from
the sample VCFs extracted from `cohort.gv.QC.vcf.gz`. PASS rows represent the
final hard-filtered call set; ALL rows show the same comparison before
restricting to PASS variants.
:::


## Exact Benchmark Commands

The public protocol above uses portable paths. The exact local paths below are
kept only to make the manuscript benchmark traceable.

<details className="cbicallCodeDetails">
<summary>Exact local commands used for the manuscript benchmark</summary>
<div>

Example `aria2` records:

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

FASTQ download:

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

HG002 downsampling:

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

HG002 CBIcall YAML:

```yaml
mode: single
pipeline: wgs
workflow_backend: bash
software_stack: gatk-4.6
cleanup_bam: true
genome: hg38
input_dir: AshkenazimTrio_fastq/HG002_NA24385_son/
```

Joint genotyping `cohort.yaml`:

```yaml
mode: cohort
pipeline: wgs
workflow_backend: bash
software_stack: gatk-4.6
genome: hg38
sample_map: ./sample_map.tsv
```

Joint genotyping `sample_map.tsv`:

```text
HG002	/media/mrueda/4TBB/GIAB/AshkenazimTrio_fastq/HG002_NA24385_son/cbicall_bash_gatk-4.6_wgs_single_hg38_178120808914892/02_varcall/HG002.hc.g.vcf.gz
HG003	/media/mrueda/4TBB/GIAB/AshkenazimTrio_fastq/HG003_NA24149_father/cbicall_bash_gatk-4.6_wgs_single_hg38_178128498815686/02_varcall/HG003.hc.g.vcf.gz
HG004	/media/mrueda/4TBB/GIAB/AshkenazimTrio_fastq/HG004_NA24143_mother/cbicall_bash_gatk-4.6_wgs_single_hg38_178128547716032/02_varcall/HG004.hc.g.vcf.gz
```

Local hap.py image:

```bash
cd /media/mrueda/4TBA/GIAB/hap.py
docker build -t local/hap.py .
```

RTG SDF:

```bash
docker run --rm \
  -v /media/mrueda/4TBB/cbicall-data/Databases/GATK_bundle/hg38:/ref \
  local/hap.py \
  /opt/hap.py/libexec/rtg-tools-install/rtg format \
  -o /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
  /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
```

HG002 single-WGS hap.py comparison:

```bash
docker run --rm \
  -v /media/mrueda/4TBB/GIAB/AshkenazimTrio_fastq/HG002_NA24385_son/cbicall_bash_gatk-4.6_wgs_single_hg38_178120808914892/02_varcall:/query \
  -v /media/mrueda/4TBB/cbicall-data/Databases/GATK_bundle/hg38:/ref \
  -v /media/mrueda/4TBA/GIAB/truthsets:/truth \
  -v /media/mrueda/4TBA/GIAB/happy_results:/out \
  local/hap.py \
  /opt/hap.py/bin/hap.py \
  /truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
  /query/HG002.hc.QC.vcf.gz \
  -f /truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
  -r /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
  --engine=vcfeval \
  --engine-vcfeval-template /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
  -o /out/HG002.single_wgs.happy
```

</div>
</details>
