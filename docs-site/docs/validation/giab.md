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

| Sample | Relationship | Benchmark role |
| --- | --- | --- |
| HG002 / NA24385 | Son | Child sample and per-sample GIAB truth comparison |
| HG003 / NA24149 | Father | Parent sample and per-sample GIAB truth comparison |
| HG004 / NA24143 | Mother | Parent sample and per-sample GIAB truth comparison |

:::info[Flow]
Download GIAB FASTQs -> check mate pairs -> downsample to approximately 30x ->
run CBIcall WGS single-sample mode for each trio member -> run CBIcall WGS
cohort mode on the three gVCFs -> extract `cohort.gv.QC.vcf.gz` by sample ->
compare each query VCF with hap.py -> check trio-level Mendelian consistency
on the joint-genotyped cohort VCF.
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
cbicall run -p HG002_wgs_single_hg38.yaml -t 12
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
cbicall run -p giab_trio_wgs_hg38.yaml -t 4
```

</div>
</details>

:::success[Keep the audit files]
Store `log.json`, `cbicall-execution-contract.json`, `run-report.json`,
`run-report.html`, and the workflow log with the GIAB results. These files
record the CBIcall configuration, selected workflow, resource bundle, runtime
context, and output fingerprints.
:::

## 4. Extract Sample VCFs

Extract one query VCF per trio member from `cohort.gv.QC.vcf.gz`. The `-c 1`
filter is applied after sample subsetting so the per-sample VCF keeps only
sites where that sample carries at least one non-reference allele. `FORMAT/AD`
is removed from the benchmark query VCFs because hap.py 0.3.15 can crash while
preprocessing decomposed multiallelic records with allele-depth vectors, and
the hap.py comparison does not use `AD`:

<details className="cbicallCodeDetails">
<summary>Extract sample VCFs from the joint VCF</summary>
<div>

```bash
COHORT_VCF=../../../4TBB/GIAB/cbicall_bash_gatk-4.6_wgs_cohort_hg38_178167960428471/02_varcall/cohort.gv.QC.vcf.gz

$BCFTOOLS view -s HG002 -c 1 "$COHORT_VCF" \
  | $BCFTOOLS annotate -x FORMAT/AD -Oz -o HG002.joint.QC.vcf.gz
$BCFTOOLS view -s HG003 -c 1 "$COHORT_VCF" \
  | $BCFTOOLS annotate -x FORMAT/AD -Oz -o HG003.joint.QC.vcf.gz
$BCFTOOLS view -s HG004 -c 1 "$COHORT_VCF" \
  | $BCFTOOLS annotate -x FORMAT/AD -Oz -o HG004.joint.QC.vcf.gz

$BCFTOOLS index -f -t HG002.joint.QC.vcf.gz
$BCFTOOLS index -f -t HG003.joint.QC.vcf.gz
$BCFTOOLS index -f -t HG004.joint.QC.vcf.gz
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
      Query source: <code>cohort.gv.QC.vcf.gz</code>, extracted into HG002, HG003, and HG004<br />
      Truth: GIAB v4.2.1 GRCh38 per-sample truth VCFs<br />
      Confident regions: GIAB per-sample <code>*_benchmark_noinconsistent.bed</code> files
    </p>
  </div>
  <div className="giabResultsStatus">F1 = harmonic mean; R = recall; P = precision</div>
</div>

<div className="giabHeatmap">
  <div className="giabHeatmapHeader giabHeatmapCorner">Sample</div>
  <div className="giabHeatmapHeader">SNP ALL</div>
  <div className="giabHeatmapHeader">SNP PASS</div>
  <div className="giabHeatmapHeader">INDEL ALL</div>
  <div className="giabHeatmapHeader">INDEL PASS</div>
  <div className="giabHeatmapSample giabHeatmapHg002">HG002 <span>joint WGS</span></div>
  <div className="giabHeatmapHg002 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 99.15%</span>
    <span><b>R</b> 99.27%</span>
    <span><b>P</b> 99.02%</span>
  </div>
  <div className="giabHeatmapHg002 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 94.93%</span>
    <span><b>R</b> 90.41%</span>
    <span><b>P</b> 99.91%</span>
  </div>
  <div className="giabHeatmapHg002 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 98.57%</span>
    <span><b>R</b> 98.46%</span>
    <span><b>P</b> 98.68%</span>
  </div>
  <div className="giabHeatmapHg002 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 95.71%</span>
    <span><b>R</b> 92.49%</span>
    <span><b>P</b> 99.15%</span>
  </div>
  <div className="giabHeatmapSample giabHeatmapHg003">HG003 <span>joint WGS</span></div>
  <div className="giabHeatmapHg003 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 99.13%</span>
    <span><b>R</b> 99.28%</span>
    <span><b>P</b> 98.98%</span>
  </div>
  <div className="giabHeatmapHg003 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 94.93%</span>
    <span><b>R</b> 90.44%</span>
    <span><b>P</b> 99.90%</span>
  </div>
  <div className="giabHeatmapHg003 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 98.60%</span>
    <span><b>R</b> 98.55%</span>
    <span><b>P</b> 98.64%</span>
  </div>
  <div className="giabHeatmapHg003 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 95.47%</span>
    <span><b>R</b> 92.07%</span>
    <span><b>P</b> 99.13%</span>
  </div>
  <div className="giabHeatmapSample giabHeatmapHg004">HG004 <span>joint WGS</span></div>
  <div className="giabHeatmapHg004 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 99.12%</span>
    <span><b>R</b> 99.27%</span>
    <span><b>P</b> 98.98%</span>
  </div>
  <div className="giabHeatmapHg004 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 94.91%</span>
    <span><b>R</b> 90.40%</span>
    <span><b>P</b> 99.90%</span>
  </div>
  <div className="giabHeatmapHg004 giabHeatmapCell giabHeatmapExcellent">
    <span><b>F1</b> 98.52%</span>
    <span><b>R</b> 98.45%</span>
    <span><b>P</b> 98.59%</span>
  </div>
  <div className="giabHeatmapHg004 giabHeatmapCell giabHeatmapGood">
    <span><b>F1</b> 95.39%</span>
    <span><b>R</b> 91.98%</span>
    <span><b>P</b> 99.07%</span>
  </div>
</div>

Values come from the per-sample `*.joint_wgs.happy.summary.csv` files. PASS
columns represent the final hard-filtered call set; ALL columns include all
compared records.

<details className="cbicallCodeDetails">
<summary>Full hap.py summary statistics</summary>
<div>

HG002:

| Type | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | Recall | Precision | Frac_NA | F1 | Truth Ti/Tv | Query Ti/Tv | Truth het/hom | Query het/hom |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| INDEL | ALL | 525,469 | 517,372 | 8,097 | 922,966 | 7,235 | 376,332 | 3,655 | 864 | 0.984591 | 0.986764 | 0.407742 | 0.985676 |  |  | 1.528276 | 1.941014 |
| INDEL | PASS | 525,469 | 486,027 | 39,442 | 823,560 | 4,353 | 311,706 | 3,391 | 474 | 0.924939 | 0.991496 | 0.378486 | 0.957062 |  |  | 1.528276 | 1.734204 |
| SNP | ALL | 3,365,127 | 3,340,463 | 24,664 | 3,989,474 | 32,921 | 615,168 | 2,859 | 3,133 | 0.992671 | 0.990244 | 0.154198 | 0.991456 | 2.100128 | 1.946593 | 1.581196 | 1.745722 |
| SNP | PASS | 3,365,127 | 3,042,498 | 322,629 | 3,253,076 | 2,599 | 207,054 | 351 | 359 | 0.904126 | 0.999147 | 0.063649 | 0.949264 | 2.100128 | 2.056778 | 1.581196 | 1.655178 |

HG003:

| Type | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | Recall | Precision | Frac_NA | F1 | Truth Ti/Tv | Query Ti/Tv | Truth het/hom | Query het/hom |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| INDEL | ALL | 504,501 | 497,196 | 7,305 | 922,071 | 7,125 | 396,750 | 3,573 | 955 | 0.985520 | 0.986437 | 0.430281 | 0.985978 |  |  | 1.489759 | 1.902552 |
| INDEL | PASS | 504,501 | 464,503 | 39,998 | 815,696 | 4,262 | 326,473 | 3,259 | 508 | 0.920718 | 0.991288 | 0.400239 | 0.954701 |  |  | 1.489759 | 1.679533 |
| SNP | ALL | 3,327,496 | 3,303,506 | 23,990 | 3,972,071 | 34,168 | 633,545 | 2,762 | 3,811 | 0.992790 | 0.989766 | 0.159500 | 0.991276 | 2.102576 | 1.944849 | 1.535137 | 1.705334 |
| SNP | PASS | 3,327,496 | 3,009,395 | 318,101 | 3,238,413 | 3,059 | 225,115 | 309 | 423 | 0.904402 | 0.998985 | 0.069514 | 0.949344 | 2.102576 | 2.054748 | 1.535137 | 1.603981 |

HG004:

| Type | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | Recall | Precision | Frac_NA | F1 | Truth Ti/Tv | Query Ti/Tv | Truth het/hom | Query het/hom |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| INDEL | ALL | 510,519 | 502,608 | 7,911 | 926,086 | 7,474 | 394,450 | 3,967 | 913 | 0.984504 | 0.985942 | 0.425932 | 0.985222 |  |  | 1.516131 | 1.951487 |
| INDEL | PASS | 510,519 | 469,556 | 40,963 | 819,085 | 4,629 | 323,924 | 3,643 | 499 | 0.919762 | 0.990652 | 0.395471 | 0.953892 |  |  | 1.516131 | 1.722373 |
| SNP | ALL | 3,346,610 | 3,322,231 | 24,379 | 4,005,943 | 34,338 | 648,579 | 2,597 | 3,958 | 0.992715 | 0.989772 | 0.161904 | 0.991242 | 2.100775 | 1.943746 | 1.586872 | 1.764840 |
| SNP | PASS | 3,346,610 | 3,025,310 | 321,300 | 3,259,928 | 3,131 | 230,687 | 259 | 445 | 0.903992 | 0.998966 | 0.070764 | 0.949109 | 2.100775 | 2.053567 | 1.586872 | 1.656305 |

</div>
</details>

## Trio Mendelian Consistency

Mendelian consistency is an orthogonal trio-level QC on the final joint-genotyped
VCF. Unlike hap.py, it does not compare each sample with an external truth VCF;
it checks whether the joint HG002/HG003/HG004 genotypes are compatible with the
expected son/father/mother relationship.

The check is performed on autosomes only, with HG002 as child, HG003 as father,
and HG004 as mother:

<details className="cbicallCodeDetails">
<summary>Run bcftools mendelian2</summary>
<div>

```bash
BCFTOOLS=/path/to/bcftools
export BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
COHORT_VCF=/path/to/cohort.gv.QC.vcf.gz
AUTOSOMES=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22

$BCFTOOLS +mendelian2 "$COHORT_VCF" \
  -r "$AUTOSOMES" \
  -p 1X:HG002,HG003,HG004 \
  -m c

$BCFTOOLS +mendelian2 "$COHORT_VCF" \
  -r "$AUTOSOMES" \
  -i 'FILTER="PASS"' \
  -p 1X:HG002,HG003,HG004 \
  -m c
```

</div>
</details>

| Call set | Good trio sites | Mendelian-error sites | Missing trio sites | Violation rate |
| --- | ---: | ---: | ---: | ---: |
| ALL | 6,251,905 | 141,293 | 299,751 | 2.210% |
| PASS | 5,258,609 | 64,094 | 92,655 | 1.204% |

`Good trio sites` are variant sites where the child genotype can be explained
by one allele inherited from each parent. `Mendelian-error sites` are sites where
the three genotypes are incompatible with the recorded family structure, for
example both parents are `0/0` while the child is `0/1`. `Missing trio sites`
have at least one missing genotype and are not used for the rate calculation.
The violation rate is calculated as `nmerr / (ngood + nmerr)`, excluding missing
trio sites from the denominator.


## Exact Benchmark Commands

The protocol above uses portable paths. The example local paths below show one
completed benchmark layout.

<details className="cbicallCodeDetails">
<summary>Example local commands used for the benchmark</summary>
<div>

Local roots:

```bash
GIAB_RUN_ROOT=/media/mrueda/4TBB/GIAB
GIAB_FASTQ_ROOT=$GIAB_RUN_ROOT/AshkenazimTrio_fastq
GIAB_WORK_ROOT=/media/mrueda/4TBA/GIAB
CBICALL_DATA_ROOT=/media/mrueda/4TBB/cbicall-data
HG38_REF_DIR=$CBICALL_DATA_ROOT/Databases/GATK_bundle/hg38
HAPPY_DIR=$GIAB_WORK_ROOT/happy_results
TRUTH_DIR=$GIAB_WORK_ROOT/truthsets
BCFTOOLS=/media/mrueda/2TBS/NGSutils/bcftools-1.21-103_x86_64/bcftools
export BCFTOOLS_PLUGINS=/media/mrueda/2TBS/NGSutils/bcftools-1.21-103_x86_64/plugins
AUTOSOMES=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22
```

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

```bash
cat > sample_map.tsv <<EOF
HG002	$GIAB_FASTQ_ROOT/HG002_NA24385_son/cbicall_bash_gatk-4.6_wgs_single_hg38_178120808914892/02_varcall/HG002.hc.g.vcf.gz
HG003	$GIAB_FASTQ_ROOT/HG003_NA24149_father/cbicall_bash_gatk-4.6_wgs_single_hg38_178128498815686/02_varcall/HG003.hc.g.vcf.gz
HG004	$GIAB_FASTQ_ROOT/HG004_NA24143_mother/cbicall_bash_gatk-4.6_wgs_single_hg38_178128547716032/02_varcall/HG004.hc.g.vcf.gz
EOF
```

Local hap.py image:

```bash
cd "$GIAB_WORK_ROOT/hap.py"
docker build -t local/hap.py .
```

RTG SDF:

```bash
docker run --rm \
  -v "$HG38_REF_DIR":/ref \
  local/hap.py \
  /opt/hap.py/libexec/rtg-tools-install/rtg format \
  -o /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
  /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
```

Extract the joint-genotyped query VCFs:

```bash
COHORT_VCF=$GIAB_RUN_ROOT/cbicall_bash_gatk-4.6_wgs_cohort_hg38_178167960428471/02_varcall/cohort.gv.QC.vcf.gz

$BCFTOOLS view -s HG002 -c 1 "$COHORT_VCF" \
  | $BCFTOOLS annotate -x FORMAT/AD -Oz -o HG002.joint.QC.vcf.gz
$BCFTOOLS view -s HG003 -c 1 "$COHORT_VCF" \
  | $BCFTOOLS annotate -x FORMAT/AD -Oz -o HG003.joint.QC.vcf.gz
$BCFTOOLS view -s HG004 -c 1 "$COHORT_VCF" \
  | $BCFTOOLS annotate -x FORMAT/AD -Oz -o HG004.joint.QC.vcf.gz

$BCFTOOLS index -f -t HG002.joint.QC.vcf.gz
$BCFTOOLS index -f -t HG003.joint.QC.vcf.gz
$BCFTOOLS index -f -t HG004.joint.QC.vcf.gz
```

Extracted joint-WGS hap.py comparisons:

```bash
for SAMPLE in HG002 HG003 HG004; do
  docker run --rm \
    -v "$HAPPY_DIR":/query \
    -v "$HAPPY_DIR":/out \
    -v "$TRUTH_DIR":/truth \
    -v "$HG38_REF_DIR":/ref \
    local/hap.py \
    /opt/hap.py/bin/hap.py \
    /truth/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
    /query/${SAMPLE}.joint.QC.vcf.gz \
    -f /truth/${SAMPLE}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    -r /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
    --engine=vcfeval \
    --engine-vcfeval-template /ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.sdf \
    --threads 12 \
    -o /out/${SAMPLE}.joint_wgs.happy
done
```

Autosomal Mendelian consistency:

```bash
$BCFTOOLS +mendelian2 "$COHORT_VCF" \
  -r "$AUTOSOMES" \
  -p 1X:HG002,HG003,HG004 \
  -m c

$BCFTOOLS +mendelian2 "$COHORT_VCF" \
  -r "$AUTOSOMES" \
  -i 'FILTER="PASS"' \
  -p 1X:HG002,HG003,HG004 \
  -m c
```

</div>
</details>
