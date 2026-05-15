# Outputs

CBIcall writes each run into a run directory named from the workflow choices:

```text
cbicall_<engine>_<pipeline>_<mode>_<genome>_<gatk-version>_<run-id>/
```

The exact files depend on the selected pipeline and mode. The tables below are derived from the workflow output definitions and the checked-in example runs.

:::tip[Most users only need these]
- WES/WGS: use the final QC VCF in `02_varcall/`.
- WES/WGS single-sample runs: keep the gVCF if you plan cohort joint genotyping.
- mtDNA: use `01_mtoolbox/mit_prioritized_variants.txt` and the browser report in `02_browser/`.
:::

## Common Run Files

| File | Meaning |
| --- | --- |
| `log.json` | Structured record of CLI arguments, resolved configuration, selected `profile`, compact `resource_bundle` provenance, and runtime parameters. |
| `run-report.json` | Compact execution report with status, elapsed time, workflow, profile, bundle fingerprint, and workflow log path. |
| `<engine>_<pipeline>_<mode>_<genome>_<gatk-version>.log` | Main wrapper log for Bash runs. |
| `logs/*.log` | Per-rule or per-step logs for Snakemake/GATK 4.6 workflows. |

Use `config.resource_bundle.fingerprint` inside `log.json` to check whether two runs used the same declared external dependency set.

## WES/WGS Single-Sample

Applies to `pipeline: wes` or `pipeline: wgs` with `mode: single`.

### Recommended Files

| File | Use |
| --- | --- |
| `02_varcall/<id>.hc.QC.vcf.gz` | Final filtered single-sample VCF. This is the primary variant file for downstream interpretation. |
| `02_varcall/<id>.hc.QC.vcf.gz.tbi` | Tabix index for the final VCF. |
| `02_varcall/<id>.hc.g.vcf.gz` | Per-sample gVCF. Use this as input for cohort joint genotyping. |
| `02_varcall/<id>.hc.g.vcf.gz.tbi` | Tabix index for the gVCF. |
| `03_stats/<id>.coverage.txt` | Coverage summary. |
| `03_stats/<id>.sex.txt` | Sex inference result from the final VCF. |

<details>
<summary>Intermediate files</summary>

| File | Meaning |
| --- | --- |
| `01_bam/<fastq-prefix>.rg.bam` | Lane-level BAM after alignment and read-group assignment. |
| `01_bam/<id>.rg.merged.bam` | BAM after merging lanes for the sample. |
| `01_bam/<id>.rg.merged.dedup.bam` | Duplicate-marked BAM. |
| `01_bam/<id>.rg.merged.dedup.metrics.txt` | Duplicate-marking metrics. |
| `01_bam/<id>.rg.merged.dedup.recal.table` | BQSR recalibration table. |
| `01_bam/<id>.rg.merged.dedup.recal.bam` | Recalibrated BAM used for variant calling. |
| `01_bam/*.bai` or `01_bam/*.bam.bai` | BAM indexes. |
| `02_varcall/<id>.hc.raw.vcf.gz` | Raw VCF from `GenotypeGVCFs`. |
| `02_varcall/<id>.hc.raw.vcf.gz.tbi` | Tabix index for the raw VCF. |

</details>

<details>
<summary>Conditional VQSR files</summary>

These appear only when the run has enough SNPs or indels to build VQSR models.

| File | Meaning |
| --- | --- |
| `02_varcall/<id>.hc.snp.recal.vcf.gz` | SNP VQSR model output. |
| `02_varcall/<id>.hc.snp.tranches.txt` | SNP VQSR tranche diagnostics. |
| `02_varcall/<id>.hc.post_snp.vcf.gz` | VCF after applying SNP VQSR. |
| `02_varcall/<id>.hc.indel.recal.vcf.gz` | INDEL VQSR model output. |
| `02_varcall/<id>.hc.indel.tranches.txt` | INDEL VQSR tranche diagnostics. |
| `02_varcall/<id>.hc.vqsr.vcf.gz` | VCF after applying SNP and INDEL VQSR. |

</details>

:::note
If VQSR is skipped because there are too few variants, the final `*.hc.QC.vcf.gz` is still produced by hard filtering.
:::

## WES/WGS Cohort

Applies to `pipeline: wes` or `pipeline: wgs` with `mode: cohort`.

### Recommended Files

| File | Use |
| --- | --- |
| `02_varcall/cohort.gv.QC.vcf.gz` | Final filtered cohort VCF. This is the primary joint-genotyped variant file. |
| `02_varcall/cohort.gv.QC.vcf.gz.tbi` | Tabix index for the final cohort VCF. |

<details>
<summary>Intermediate files</summary>

| File | Meaning |
| --- | --- |
| `02_varcall/cohort.genomicsdb.<run-id>/` | GenomicsDB workspace used by `GenomicsDBImport`. |
| `02_varcall/genomicsdbimport.done` | Snakemake marker showing that GenomicsDB import completed. |
| `02_varcall/cohort.gv.raw.vcf.gz` | Raw cohort VCF from `GenotypeGVCFs`. |
| `02_varcall/cohort.gv.raw.vcf.gz.tbi` | Tabix index for the raw cohort VCF. |
| `logs/01_genomicsdbimport.log` | GenomicsDB import log. |
| `logs/02_genotype_gvcfs.log` | Cohort genotyping log. |
| `logs/03_vqsr_and_qc.log` | VQSR and final filtering log. |

</details>

<details>
<summary>Conditional VQSR files</summary>

| File | Meaning |
| --- | --- |
| `02_varcall/cohort.snp.recal.vcf.gz` | SNP VQSR model output. |
| `02_varcall/cohort.snp.tranches.txt` | SNP VQSR tranche diagnostics. |
| `02_varcall/cohort.post_snp.vcf.gz` | VCF after applying SNP VQSR. |
| `02_varcall/cohort.indel.recal.vcf.gz` | INDEL VQSR model output. |
| `02_varcall/cohort.indel.tranches.txt` | INDEL VQSR tranche diagnostics. |
| `02_varcall/cohort.vqsr.vcf.gz` | VCF after applying SNP and INDEL VQSR. |

</details>

## mtDNA Single-Sample

Applies to `pipeline: mit` with `mode: single`.

### Recommended Files

| File | Use |
| --- | --- |
| `01_mtoolbox/mit_prioritized_variants.txt` | Final prioritized mtDNA variant report with `GT`, `DP`, and heteroplasmic fraction columns appended by CBIcall. |
| `01_mtoolbox/VCF_file.vcf` | mtDNA VCF from MToolBox. |
| `02_browser/<run-id>.html` | Interactive HTML report. |
| `02_browser/mit.json` | JSON used by the browser report. |
| `02_browser/README.txt` | Local instructions for opening the browser report. |

<details>
<summary>Intermediate files</summary>

| File | Meaning |
| --- | --- |
| `01_mtoolbox/<id>-DNA_MIT.bam` | Extracted mitochondrial BAM used as MToolBox input. |
| `01_mtoolbox/<id>-DNA_MIT.bam.bai` | BAM index. |
| `01_mtoolbox/prioritized_variants.txt` | Raw MToolBox prioritized variant list before CBIcall appends genotype/depth/HF fields. |
| `01_mtoolbox/mit.raw.json` | Raw JSON conversion of the final prioritized report. |
| `01_mtoolbox/mt_classification_best_results.csv` | MToolBox haplogroup/classification output. |
| `01_mtoolbox/processed_fastq.tar.gz` | MToolBox processed FASTQ archive. |
| `01_mtoolbox/summary_*.txt` | MToolBox run summary. |
| `01_mtoolbox/OUT_*/` | MToolBox working directory with alignment, pileup, coverage, and annotation intermediates. |

</details>

## mtDNA Cohort

Applies to `pipeline: mit` with `mode: cohort`.

The cohort workflow uses the same output directories as mtDNA single-sample mode, but extracts mtDNA BAMs from all matching sibling sample directories before running MToolBox jointly.

### Recommended Files

| File | Use |
| --- | --- |
| `01_mtoolbox/mit_prioritized_variants.txt` | Final joint mtDNA variant report with per-sample `GT`, `DP`, and heteroplasmic fraction fields. |
| `01_mtoolbox/VCF_file.vcf` | Joint mtDNA VCF from MToolBox. |
| `02_browser/<run-id>.html` | Interactive cohort HTML report. |
| `02_browser/mit.json` | JSON used by the browser report. |
| `02_browser/README.txt` | Local instructions for opening the browser report. |

<details>
<summary>Intermediate files</summary>

| File | Meaning |
| --- | --- |
| `01_mtoolbox/<sample-id>-DNA_MIT.bam` | Extracted mitochondrial BAM for each cohort sample. |
| `01_mtoolbox/<sample-id>-DNA_MIT.bam.bai` | BAM index for each extracted mitochondrial BAM. |
| `01_mtoolbox/prioritized_variants.txt` | Raw MToolBox prioritized variant list. |
| `01_mtoolbox/missing_variants.txt` | Temporary variant list used while appending cohort genotype/depth/HF fields. |
| `01_mtoolbox/mit.raw.json` | Raw JSON conversion of the final prioritized report. |
| `01_mtoolbox/OUT_*/` | MToolBox working directories and intermediate outputs. |

</details>
