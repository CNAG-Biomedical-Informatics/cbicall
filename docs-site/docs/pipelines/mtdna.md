import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# mtDNA Pipelines

These pipelines extract mitochondrial reads from exome data and run **MToolBox** to generate mtDNA variant calls, annotations, and heteroplasmy estimates.

There are two processing modes:

- **Single-sample analysis**: `mit_single`
- **Cohort / family analysis**: `mit_cohort`

Both consume **WES single-sample outputs** and assume [this nomenclature](../help/naming-conventions).

---

## Choosing a Pipeline

| Use Case | Pipeline | Description |
|----------|----------|-------------|
| Analyze **one individual** | `mit_single` | Fast, sample-specific mtDNA variant calling + HF/DP/GT extraction |
| Analyze **a full family or cohort** | `mit_cohort` | Joint variant calling across samples; useful for transmission and segregation checks |

---

## Workflow Details

<Tabs groupId="workflow-mode">
<TabItem value="single-sample-mit-single" label="Single-Sample: mit_single" default>

## mtDNA Single-Sample Pipeline

**Source:** [View source](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-3.5/mit_single.sh)

### Workflow Diagram

![mtDNA single-sample workflow](/img/diagram-mtdna-single.svg)

### Summary

The `mit_single` pipeline processes **one individual** at a time.  
It extracts mtDNA reads, runs MToolBox, and enriches the prioritized variants with:

- **GT** — genotype  
- **DP** — depth  
- **HF** — heteroplasmic fraction  

### Inputs

- Run inside directory: `cbicall_bash_mit_single_*`
- Needs WES single-sample directory:
  `../../cbicall_bash_wes_single*/01_bam/input.merged.filtered.realigned.fixed.bam`
- `env.sh` provides:
  - `REF`
  - `SAM` path
  - `MTOOLBOXDIR`

### Outputs

| File | Description |
|------|-------------|
| `VCF_file.vcf` | mtDNA VCF from MToolBox |
| `prioritized_variants.txt` | Raw prioritized list |
| `mit_prioritized_variants.txt` | Final prioritized list with GT/DP/HF |

</TabItem>
<TabItem value="cohort-mit-cohort" label="Cohort: mit_cohort">

## mtDNA Cohort Pipeline

**Source:** [View source](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-3.5/mit_cohort.sh)

### Workflow Diagram

![mtDNA cohort workflow](/img/diagram-mtdna-cohort.svg)

### Summary

`mit_cohort` processes **all samples together**, generating a **joint mtDNA variant table** useful for:

- Family mtDNA transmission studies  
- Maternal-lineage analysis  
- Variant segregation and consistency checks  

### Inputs

- Run inside: `cbicall_bash_mit_cohort_*`
- Expects WES single-sample directories:
  `../../??????????_ex/cbicall_bash_wes_single*/01_bam/`
- `env.sh` defines `REF`, `SAM`, `MTOOLBOXDIR`

### Outputs

| File | Description |
|------|-------------|
| `VCF_file.vcf` | mtDNA VCF for full cohort |
| `prioritized_variants.txt` | MToolBox-prioritized list |
| `mit_prioritized_variants.txt` | Joint variant table with GT/DP/HF |

</TabItem>
</Tabs>
## When to Use Each Pipeline

### Use `mit_single` when:

- You need results for **one individual**.
- You are adding or reprocessing a single relative.
- You want faster turnaround for a standalone case.

### Use `mit_cohort` when:

- You want a **joint variant table** across multiple individuals.
- You are analyzing **mtDNA inheritance** within a family.
- You need a cohort-level mtDNA table for downstream review or comparison.

---

## Background Information

**CBIcall mtDNA** builds on [MToolbox v1.0](https://github.com/mitoNGS/MToolBox) and performs:

---
### 1. Preprocessing (PicardTools)
Converts BAM → FASTQ using:
- *SortSam.jar*
- *MarkDuplicates.jar*
- *SamFormatConverter.jar*
([PicardTools](http://picard.sourceforge.net))

---
### 2. Alignment
- Aligns reads to **RSRS** via `mapExome.py`
- Uses **GSNAP** ([2015-12-31.v7](http://research-pub.gene.com/gmap/))

---
### 3. Variant Calling & Annotation (MToolBox)
Pipeline steps include:

- *mpileup* (SAMtools)
- *mtVariantCaller.py*
- *VCFoutput.py* (with [PyVCF](https://github.com/jamescasbon/PyVCF))
- *mt-classifier.py* (haplogroup prediction)
- *variants_functional_annotation.py*
- *prioritization.py*
- *summary.py*

---
### Reference

1. Calabrese C. *et al.*  
   **MToolBox: a highly automated pipeline for heteroplasmy annotation and prioritization analysis of human mitochondrial variants in high-throughput sequencing.**  
   *Bioinformatics* (2014).  
   [Read paper](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu483)
