# Overview

These pipelines extract mitochondrial reads from exome data and run **MToolBox** to generate mtDNA variant calls, annotations, and heteroplasmy estimates.

There are two processing modes:

- **Single-sample analysis** (`mit_single`)
- **Cohort / family analysis** (`mit_cohort`)

Both consume **WES single-sample outputs** and assume [this nomenclature](../help/naming-conventions.md).

---

## Choosing a Pipeline

| Use Case | Pipeline | Description |
|----------|----------|-------------|
| Analyze **one individual** | `mit_single` | Fast, sample-specific mtDNA variant calling + HF/DP/GT extraction |
| Analyze **a full family or cohort** | `mit_cohort` | Joint variant calling across samples; useful for transmission and segregation checks |

---

# mtDNA Pipelines

=== "Single-Sample (`mit_single`)"

    ## mtDNA Single-Sample Pipeline

    ??? Example "See Bash pipeline"
        ```bash
        --8<-- "https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/workflows/bash/gatk-3.5/mit_single.sh"
        ```

    ### Workflow Diagram

    ??? Example "View diagram"
        ```mermaid
        flowchart TD
            A["Input WES BAM from single sample pipeline"]
            B["Extract mitochondrial reads (chrM/MT)"]
            C["Create mtDNA-only BAM"]
            D["Run MToolBox (RSRS reference)"]
            E["Get VCF + prioritized variants"]
            F["Extract GT/DP/HF for each variant"]
            G["Write mit_prioritized_variants.txt"]
            A --> B --> C --> D --> E --> F --> G
        ```

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
    - `parameters.sh` provides:
      - `REF`
      - `SAM` path
      - `MTOOLBOXDIR`

    ### Outputs

    | File | Description |
    |------|-------------|
    | `VCF_file.vcf` | mtDNA VCF from MToolBox |
    | `prioritized_variants.txt` | Raw prioritized list |
    | `mit_prioritized_variants.txt` | Final prioritized list with GT/DP/HF |

=== "Cohort (`mit_cohort`)"

    ## mtDNA Cohort Pipeline

    ??? Example "See Bash pipeline"
        ```bash
        --8<-- "https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/workflows/bash/gatk-3.5/mit_cohort.sh"
        ```

    ### Workflow Diagram

    ??? Example "View diagram"
        ```mermaid
        flowchart TD
            A["Input all WES single-sample directories"]
            B["Extract mtDNA BAMs for each sample"]
            C["mtDNA BAM collection"]
            D["Run MToolBox jointly on cohort"]
            E["Get cohort VCF + prioritized variants"]
            F["Extract GT/DP/HF for all samples"]
            G["Write mit_prioritized_variants.txt (cohort)"]
            A --> B --> C --> D --> E --> F --> G
        ```

    ### Summary

    `mit_cohort` processes **all samples together**, generating a **joint mtDNA variant table** useful for:

    - Family mtDNA transmission studies  
    - Maternal-lineage analysis  
    - Variant segregation and consistency checks  

    ### Inputs

    - Run inside: `cbicall_bash_mit_cohort_*`
    - Expects WES single-sample directories:
      `../../??????????_ex/cbicall_bash_wes_single*/01_bam/`
    - `parameters.sh` defines `REF`, `SAM`, `MTOOLBOXDIR`

    ### Outputs

    | File | Description |
    |------|-------------|
    | `VCF_file.vcf` | mtDNA VCF for full cohort |
    | `prioritized_variants.txt` | MToolBox-prioritized list |
    | `mit_prioritized_variants.txt` | Joint variant table with GT/DP/HF |

# When to Use Each Pipeline

### Use `mit_single` when:

- You need results for **one individual**.
- You are adding or reprocessing a single relative.
- You want faster turnaround for a standalone case.

### Use `mit_cohort` when:

- You want a **joint variant table** across multiple individuals.
- You are analyzing **mtDNA inheritance** within a family.
- You need cohort-level prioritization or downstream analysis.

---

# Background Information

**SG-ADVISER mtDNA** builds on [MToolbox v1.0](https://github.com/mitoNGS/MToolBox) and performs:

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


