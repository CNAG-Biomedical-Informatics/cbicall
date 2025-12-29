# Frequently Asked Questions

## WES / WGS 

??? question "What are the reference genomes used?"

    GRCh37 (b37) - GATK-compatible reference genome

    GRCh38 (hg38) - GATK-compatible reference genome

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)


??? question "What are the capture kits for WES?"

    * For GATK version 3.5: Exome capture is based on Agilent SureSelect.

    * For GATK version 4.6: Exome and WGS reference is based on the GATK bundle (b37).

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)


## mtDNA (MToolBox)

??? question "What is the reference genome used?"

    RSRS (rsrs) - Reconstructed Sapiens Reference Sequence

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)


??? question "What does GT=1 mean in results?"

    In variant reports, the **Genotype (GT)** field shows the observed
    allele using **VCF allele indices**:
    
    -   `0` = reference allele
    -   `1` = first alternate (ALT) allele
    -   `2`, `3`, ... = additional ALT alleles (multiallelic)
    
    For **chrM/MT (mtDNA)**, callers typically encode genotypes as
    **haploid** (not allele pairs).
    
    ## Meaning
    
    -   **`GT = 1` → ALT allele detected in that sample**
    -   No `/` or `|` separator because only one allele index is stored
    -   Biological interpretation relies on:
        -   **`HF`** → heteroplasmy fraction (molecules supporting ALT)
        -   **`DP`** → read depth (total support)
    
    ## Examples
    
    | GT | Interpretation (mtDNA) |
    |---|---|
    | `0` | Only reference allele observed |
    | `1` | ALT allele present (homoplasmic or heteroplasmic, check `HF` + `DP`) |
    | `0/1`, `1/2` *(rare)* | Multiallelic call, still haploid encoding — not diploid zygosity |

    > **TL;DR:** `GT = 1` = ALT detected. Check `HF` and `DP` for biology.
    
    !!! info "Tip" 

        For mtDNA, **`GT` tells you *which allele*, not *how much***.

        Use **`HF` + `DP`** to interpret heteroplasmy or homoplasmy.

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)

## General 

??? question "How do I set up `cbicall` to work on an HPC system?"

    To adapt `cbicall` for your HPC environment, update the file  
    `workflows/bash/gatk_3.5/parameters.sh` so that it reflects your local module
    system, paths, and resource settings.

    Below is the configuration used at **CNAG-HPC**, which you can use as a template:

    ```bash
    --8<-- "https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/workflows/bash/gatk-3.5/parameters.sh"
    ```

??? question "Do you have an example in how to run `cbicall` in **Slurm** HPC?"

    ```bash
    --8<-- "https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/examples/scripts/run_cbicall_slurm.sh"
    ```

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)

??? question "How do I cite **CBICall**?"

    You can cite the **CBICall** paper. Thx!

    !!! Note "Citation"

        CBICall: a pipeline-agnostic framework for variant calling in large DNA-seq cohorts. _Manuscript In preparation._

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)
