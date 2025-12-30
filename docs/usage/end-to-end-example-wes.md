# End-to-end examples (GATK 4.6)

!!! note "Prerequisites"
    Installation, reference bundles, and all dependencies must be completed beforehand.  

    [➡️ Installation](../installation/non-containerized.md){ .md-button .md-button--primary }

---

=== "WES single-sample run"

    This example demonstrates how to run CBIcall on a real WES sample from FASTQ files through final VCF and QC outputs.

    ## 1. Prepare your FASTQ files
    
    
    CBIcall expects paired-end FASTQ files with a shared prefix, for example:
    
    ```
    # Project    / Sample (Proband WES)
    CNAG999_exome/CNAG99901P_ex/
      CNAG99901P_ex_S1_L001_R1_001.fastq.gz
      CNAG99901P_ex_S1_L001_R2_001.fastq.gz
    ```
    
    !!! Note "Note on nomenclature"
    
        Please see [this page](../help/naming-conventions.md).
    ---
    
    ## 2. Create a parameters file
    
    Create a YAML file, e.g. `wes_single.yaml`:
    
    ```yaml
    mode:            single
    pipeline:        wes
    workflow_engine: bash
    gatk_version:    gatk-4.6
    sample:          CNAG999_exome/CNAG99901P_ex
    genome:          b37
    cleanup_bam:     false
    ```
    
    Notes:
    
    - `mode` selects single-sample or cohort (joint genotyping).  
    - `pipeline` switches between WES, WGS or mtDNA.  
    - `workflow_engine` chooses the backend (bash or snakemake).  
    
    ---
    
    ## 3. Run CBIcall
    
    ```bash
    bin/cbicall -p wes_single.yaml -t 4
    ```
    
    - `-p` selects the YAML parameters file  
    - `-t` sets the number of threads
    
    ---
    
    ## 4. Inspect outputs
    
    After completion, you will find:
    
    ```
    CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_gatk-4.6_*/
      01_bam/
      02_varcall/
      03_stats/
      logs/
    ```
    
    Where:
    
    - VCF files are stored in `02_varcall/`  
    - QC metrics (coverage, sample stats, sex prediction) are in `03_stats`  
    - Logs for all pipeline steps are under `logs/`  
    
    These files are ready for downstream analysis, annotation or integration with cohort-level studies.
    
    ---
    
    For advanced parameters, multi-sample analyses, mtDNA workflows and troubleshooting, see the **Usage** and **FAQ** sections.

=== "WES cohort run"

    
    !!! Warning "Important"

        In order to run a `cohort` based calculation you first have to create `GVCF` for each sample. This is being done by running `wes` mode `single`.

 
    ## 1. Create a sample map file like the one we display below:
    
    
    ```bash
    --8<-- "https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/examples/input/sample_map.tsv"
    ```
 
    **GATK** needs **absolute paths** for the files.
    
    ---
    
    ## 2. Create a parameters file
    
    Create a YAML file, e.g. `wes_cohort.yaml`:
    
    ```yaml
    mode:            cohort
    pipeline:        wes
    workflow_engine: bash
    gatk_version:    gatk-4.6
    genome:          b37
    sample_map:      ./sample_map.tsv
    ```
    
    ---
    
    ## 3. Run CBIcall
    
    ```bash
    bin/cbicall -p wes_cohort.yaml -t 4
    ```
    
    - `-p` selects the YAML parameters file  
    - `-t` sets the number of threads
    
    ---
    
    ## 4. Inspect outputs
    
    After completion, you will find:
    
    ```
    cbicall_bash_wes_cohort_gatk-4.6_*/
      02_varcall/
      logs/
    ```
    
    Where:
    
    - Final VCF files are stored in `02_varcall/`  
    - Logs for all pipeline steps are under `logs/`  
