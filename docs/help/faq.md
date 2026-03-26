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

??? question "VCF vs. prioritized variants allele notation"

    In rare cases, the allele reported in `prioritized_variants.txt` may differ from the ALT allele reported in the VCF. The `Variant_Allele` column is generated during annotation and prioritization and does not always follow VCF semantics, where ALT is defined relative to the mapping reference (e.g. RSRS).

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)
??? question "What does GT=1 mean in results?"

    In variant reports, the **Genotype (GT)** field shows the observed
    allele using **VCF allele indices**:
    
    -   `0` = reference allele
    -   `1` = first alternate (ALT) allele
    -   `2`, `3`, ... = additional ALT alleles (multiallelic)
    
    ## Meaning
    
    -   **`GT = 1` → ALT allele detected in that sample**
    -   Biological interpretation relies on:
        -   **`HF`** → heteroplasmy fraction (molecules supporting ALT)
        -   **`DP`** → read depth (total support)
    
    !!! info "Tip" 

        For mtDNA, **`GT` tells you *which allele*, not *how much***.

        Use **`HF` + `DP`** to interpret heteroplasmy or homoplasmy.

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)

## General 

??? question "How do I set up `cbicall` to work on an HPC system?"

    On most HPC systems, Docker is not available. Instead, `cbicall` is designed to run
    using **Apptainer (formerly Singularity)**, which is the recommended approach.

    Apptainer can execute Docker images directly and is well suited for HPC environments
    because it requires no root privileges and integrates cleanly with batch schedulers.

    In this setup:

    - the container image is **read-only**
    - all configuration files and workflows are stored in a **writable host directory**
    - external databases are downloaded outside the container and bind-mounted at runtime

    The recommended workflow is:

    1. Pull the CBIcall container image using Apptainer  
    2. Download the required databases on the host filesystem  
    3. Create a writable copy of the CBIcall workflow directory  
    4. Run the pipeline by bind-mounting the writable copy and data directory  

    This approach avoids manual dependency installation, improves reproducibility,
    and works on both interactive and batch-based HPC systems.

    Step-by-step instructions are provided in:

    **⬇️ Installation → HPC (Apptainer / Singularity)**

??? question "Do you have an example in how to run `cbicall` in **Slurm** HPC with **apptainer**?"

    ```bash
    --8<-- "https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/examples/scripts/run_cbicall_apptainer_slurm.sh"
    ```

    ##### last change 2026-01-14 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)

??? question "How do I cite **CBIcall**?"

    You can cite the **CBIcall** paper. Thx!

    !!! Note "Citation"
        CBIcall: a configuration-driven framework for variant calling in large sequencing cohorts. [Preprint DOI](https://doi.org/10.64898/2026.03.23.713646).

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)
