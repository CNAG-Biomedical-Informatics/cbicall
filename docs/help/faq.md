# Frequently Asked Questions

## WES / WGS 

??? question "What is the reference genome used?"

    GRCh37 (b37/hs37d5) - GATK-compatible reference genome

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)


??? question "What are the capture kits for WES?"

    * For GATK version 3.5: Exome capture is based on Agilent SureSelect.

    * For GATK version 4.6: Exome and WGS reference is based on the GATK bundle (b37).

    ##### last change 2025-10-15 by Manuel Rueda [:fontawesome-brands-github:](https://github.com/mrueda)


## mtDNA 

## General 

??? question "How do I set up `cbicall` to work on an HPC system?"

    To adapt `cbicall` for your HPC environment, update the file  
    `workflows/bash/gatk_3.5/parameters.sh` so that it reflects your local module
    system, paths, and resource settings.

    Below is the configuration used at **CNAG-HPC**, which you can use as a template:

    ```bash
    --8<-- "https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/workflows/bash/gatk-3.5/cnag-hpc-parameters.sh"
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
