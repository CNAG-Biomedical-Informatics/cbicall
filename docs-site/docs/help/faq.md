# Frequently Asked Questions

## WES / WGS

<details>
<summary>What reference genomes are used?</summary>

CBIcall supports:

- **GRCh37 / b37**: GATK-compatible reference genome
- **GRCh38 / hg38**: GATK-compatible reference genome

</details>

<details>
<summary>What capture kits are used for WES?</summary>

- **GATK 3.5**: exome capture is based on Agilent SureSelect.
- **GATK 4.6**: exome and WGS references are based on the GATK bundle.

</details>

## mtDNA (MToolBox)

<details>
<summary>What reference genome is used?</summary>

mtDNA workflows use **RSRS** (`rsrs`), the Reconstructed Sapiens Reference Sequence.

</details>

<details>
<summary>VCF vs. prioritized variants allele notation</summary>

:::note

In rare cases, the allele reported in `prioritized_variants.txt` may differ from the ALT allele reported in the VCF.

The `Variant_Allele` column is generated during annotation and prioritization and does not always follow VCF semantics, where ALT is defined relative to the mapping reference, such as RSRS.

:::

</details>

<details>
<summary>What does <code>GT=1</code> mean in results?</summary>

In variant reports, the **Genotype (`GT`)** field shows the observed allele using **VCF allele indices**:

- `0` = reference allele
- `1` = first alternate allele
- `2`, `3`, ... = additional alternate alleles in multiallelic records

For CBIcall mtDNA reports, **`GT=1` means an ALT allele was detected in that sample**.

Biological interpretation should use:

- **`HF`**: fraction of reads supporting the ALT allele
- **`DP`**: total read depth at the variant position

:::tip

For mtDNA, `GT` tells you **which allele** was detected, not **how much** of it was detected. Use `HF` and `DP` to interpret heteroplasmy or homoplasmy.

:::

</details>

## General

<details>
<summary>What is the difference between native and external workflows?</summary>

**Native CBIcall workflows** are maintained in this repository and launched through a supported workflow backend: Bash, Snakemake, Nextflow, or Cromwell. They use the CBIcall project layout, provenance files, run reports, and usually the CBIcall resource bundle.

**External workflows** are third-party workflows registered in CBIcall. Today this means selected nf-core workflows launched through Nextflow. CBIcall validates the YAML contract, pins the registered workflow, and writes provenance and run reports, while nf-core keeps its own output layout, profiles, containers, and reference-resource assumptions.

See [Native Backends](../backends/native), [nf-core Provider](../backends/nf-core), and [Resource Validation](../usage/resource-validation).

</details>

<details>
<summary>How do I set up <code>cbicall</code> on an HPC system?</summary>

On most HPC systems, Docker is not available. CBIcall is designed to run with **Apptainer**, formerly Singularity, which is the recommended approach for HPC environments.

Apptainer can execute Docker images directly, requires no root privileges, and integrates cleanly with batch schedulers.

In this setup:

- the container image is read-only
- configuration files and workflows are stored in a writable host directory
- native CBIcall resources, when needed, are downloaded outside the container and bind-mounted at runtime

Recommended workflow:

1. Pull the CBIcall container image using Apptainer.
2. Create a writable copy of the CBIcall workflow directory.
3. Choose `workflow_provider: nf-core` for a no-bundle first run, or download the CBIcall resource bundle for native WES/WGS/mtDNA workflows.
4. Run the pipeline with the writable copy, adding the data bind only for native workflows.

See [HPC with Apptainer / Singularity](../installation/apptainer).

</details>

<details>
<summary>Do you have a Slurm + Apptainer example?</summary>

Yes. See the example script:

[run_cbicall_apptainer_slurm.sh](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/examples/scripts/run_cbicall_apptainer_slurm.sh)

</details>

<details>
<summary>How do I cite CBIcall?</summary>

Please cite:

CBIcall: a configuration-driven framework for variant calling in large sequencing cohorts. [Preprint DOI](https://doi.org/10.64898/2026.03.23.713646).

</details>
