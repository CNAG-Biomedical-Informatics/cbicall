# CBIcall

![CBIcall](/img/cbicall-logo.png)

<p align="center"><em>Configuration-driven framework for reproducible cohort-scale variant calling</em></p>

**CBIcall** (**C**NAG **B**iomedical **I**nformatics framework for variant **call**ing) is a configuration-driven framework for reproducible variant calling in large sequencing cohorts.

:::tip[In one sentence]
CBIcall validates a parameters YAML file against approved workflow and resource contracts, creates a deterministic run directory, and records structured reports for audit and run comparison.
:::


## Execution Contract

CBIcall separates routine user input from framework-managed implementation and
resource definitions.

| Contract layer | Purpose |
| --- | --- |
| Parameters YAML | User analysis intent: inputs, pipeline, mode, genome, backend, and runtime options. |
| Workflow registry | Approved workflow implementations and backend-specific entrypoints. |
| Resource catalog | External references, tool/resource bundles, compatibility rules, and resource identity. |

This separation is the basis for validation, provenance capture, and repeatable
run comparison across computing environments.

## What CBIcall Does

CBIcall does not re-implement alignment or variant-calling algorithms; those
steps remain in curated pipeline implementations built from established tools such as BWA,
GATK, MToolBox, Snakemake, Nextflow, Cromwell, and selected nf-core pipelines.

- validates the **parameters YAML** and compatibility contract before launch
- resolves **Bash**, **Snakemake**, **Nextflow**, or **Cromwell** workflow backends
- checks the selected workflow against the **resource catalog** when resources are required
- records **logs, run reports, workflow fingerprints, resource identity, output inventories, normalized VCF hashes, and run comparisons** when available

## Pipelines and Backends

CBIcall separates the **analysis pipeline** from the **workflow backend** that
executes it.

| Layer | Current scope |
| --- | --- |
| Out-of-box analysis pipelines | WES, WGS, and mtDNA variant calling |
| Workflow backends | Bash, Snakemake, Nextflow, and Cromwell |
| External provider entries | Registered nf-core pipelines launched through the Nextflow backend |

The native CBIcall WES/WGS/mtDNA pipeline implementations use the CBIcall-provided resource
bundle. Registered nf-core provider entries can be validated and launched without that
bundle; Nextflow/nf-core manages the external pipeline's own test data,
containers, and references.

Use [Included Pipelines](pipelines/overview) for the shipped analyses and [Native Backends](backends/native) for the supported workflow backends.

## Installation at a Glance

| Use case | Method |
| --- | --- |
| Local workstation or server | [PyPI](installation/non-containerized), [Docker](installation/docker), or a source checkout for development |
| HPC cluster | [PyPI](installation/non-containerized), [Apptainer / Singularity](installation/apptainer), or site modules |
| Cloud reproducibility check | [Google Cloud](installation/google-cloud-docker) |
| Development or debugging | [Non-containerized source install](installation/non-containerized) |

## Where to Go Next

| Goal | Page |
| --- | --- |
| Try nf-core without the CBIcall bundle | [Quickstart](usage/quickstart) |
| Run the native shipped test data | [Quickstart](usage/quickstart) |
| Run WES/WGS data | [WES Example](usage/end-to-end-example-wes) |
| Run mtDNA analysis | [mtDNA Example](usage/end-to-end-example-mit) |
| See included pipelines | [Included Pipelines](pipelines/overview) |
| Understand workflow backends | [Native Backends](backends/native) |
| Understand output files | [Outputs](help/outputs) |
| See the system design | [Architecture](technical-details/architecture) |
