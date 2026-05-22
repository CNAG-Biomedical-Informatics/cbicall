# CBIcall

![CBIcall](/img/cbicall-logo.png)

<p align="center"><em>Configuration-driven framework for reproducible cohort-scale variant calling</em></p>

**CBIcall** (**C**NAG **B**iomedical **I**nformatics framework for variant **call**ing) is a configuration-driven framework for reproducible variant calling in large sequencing cohorts.

:::tip[In one sentence]
CBIcall validates a parameters YAML file, resolves an approved workflow backend, creates a deterministic run directory, and records structured provenance for audit and run comparison.
:::

## What CBIcall Does

CBIcall does not re-implement alignment or variant-calling algorithms; those
steps remain in curated workflows built from established tools such as BWA,
GATK, MToolBox, Snakemake, Nextflow, Cromwell, and selected nf-core pipelines.

- validates the **parameters YAML** and compatibility contract before launch
- resolves **Bash**, **Snakemake**, **Nextflow**, or **Cromwell** workflow backends
- records **logs, provenance, run reports, output fingerprints, and run comparisons** when available

## Pipelines and Backends

CBIcall separates the **analysis pipeline** from the **workflow backend** that
executes it.

| Layer | Current scope |
| --- | --- |
| Out-of-box analysis pipelines | WES, WGS, and mtDNA variant calling |
| Workflow backends | Bash, Snakemake, Nextflow, and Cromwell |
| External workflows | Registered nf-core workflows launched through the Nextflow backend |

The native CBIcall WES/WGS/mtDNA workflows use the CBIcall-provided resource
bundle. Registered nf-core workflows can be validated and launched without that
bundle; Nextflow/nf-core manages the external workflow's own test data,
containers, and references.

Use [Included Pipelines](pipelines/overview) for the shipped analyses and
[Native Backends](backends/native) for the supported workflow backends.

## Installation at a Glance

| Use case | Method |
| --- | --- |
| Local workstation or server | [Docker](installation/docker) |
| HPC cluster | [Apptainer / Singularity](installation/apptainer) |
| Cloud reproducibility check | [Google Cloud](installation/google-cloud-docker) |
| Development or debugging | [From Source](installation/non-containerized) |

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
