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
GATK, MToolBox, Snakemake, Nextflow, and selected nf-core pipelines.

- validates the **parameters YAML** and compatibility contract before launch
- resolves **Bash**, **Snakemake**, or **Nextflow** workflow backends
- records **logs, provenance, run reports, output fingerprints, and run comparisons** when available

## Pipelines and Backends

CBIcall separates the **analysis pipeline** from the **workflow backend** that
executes it.

| Layer | Current scope |
| --- | --- |
| Out-of-box analysis pipelines | WES, WGS, and mtDNA variant calling |
| Workflow backends | Bash, Snakemake, and Nextflow |
| External workflows | Registered nf-core workflows launched through the Nextflow backend |

Use [Included Pipelines](pipelines/overview) for the shipped analyses and
[Native Backends](backends/native) for the supported workflow backends.

## Installation at a Glance

| Use case | Method |
| --- | --- |
| Local workstation or server | [Docker](installation/docker) |
| HPC cluster | [Apptainer / Singularity](installation/apptainer) |
| Development or debugging | [From Source](installation/non-containerized) |
| Cloud reproducibility check | [Google Cloud](installation/google-cloud-docker) |

## Where to Go Next

| Goal | Page |
| --- | --- |
| Run the shipped test data | [Quickstart](usage/quickstart) |
| Run WES/WGS data | [WES Example](usage/end-to-end-example-wes) |
| Run mtDNA analysis | [mtDNA Example](usage/end-to-end-example-mit) |
| See included pipelines | [Included Pipelines](pipelines/overview) |
| Understand workflow backends | [Native Backends](backends/native) |
| Understand output files | [Outputs](help/outputs) |
| See the system design | [Architecture](technical-details/architecture) |
