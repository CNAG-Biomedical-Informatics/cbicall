# CBIcall

![CBIcall](/img/cbicall-logo.png)

<p align="center"><em>Reproducible germline variant calling for Illumina DNA sequencing</em></p>

**CBIcall** (**C**NAG **B**iomedical **I**nformatics framework for variant **call**ing) gives users a single command-line entry point for curated WES, WGS, and mtDNA workflows, while keeping run configuration, logs, and outputs traceable.

:::tip[In one sentence]
CBIcall validates a parameters YAML file, resolves the requested workflow, creates a run directory, and launches the matching Bash, Snakemake, or Nextflow pipeline.
:::

## What CBIcall Does

CBIcall is an orchestrator. It does not re-implement alignment or variant
calling algorithms; it runs curated workflows built from established tools such
as BWA, GATK, and MToolBox.

- validates the **parameters YAML** before launch
- dispatches **Bash**, **Snakemake**, or **Nextflow** workflows
- records **logs, provenance, run reports, and output fingerprints** when available

## Pipelines and Backends

CBIcall separates the **analysis pipeline** from the **workflow backend** that
executes it.

| Layer | Current scope |
| --- | --- |
| Out-of-box analysis pipelines | WES, WGS, and mtDNA variant calling |
| Workflow backends | Bash, Snakemake, and Nextflow |
| External workflows | Registered nf-core workflows launched through the Nextflow backend |

Use [Included Pipelines](pipelines/overview) for the shipped analyses and
[Native Backends](backends/native) for how those workflows are dispatched.

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
