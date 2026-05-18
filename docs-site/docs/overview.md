# CBIcall

![CBIcall](/img/cbicall-logo.png)

<p align="center"><em>Reproducible germline variant calling for Illumina DNA sequencing</em></p>

**CBIcall** (**C**NAG **B**iomedical **I**nformatics framework for variant **call**ing) gives users a single command-line entry point for curated WES, WGS, and mtDNA workflows, while keeping run configuration, logs, and outputs traceable.

:::tip[In one sentence]
CBIcall validates a parameters YAML file, resolves the requested workflow, creates a run directory, and launches the matching Bash or Snakemake pipeline.
:::

## What CBIcall Does

CBIcall is an orchestrator. It does not re-implement alignment or variant
calling algorithms; it runs curated workflows built from established tools such
as BWA, GATK, and MToolBox.

- validates the **parameters YAML** before launch
- dispatches **Bash** or **Snakemake** workflows
- records **logs, provenance, run reports, and output fingerprints** when available

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
| Understand output files | [Outputs](help/outputs) |
| See the system design | [Architecture](technical-details/architecture) |
