# CBIcall

![CBIcall](/img/cbicall-logo.png)

<p align="center"><em>Reproducible germline variant calling for Illumina DNA sequencing</em></p>

**CBIcall** (**C**NAG **B**iomedical **I**nformatics framework for variant **call**ing) gives users a single command-line entry point for curated WES, WGS, and mtDNA workflows, while keeping run configuration, logs, and outputs traceable.

:::tip[In one sentence]
CBIcall validates a YAML file, resolves the requested workflow, creates a run directory, and launches the matching Bash or Snakemake pipeline.
:::

## What CBIcall Does

CBIcall is an orchestrator. It does not re-implement alignment or variant calling algorithms; it runs curated workflows built from established tools such as BWA, GATK, and MToolBox.

| Responsibility | What it means in practice |
| --- | --- |
| Configuration validation | Invalid engines, genomes, GATK versions, and pipeline/mode combinations fail before execution. |
| Workflow resolution | Requested workflows are resolved through a versioned registry. |
| Run isolation | Each execution gets its own timestamped run directory. |
| Logging and traceability | `log.json` records CLI arguments, resolved configuration, and user parameters. |
| Output organization | Workflow results are written into predictable directories such as `01_bam/`, `02_varcall/`, `03_stats/`, and `logs/`. |

## Supported Workflows

| Input type | Start with | Main output |
| --- | --- | --- |
| WES FASTQs | `pipeline: wes`, `mode: single` | Per-sample QC VCF and gVCF |
| WGS FASTQs | `pipeline: wgs`, `mode: single` | Per-sample QC VCF and gVCF |
| Existing gVCFs | `mode: cohort` | Joint cohort QC VCF |
| WES/WGS BAMs | `pipeline: mit` | mtDNA prioritized variant report and browser output |

## Compatibility

| Pipeline | Mode | Genome | GATK version | Notes |
|---------|------|--------|--------------|------|
| WES | `single` | `b37` | `gatk-3.5`, `gatk-4.6` | Supported |
| WES | `cohort` | `b37` | `gatk-3.5`, `gatk-4.6` | Supported |
| WGS | `single` | `b37`, `hg38` | `gatk-4.6` | Supported |
| WGS | `cohort` | `b37`, `hg38` | `gatk-4.6` | Supported |
| mtDNA | `single` | `rsrs` | `gatk-3.5` | x86_64 only |
| mtDNA | `cohort` | `rsrs` | `gatk-3.5` | x86_64 only |

:::info[Workflow engines]
`bash` is the broadest supported engine. `snakemake` is supported for GATK 4.6 WES/WGS workflows. mtDNA workflows currently use Bash.
:::

## How a Run Is Traced

Every run creates:

- a generated run identifier
- a dedicated output directory
- a workflow log
- `log.json` with CLI arguments, resolved configuration, and user parameters
- pipeline-specific outputs

This makes completed runs easier to inspect, compare, debug, and reproduce.

## Where to Go Next

| Goal | Page |
| --- | --- |
| Decide what to run | [Choose Your Path](usage/choose-your-path) |
| Run the shipped test data | [Quickstart](usage/quickstart) |
| Configure a real run | [Configuration Reference](help/configuration-reference) |
| Understand output files | [Outputs](help/outputs) |
| See the system design | [Architecture](technical-details/architecture) |
