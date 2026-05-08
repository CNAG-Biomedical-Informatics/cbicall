# Choose Your Path

Use this page to pick the right installation method and workflow before opening the detailed examples.

:::tip[Short answer]
Use Docker for local runs, Apptainer on HPC, `pipeline: wes` or `pipeline: wgs` for FASTQs, and `pipeline: mit` only after WES/WGS BAMs already exist.
:::

## Fast Recommendations

| Situation | Start here | Then run |
|---|---|---|
| Local workstation or server with Docker | [Docker](../installation/docker) | [Quickstart](quickstart) |
| HPC cluster or Slurm environment | [Apptainer](../installation/apptainer) | WES/WGS or mtDNA examples |
| Development, debugging, or workflow editing | [Non-containerized](../installation/non-containerized) | Any supported workflow |
| First time testing CBIcall | [Quickstart](quickstart) | WES single-sample smoke test |
| Real WES or WGS data from FASTQ | [WES example](end-to-end-example-wes) | `pipeline: wes` or `pipeline: wgs` |
| mtDNA analysis | [mtDNA example](end-to-end-example-mit) | `pipeline: mit` after WES/WGS BAMs exist |

## 1. Pick an Installation Route

:::tip[Most users]

Use **Docker** on a local workstation or server. Use **Apptainer** on HPC.

:::

### Docker

Choose Docker when:

- you are running locally or on a server where Docker is available
- you want the fastest reproducible setup
- you do not need direct scheduler integration

Go to [Docker installation](../installation/docker).

### Apptainer / Singularity

Choose Apptainer when:

- you are on an HPC cluster
- Docker is not available or not permitted
- you need a container workflow compatible with Slurm or another scheduler

Go to [HPC with Apptainer / Singularity](../installation/apptainer).

### Non-containerized

Choose a source installation when:

- you are developing or debugging CBIcall
- you want to inspect or modify workflow scripts
- you already manage Python, Java, Snakemake, and external tools on the host

Go to [Non-containerized installation](../installation/non-containerized).

## 2. Pick a Pipeline

| Pipeline | Input | Main use | Notes |
|---|---|---|---|
| `wes` | WES FASTQs | Exome germline calling | Supports single and cohort modes |
| `wgs` | WGS FASTQs | Genome-wide germline calling | Supports single and cohort modes with GATK 4.6 |
| `mit` | BAMs from previous WES/WGS runs | mtDNA calling with MToolBox | Does not start from FASTQ |

Use [WES/WGS single-sample](../pipelines/wes-wgs-single) first when starting from FASTQ.

Use [WES/WGS cohort](../pipelines/wes-wgs-cohort) when per-sample gVCFs already exist and you want joint genotyping.

Use [mtDNA pipelines](../pipelines/mtdna) only after the required WES/WGS BAM outputs are available.

## 3. Pick `single` or `cohort`

### `mode: single`

Use this when:

- processing one individual
- creating a per-sample VCF or gVCF
- preparing inputs for a later cohort run

This is the normal first step for WES/WGS.

### `mode: cohort`

Use this when:

- joint genotyping multiple WES/WGS samples
- combining existing per-sample gVCFs
- analyzing mtDNA across a family or cohort after WES/WGS BAMs exist

:::warning

For WES/WGS, run `single` first for each sample. Cohort mode expects per-sample gVCF inputs.

:::

## 4. Pick a Workflow Engine

| Engine | Use when | Limitations |
|---|---|---|
| `bash` | You want the broadest supported path or mtDNA workflows | No partial workflow-rule starts |
| `snakemake` | You want Snakemake orchestration for WES/WGS GATK 4.6 workflows | Not supported for mtDNA or GATK 3.5 |

When unsure, start with:

```yaml
workflow_engine: bash
gatk_version: gatk-4.6
```

## Recommended Next Pages

- [Quickstart](quickstart): shortest validated local test
- [End-to-end Example: WES](end-to-end-example-wes): realistic WES run from FASTQ to outputs
- [End-to-end Example: mtDNA](end-to-end-example-mit): mtDNA workflow after WES/WGS BAM generation
- [Configuration Reference](../help/configuration-reference): YAML keys and supported combinations
- [General Usage](/docs/usage): command syntax and execution patterns
