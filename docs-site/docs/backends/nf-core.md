# nf-core Backend

CBIcall can launch registered [nf-core](https://nf-co.re/) workflows through
[Nextflow](https://www.nextflow.io/docs/latest/). CBIcall validates the YAML,
pins the registered nf-core release, writes the Nextflow params/config files,
and records the run metadata. The nf-core workflow keeps its own native output
layout and runtime behavior.

nf-core is an external workflow adapter in CBIcall, not one of the shipped
WES/WGS/mtDNA analysis pipelines. Use this page when testing registered nf-core
examples on a workstation or on a cluster.

:::note[Run directory location]
External nf-core run directories are created where `cbicall run` is launched,
for example `cbicall_nextflow_sarek_cohort_external_nf-core_<run-id>/`.
This differs from native CBIcall WES/WGS/mtDNA pipelines, whose run directories
are created under the discovered sample/input directory.
:::

## How CBIcall Complements nf-core

CBIcall is not a replacement for nf-core. nf-core provides community-maintained
Nextflow workflows; CBIcall adds a project-level execution contract around
selected registered workflows.

| Aspect | nf-core / Nextflow | CBIcall |
| --- | --- | --- |
| Main role | Community workflow ecosystem | Project-level execution contract |
| Validation | Parameter and samplesheet schemas | Pipeline/mode/genome/backend/resource compatibility |
| Configuration | Samplesheets, params, profiles, configs | User YAML plus controlled registries |
| Outputs | Pipeline-native layout | Deterministic run directory and output inventory |
| Provenance | Nextflow reports and logs | `log.json`, `run-report.json`, fingerprints, VCF hashes |
| Offline/HPC use | Requires staged code, containers, references, and caches | Records the resolved runtime/resource contract |

## Demo Smoke Test

Use [`nf-core/demo`](https://nf-co.re/demo) as a lightweight external
Nextflow smoke test. The checked-in `examples/input/demo.yaml` uses
`test,singularity` because that is the practical profile on HPC systems where
Docker is not available.

For an HPC or Apptainer/Singularity test:

```yaml
mode:             single
pipeline:         demo
workflow_engine:  nextflow
workflow_version: nf-core
resource:         nf-core-demo-managed-resources-v1
nextflow_profile: test,singularity
nextflow_args:    {}
```

Run it from the directory containing `demo.yaml`:

```bash
../../bin/cbicall validate-param -p demo.yaml --no-color
../../bin/cbicall run -p demo.yaml -t 4 --no-color
```

:::note[How nf-core profiles combine]
The `test` profile provides nf-core's built-in test inputs and smoke-test
settings. It does not select the software runtime by itself. Combine it with a
runtime profile such as `singularity` or `docker`: for example,
`test,singularity`. If `nextflow_args.input` is set, that explicit input
overrides the input supplied by the `test` profile. See the
[nf-core profile documentation](https://nf-co.re/docs/usage/getting_started/configuration)
for generic profile behavior.
:::

:::tip[Workstation Docker alternative]
On an x86_64 workstation with Docker, change `nextflow_profile` to
`test,docker`. On HPC, prefer `test,singularity`; Docker is usually unavailable
inside cluster jobs.
:::

:::note[CPU requirement]
`nf-core/demo` includes steps that request more than one CPU. Use at least
`-t 4` and request matching scheduler CPUs. The tiny Sarek chr22 smoke test can
run with `-t 1`, but demo generally cannot.
:::

## Sarek Example

[`nf-core/sarek`](https://nf-co.re/sarek) is launched as an external nf-core
workflow. CBIcall passes `nextflow_args` to the generated Nextflow params file.

```yaml
mode:             cohort
pipeline:         sarek
workflow_engine:  nextflow
workflow_version: nf-core
resource:         nf-core-sarek-managed-resources-v1
nextflow_profile: singularity
nextflow_singularity_cache_dir: nxf-singularity-cache
nextflow_args:
  input: sarek_samplesheet.csv
  genome: GATK.GRCh38
  tools: haplotypecaller
  skip_tools: haplotypecaller_filter
  wes: true
  intervals: ../../workflows/nf-core/sarek/grch38_chr22_test.bed
  max_memory: 30.GB
```

Use `max_memory` to cap nf-core/Sarek process memory requests. CBIcall writes it
to the generated nf-core params file and mirrors it in Nextflow
`process.resourceLimits`, so the scheduler sees the same memory ceiling.
The example interval file is a tiny GRCh38 chr22 BED stored with the registered
Sarek workflow support files.
The example also skips Sarek's HaplotypeCaller filtering step because tiny
single-chromosome smoke tests may not contain enough overlapping resource
variants for GATK `FilterVariantTranches`.

On HPC, use `nextflow_singularity_cache_dir` when the default shared container
cache contains unreadable images or when you want the cache path recorded in the
generated Nextflow config. CBIcall inherits any `NXF_*` variables already set by
the shell or scheduler bootstrap.

The Sarek registry entry also declares its canonical HaplotypeCaller VCF output:

```text
sarek/variant_calling/haplotypecaller/*/*.haplotypecaller.vcf.gz
```

When this file exists, `run-report.json` records a normalized VCF hash so
`cbicall compare-runs` can audit repeated Sarek runs.

![Screenshot of a CBIcall HTML run comparison report for repeated nf-core/Sarek runs, showing matching framework, pipeline, resource, and canonical HaplotypeCaller VCF fingerprints.](/img/compare-runs-nfcore-html-report.png)

For germline HaplotypeCaller testing, the Sarek samplesheet must include a
header and mark the sample as normal:

```csv
patient,sample,lane,fastq_1,fastq_2,status
CNAG99901P_ex,CNAG99901P_ex,L001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,0
```

The `status` column uses Sarek's convention:

| Value | Meaning |
| --- | --- |
| `0` | Normal |
| `1` | Tumor |

If the header is missing, Sarek may interpret the sample incorrectly and report
that only tumor samples were provided.

## HPC Notes

On a cluster, load the runtime used by the selected Nextflow profile. Check
available modules:

```bash
module avail Nextflow/25.10.2
module avail apptainer
module avail singularity
```

For Apptainer/Singularity:

```bash
module load Nextflow/25.10.2
module load apptainer   # or: module load singularity
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity/cache
export NXF_SINGULARITY_LIBRARYDIR=/path/to/singularity/library
```

At CNAG, source the example bootstrap before running CBIcall from a source
checkout:

```bash
source /software/biomed/cbicall/examples/scripts/cnag_hpc_cbicall_env.sh
```

This loads the Python runtime required by CBIcall before the Python entrypoint
starts. It also tries to load Nextflow for nf-core workflows.

:::tip[Container profile name]
Many nf-core configs still use the profile name `singularity`, even when the
cluster executable is Apptainer. If `test,apptainer` is not accepted by the
pipeline, use `test,singularity`.
:::

The thread value passed to CBIcall limits the workflow, but it does not allocate
CPUs from the scheduler. Request enough CPUs in the batch job and pass the same
value to CBIcall:

```bash
#SBATCH --cpus-per-task=6

../../bin/cbicall run -p demo.yaml -t "$SLURM_CPUS_PER_TASK"
```

If Nextflow reports:

```text
Process requirement exceeds available CPUs -- req: 4; avail: 1
```

the scheduler job was started with too few CPUs. Request at least four CPUs for
the nf-core/demo smoke test. The tiny Sarek chr22 smoke test can run with
`-t 1`, but demo generally needs more than one CPU.

:::info[Nextflow syntax parser]
If Nextflow 26.04 or newer fails while parsing an older nf-core config, run with:

```bash
export NXF_SYNTAX_PARSER=v1
```

This is a Nextflow compatibility setting, not a CBIcall parameter.
:::

:::tip[Container cache]
If Apptainer/Singularity cannot read an image from a shared Nextflow cache, use a
user-owned or project-owned cache and point CBIcall to it with
`nextflow_singularity_cache_dir`. For general container and offline setup, see
the [nf-core documentation](https://nf-co.re/docs/usage/getting_started/configuration)
and the [Nextflow container documentation](https://www.nextflow.io/docs/latest/container.html).
:::

## Logs

CBIcall writes the main Nextflow launcher log in the run directory:

```text
cbicall_nextflow_<pipeline>_<mode>_external_nf-core_<run_id>/
  nextflow_<pipeline>_<mode>_external_nf-core.log
  .nextflow.log
  cbicall_external_nextflow.params.yaml
  cbicall_external_nextflow.config
```

For task-level failures, start with `.nextflow.log` and the main CBIcall
launcher log. Nextflow often reports the task work directory in the error text.

:::tip[Debugging failed nf-core tasks]
When a task reaches execution, Nextflow usually writes `.command.*` files under
`work/<hash>/<task>/`. A compact first check is:

```bash
cd work/<hash>/<task>
ls -la .command*
```

If those files are absent, the task likely failed before execution, commonly
because the scheduler allocation was smaller than the task request. See the
[Nextflow troubleshooting documentation](https://www.nextflow.io/docs/latest/troubleshooting.html)
for generic engine-level debugging.
:::
