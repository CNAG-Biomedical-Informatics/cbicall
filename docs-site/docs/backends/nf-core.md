# nf-core Backend

CBIcall can launch registered nf-core workflows through Nextflow. CBIcall
validates the YAML, pins the registered nf-core release, writes the Nextflow
params/config files, and records the run metadata. The nf-core workflow keeps
its own native output layout and runtime behavior.

nf-core is an external workflow adapter in CBIcall, not one of the shipped
WES/WGS/mtDNA analysis pipelines. Use this page when testing registered nf-core
examples on a workstation or on a cluster.

## How CBIcall Complements nf-core

CBIcall is not a replacement for nf-core. nf-core provides community-maintained
Nextflow workflows; CBIcall adds a project-level execution contract around
selected workflows.

| Aspect | nf-core / Nextflow | CBIcall |
| --- | --- | --- |
| Main role | Community workflow ecosystem | Project-level execution contract |
| Validation | Parameter and samplesheet schemas | Pipeline/mode/genome/backend/resource compatibility |
| Configuration | Samplesheets, params, profiles, configs | User YAML plus controlled registries |
| Outputs | Pipeline-native layout | Deterministic run directory and output inventory |
| Provenance | Nextflow reports and logs | `log.json`, `run-report.json`, fingerprints, VCF hashes |
| Offline/HPC use | Requires staged code, containers, references, and caches | Records the resolved runtime/resource contract |

## Workstation Smoke Test

Use `nf-core/demo` as a lightweight external Nextflow smoke test.

For an x86_64 workstation with Docker:

```yaml
mode:             single
pipeline:         demo
workflow_engine:  nextflow
workflow_version: nf-core
resource:         nf-core-demo-managed-resources-v1
nextflow_profile: test,docker
nextflow_args:    {}
```

Run it from the directory containing `demo.yaml`:

```bash
../../bin/cbicall validate-param -p demo.yaml --no-color
../../bin/cbicall run -p demo.yaml -t 6 --no-color
```

:::tip[Apple Silicon and Linux ARM64]
Some nf-core containers are published primarily for AMD64. On Apple Silicon or
Linux ARM64, Docker may require container emulation. On HPC, prefer the
Apptainer/Singularity profile.
:::

## Sarek Example

Sarek is launched as an external nf-core workflow. CBIcall passes
`nextflow_args` to the generated Nextflow params file.

```yaml
mode:             cohort
pipeline:         sarek
workflow_engine:  nextflow
workflow_version: nf-core
resource:         nf-core-sarek-managed-resources-v1
nextflow_profile: singularity
nextflow_singularity_cache_dir: /path/to/project/nxf-singularity-cache
nextflow_args:
  input: sarek_samplesheet.csv
  genome: GATK.GRCh38
  tools: haplotypecaller
  max_memory: 30.GB
```

Use `max_memory` to cap nf-core/Sarek process memory requests. CBIcall passes it
through as an nf-core parameter; it is not a CBIcall-specific resource field.

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
the nf-core/demo smoke test.

:::info[Nextflow syntax parser]
If Nextflow 26.04 or newer fails while parsing an older nf-core config, run with:

```bash
export NXF_SYNTAX_PARSER=v1
```

This is a Nextflow compatibility setting, not a CBIcall parameter.
:::

## Container Cache

Set a shared Singularity/Apptainer cache when the HPC allows it:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/nxf-singularity-library
export NXF_SINGULARITY_LIBRARYDIR=/shared/nxf-singularity-library
```

Without this variable, Nextflow stores pulled images under the run work
directory, which can make repeated tests slower.

If a shared cache contains an image that your user cannot read, Apptainer may
fail with an error like:

```text
could not open image /shared/nxf-singularity-library/...img
... is not readable by the current user, check permissions
```

Use a user-owned or project-owned cache instead:

```bash
mkdir -p "$HOME/.nextflow-singularity-cache"
export NXF_SINGULARITY_CACHEDIR="$HOME/.nextflow-singularity-cache"
export NXF_SINGULARITY_LIBRARYDIR="$HOME/.nextflow-singularity-cache"
```

or:

```bash
mkdir -p /path/to/project/nxf-singularity-cache
export NXF_SINGULARITY_CACHEDIR=/path/to/project/nxf-singularity-cache
export NXF_SINGULARITY_LIBRARYDIR=/path/to/project/nxf-singularity-cache
```

Then rerun the same CBIcall command. This lets Nextflow pull readable container
images for the current user.

You can also pin this cache in the CBIcall YAML so it is written into the
generated Nextflow config as both the cache and library lookup directory:

```yaml
nextflow_singularity_cache_dir: /path/to/project/nxf-singularity-cache
```

CBIcall writes this as Singularity and Apptainer cache/library settings in
`cbicall_external_nextflow.config`, which is useful when a site-level Nextflow
configuration otherwise points to a shared cache with unreadable images.

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
launcher log. Nextflow reports the task work directory in many errors, but the
task files may not exist if the failure happened before the task was submitted
to the executor.

If the work directory contains task files, inspect them:

```bash
cd work/<hash>/<task>
ls -la
test -f .command.err && cat .command.err
test -f .command.out && cat .command.out
test -f .command.log && cat .command.log
test -f .command.sh && cat .command.sh
```

If `.command.log` and related files are absent, the task usually failed before
execution, for example because requested CPUs exceeded the CPUs allocated by the
scheduler.
