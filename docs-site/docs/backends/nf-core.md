# nf-core Provider

CBIcall can launch registered [nf-core](https://nf-co.re/) pipelines through
[Nextflow](https://www.nextflow.io/docs/latest/). CBIcall validates the YAML,
pins the registered nf-core release, writes the Nextflow params/config files,
and records the run metadata. The nf-core pipeline keeps its own native output
layout and runtime behavior.

In CBIcall terms, nf-core is a **workflow provider** and Nextflow is the
**workflow backend**. Registered nf-core entries are launched through the
Nextflow backend, but they are not one of the shipped CBIcall WES/WGS/mtDNA
native CBIcall implementations. Use this page when testing registered nf-core
examples on a workstation or on a cluster.

:::note[Run directory location]
nf-core run directories are created where `cbicall run` is launched,
for example `cbicall_nextflow_nf-core_sarek_cohort_GATK.GRCh38_<run-id>/`.
This differs from native CBIcall WES/WGS/mtDNA pipelines, whose run directories
are created under the discovered sample/input directory.
:::

## Profiles and Container Runtimes

CBIcall passes `nfcore_profile` to nf-core as the Nextflow `-profile` value.
Profiles are comma-separated and combined by Nextflow; each profile adds a
different part of the runtime configuration.

| Profile | What it does |
| --- | --- |
| `test` | Uses the nf-core pipeline's built-in test data and smoke-test settings. It is not a software runtime. |
| `docker` | Runs tasks with Docker. Images are pulled into the local Docker image store. |
| `singularity` | Runs tasks with Singularity/Apptainer. Images are pulled into the Nextflow/Singularity cache. |
| `conda` | Builds Conda environments when the workflow supports it. Useful in some setups, but not the recommended CBIcall HPC path. |

For smoke tests, combine `test` with a runtime profile, for example
`test,docker` on a Docker workstation or `test,singularity` on HPC. If
`nfcore_parameters.input` is set, that explicit input overrides the test input
provided by the nf-core `test` profile.

For Singularity/Apptainer, image placement is controlled by the environment
variables `NXF_SINGULARITY_CACHEDIR` and `NXF_SINGULARITY_LIBRARYDIR`, or by
CBIcall's `nfcore_singularity_cache_dir` setting.

## How CBIcall Complements nf-core

CBIcall is not a replacement for nf-core. nf-core provides community-maintained
pipelines run by Nextflow; CBIcall adds a project-level execution contract
around selected registered provider entries.

| Aspect | nf-core / Nextflow | CBIcall |
| --- | --- | --- |
| Main role | Community pipeline ecosystem | Project-level execution contract |
| Validation | Parameter and samplesheet schemas | Pipeline/mode/genome/backend/resource compatibility |
| Configuration | Samplesheets, params, profiles, configs | User YAML plus controlled registries |
| Outputs | Pipeline-native layout | Deterministic run directory and output inventory |
| Provenance | Nextflow reports and logs | `log.json`, `run-report.json`, fingerprints, VCF hashes |
| Offline/HPC use | Requires staged code, containers, references, and caches | Records the resolved runtime/resource contract |

## Demo Smoke Test

Use [`nf-core/demo`](https://nf-co.re/demo) as a lightweight nf-core provider
smoke test through the Nextflow backend. The checked-in
`examples/input/nf-core-demo.yaml` uses `test,singularity` because that is the
practical profile on HPC systems where Docker is not available.

For an HPC or Apptainer/Singularity test:

```yaml
mode:             single
pipeline:         demo
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-demo-managed-resources-v1
nfcore_profile: test,singularity
```

Run it from the directory containing `nf-core-demo.yaml`:

```bash
../../bin/cbicall validate-parameters -p nf-core-demo.yaml --no-color
../../bin/cbicall run -p nf-core-demo.yaml -t 4 --no-color
```

:::tip[Workstation Docker alternative]
On an x86_64 workstation with Docker installed, change `nfcore_profile` to
`test,docker`. On HPC, prefer `test,singularity`; Docker is usually unavailable
inside cluster jobs.
:::

:::note[Default Singularity cache]
For nf-core runs using `singularity`, if neither `NXF_SINGULARITY_CACHEDIR` nor
`nfcore_singularity_cache_dir` is configured, Nextflow stores pulled images
inside the run directory, usually under `work/singularity`, and may print:

```text
WARN: Singularity cache directory has not been defined -- Remote image will be stored in the path: <run-dir>/work/singularity
```

For repeated or production runs, configure the cache explicitly through the
shell, scheduler environment, or `nfcore_singularity_cache_dir`.
:::

:::note[CPU requirement]
`nf-core/demo` includes steps that request more than one CPU. Use at least
`-t 4` and request matching scheduler CPUs. The tiny Sarek chr22 smoke test can
run with `-t 1`, but demo generally cannot.
:::

## Sarek Examples

[`nf-core/sarek`](https://nf-co.re/sarek) is launched as an nf-core provider
entry through the Nextflow backend. CBIcall passes `nfcore_parameters` to the
generated Nextflow params file.

### Built-In nf-core Test

The smallest Sarek check uses nf-core's own `test` profile. It verifies that
CBIcall can launch the registered Sarek release, write the params/config files,
pull containers, parse Nextflow reports, and record software versions.

```yaml
mode:             cohort
pipeline:         sarek
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-sarek-managed-resources-v1
nfcore_profile: test,singularity
```

This is an execution smoke test, not the HaplotypeCaller output example. The
Sarek `test` profile supplies its own inputs and defaults to Strelka, so the
CBIcall canonical HaplotypeCaller VCF pattern is expected to be missing.

### CBIcall Test Data

Use `examples/input/nf-core-sarek.yaml` when you want the CBIcall test-data
example. It points Sarek to `examples/input/sarek_samplesheet.csv`, uses the
tiny GRCh38 chr22 interval file, and asks for HaplotypeCaller output.

```yaml
mode:             cohort
pipeline:         sarek
workflow_backend:  nextflow
workflow_provider: nf-core
resource:         nf-core-sarek-managed-resources-v1
nfcore_profile: singularity
# nfcore_singularity_cache_dir: nxf-singularity-cache
nfcore_parameters:
  input: sarek_samplesheet.csv
  genome: GATK.GRCh38
  tools: haplotypecaller
  skip_tools: haplotypecaller_filter
  wes: true
  intervals: ../../workflows/nf-core/sarek/grch38_chr22_test.bed
  max_memory: 30.GB
```

Use `max_memory` to cap nf-core/Sarek process memory requests. CBIcall writes it
to the generated nf-core parameters file and mirrors it in Nextflow
`process.resourceLimits`, so the scheduler sees the same memory ceiling.
The example interval file is a tiny GRCh38 chr22 BED stored with the registered
Sarek support files.
The example also skips Sarek's HaplotypeCaller filtering step because tiny
single-chromosome smoke tests may not contain enough overlapping resource
variants for GATK `FilterVariantTranches`.

On HPC, use `nfcore_singularity_cache_dir` when the default shared container
cache contains unreadable images or when you want the cache path recorded in the
generated Nextflow config. CBIcall inherits any `NXF_*` variables already set by
the shell or scheduler bootstrap.

The Sarek registry entry also declares its canonical HaplotypeCaller VCF output:

```text
sarek/variant_calling/haplotypecaller/*/*.haplotypecaller.vcf.gz
```

When this file exists, `run-report.json` records a normalized VCF hash so
`cbicall compare-runs` can audit repeated Sarek runs.

![Screenshot of a CBIcall HTML run comparison report for repeated nf-core/Sarek runs, showing the combined status KPI and legend row plus matching canonical HaplotypeCaller VCF fingerprints.](/img/compare-runs-nfcore-html-report.png)

For the germline HaplotypeCaller smoke test, the Sarek samplesheet can use the
minimal FASTQ layout shown in `examples/input/sarek_samplesheet.csv`:

```csv
patient,sample,lane,fastq_1,fastq_2
CNAG99901P_ex,CNAG99901P_ex,L001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
```

Keep the header. Without it, Sarek may parse the first sample row as column
names and fail before launching the expected pipeline steps. For the full
samplesheet schema and advanced use cases, use the
[nf-core/Sarek documentation](https://nf-co.re/sarek).

## HPC Notes

On a cluster, load the runtime used by the selected nf-core profile. Check
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

When running from a source checkout, make sure the SLURM script loads the Python
runtime needed by CBIcall before calling `bin/cbicall`. For example:

```bash
module load Python/3.10.8-GCCcore-12.2.0
export PYTHONPATH=/software/biomed/cbi_py3/lib/python3.10/site-packages:$PYTHONPATH
```

The bundled SLURM examples set this directly in the generated job script rather
than requiring a separate sourced bootstrap.

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

../../bin/cbicall run -p nf-core-demo.yaml -t "$SLURM_CPUS_PER_TASK"
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
`nfcore_singularity_cache_dir`. For general container and offline setup, see the
[nf-core documentation](https://nf-co.re/docs/usage/getting_started/configuration)
and the [Nextflow container documentation](https://www.nextflow.io/docs/latest/container.html).
:::

## Logs

CBIcall writes the main Nextflow launcher log in the run directory:

```text
cbicall_nextflow_nf-core_<pipeline>_<mode>_<display-genome>_<run_id>/
  nf-core_<pipeline>_<mode>.log
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
for generic backend-level debugging.
:::
