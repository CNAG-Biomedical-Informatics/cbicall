# Adding a Pipeline

Start by deciding whether the workflow will be **CBIcall-native** or
**external**. This is the most important design choice.

| Integration level | Meaning | What CBIcall can audit |
| --- | --- | --- |
| **CBIcall-native** | The workflow writes the expected CBIcall run layout inside the generated `cbicall_*` directory. | Full run report, output inventory, CBIcall logs, final-output fingerprints, and `compare-runs`. |
| **External** | CBIcall launches and audits the workflow, but the workflow keeps its upstream output layout. | Launch contract, parameters, resources, backend logs, output inventory, and declared `canonical_outputs` when configured. |

In short: **native means CBIcall output contract**, not workflow language. A
Snakemake, Nextflow, or Cromwell workflow can be native if it writes or promotes
outputs into the CBIcall layout. An nf-core workflow is external because it keeps
the nf-core output layout. nf-core is the first external provider supported by
CBIcall, but the same registry model can be used for other workflow providers
that expose a stable source, version, launch command, and final-output patterns.

## Recommended Path

For a new CBIcall-native workflow, usually do this:

1. Add the workflow entrypoint under `workflows/<backend>/<software_stack>/`.
2. Make the workflow write the CBIcall output layout: `01_bam/`, `02_varcall/`,
   `03_stats/`, `logs/`, and any expected final files.
3. Register it in `workflows/registry/cbicall-workflow-registry.yaml`.
4. Add or update resource-catalog compatibility if it needs a resource bundle.
5. Run `bin/cbicall validate-registry` and a small test run.

For an external workflow, usually do this:

1. Register the upstream workflow source and pinned release.
2. Decide which runtime parameters CBIcall should pass through.
3. Add `canonical_outputs` for final deliverables that should be fingerprinted.
4. Add a lightweight resource entry if the upstream workflow manages references
   and containers itself.
5. Run a smoke test and inspect `run-report.json`.

## What You Need To Define

| Decision | Examples | Why it matters |
| --- | --- | --- |
| Integration level | native, external | Determines whether the workflow must follow the CBIcall output contract. |
| Backend | `bash`, `snakemake`, `nextflow`, `cromwell` | Determines how CBIcall launches the workflow. |
| Provider | `cbicall`, `nf-core` | Identifies whether the workflow is maintained by CBIcall or an external ecosystem. |
| Software stack | `gatk-3.5`, `gatk-4.6`, `nf-core` | Selects implementation files or an external workflow source. |
| Pipeline and mode | `wes`/`single`, `wgs`/`cohort` | Defines the user-facing YAML values and registry path. |
| Inputs | `input_dir`, `sample_map`, backend parameters | Defines how users provide samples and workflow-specific options. |

Example native YAML:

```yaml
mode: single
pipeline: wes
workflow_provider: cbicall
workflow_backend: bash
software_stack: gatk-4.6
genome: b37
input_dir: SAMPLE01
```

Example external nf-core YAML:

```yaml
mode: cohort
pipeline: sarek
workflow_provider: nf-core
workflow_backend: nextflow
resource: nf-core-sarek-managed-resources-v1
nfcore_profile: singularity
nfcore_parameters:
  input: sarek_samplesheet.csv
  genome: GATK.GRCh38
  tools: haplotypecaller
```

## How Execution Works

At runtime, CBIcall reads the parameters YAML, validates compatibility rules,
resolves the workflow registry entry, creates a run directory, launches the
selected backend from that directory, and writes `log.json`,
`cbicall-execution-contract.json`, and `run-report.json`.

The main implementation files are:

| File | Responsibility |
| --- | --- |
| `src/cbicall/config.py` | Parameter defaults and semantic validation. |
| `src/cbicall/workflow_registry.py` | Registry loading and workflow resolution. |
| `src/cbicall/execution.py` | Backend-specific launch logic. |
| `workflows/registry/cbicall-workflow-registry.yaml` | Developer-facing workflow map. |
| `workflows/schema/cbicall-workflow-registry.schema.json` | Registry schema. |
| `resources/cbicall-resource-catalog.json` | Resource-bundle compatibility and identity. |

After editing the registry, run:

```bash
bin/cbicall validate-registry
```

## 1. Add the Workflow Entrypoint

For CBIcall-native workflows, entrypoint filenames should make the pipeline and
mode obvious. Use the backend extension that matches the workflow language:

```text
{pipeline}_{mode}.sh   # Bash
{pipeline}_{mode}.smk  # Snakemake
{pipeline}_{mode}.nf   # Nextflow
{pipeline}_{mode}.wdl  # Cromwell/WDL
```

### Bash Layout

```text
workflows/bash/gatk-4.6/
  env.sh
  mypipe_single.sh
  mypipe_cohort.sh
```

Bash workflows are executed directly. CBIcall sets `GENOME` in the environment and launches the script from inside the generated run directory.

:::info[Entrypoint location]
CBIcall does not copy Bash workflow scripts into the run directory. It launches the registered script from `workflows/bash/...` while setting the process working directory to the generated run directory.

This keeps workflow code centralized and keeps helper paths such as `env.sh` stable through `BASH_SOURCE[0]`. The tradeoff is that workflow `.sh` files should not be edited while jobs are running.
:::

Minimal Bash example:

```bash
#!/usr/bin/env bash
set -eu

BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNDIR="$(pwd)"

source "$BINDIR/env.sh"

mkdir -p logs results

echo "Running mypipe in $RUNDIR" | tee logs/mypipe.log

for R1 in ../*_R1_*fastq.gz; do
  R2="${R1/_R1_/_R2_}"
  echo "Pair: $R1 $R2" >> logs/mypipe.log
done

echo "done" > results/mypipe.done
```

Make it executable:

```bash
chmod +x workflows/bash/gatk-4.6/mypipe_single.sh
```

### Snakemake Layout

```text
workflows/snakemake/gatk-4.6/
  config.yaml
  mypipe_single.smk
  mypipe_cohort.smk
```

CBIcall launches Snakemake with:

- the resolved Snakefile through `-s`
- the shared config file through `--configfile`
- `genome` through `--config`
- `sample_map` and `workspace` for cohort workflows when needed
- `snakemake_parameters.target` as the Snakemake target when targeted execution is requested

### Native Nextflow Layout

Native Nextflow workflows live under the selected software stack:

```text
workflows/nextflow/gatk-4.6/
  config.yaml
  wes_single.nf
  wgs_single.nf -> wes_single.nf
  wes_cohort.nf
  wgs_cohort.nf -> wes_cohort.nf
```

CBIcall launches the resolved `.nf` file with CBIcall-managed parameters,
including genome, threads, sample map or input directory, helper scripts, and
workspace paths. To be CBIcall-native, the workflow must write or promote final
outputs into the standard CBIcall run layout. Native Nextflow registry entries
also require the helper paths used by post-processing and run comparison.

### Native Cromwell Layout

Native Cromwell workflows live under the selected software stack:

```text
workflows/cromwell/gatk-4.6/
  config.yaml -> ../../snakemake/gatk-4.6/config.yaml
  wes_single.wdl
  wgs_single.wdl -> wes_single.wdl
  wes_cohort.wdl
  wgs_cohort.wdl -> wes_cohort.wdl
```

If `WOMTOOL_JAR` is set, `bin/cbicall validate-registry` validates registered
Cromwell WDL syntax before any run is launched.

CBIcall writes `cbicall_cromwell.inputs.json`,
`cbicall_cromwell.options.json`, and `cbicall_cromwell.metadata.json` in the run
directory. The WDL receives absolute paths for inputs, tools, references, and
helper scripts. Final outputs are promoted back into the standard CBIcall
`01_bam/`, `02_varcall/`, `03_stats/`, and `logs/` layout so run reports and
`compare-runs` work like the other native backends.

### External Workflow Integration

External workflows can be registered when CBIcall should validate, launch, log,
and report them without adapting their internal output layout to the CBIcall
folder contract. The currently bundled external provider examples are nf-core
workflows launched through the Nextflow backend. The concept is not limited to
nf-core; another provider can be added when CBIcall can declare how to launch it
and which outputs should be audited.

For nf-core entries, `workflow_provider` is `nf-core`, `software_stack` is
`nf-core`, the registry entry stores the external `source` and `release`, and
CBIcall runs `nextflow run <source> -r <release>`.

The current lightweight smoke-test example is nf-core/demo:

```yaml
workflows:
  nextflow:
    base_dir: "workflows/nextflow"
    software_stacks:
      nf-core:
        pipelines:
          demo:
            single:
              default_registry_version: "v1"
              registry_versions:
                v1:
                  provider: "nf-core"
                  source: "nf-core/demo"
                  release: "1.1.0"
                  default_outdir: "demo"
          sarek:
            cohort:
              default_registry_version: "v1"
              registry_versions:
                v1:
                  provider: "nf-core"
                  source: "nf-core/sarek"
                  release: "3.8.1"
                  default_outdir: "sarek"
                  canonical_outputs:
                    - name: "haplotypecaller_vcf"
                      type: "vcf"
                      pattern: "variant_calling/haplotypecaller/*/*.haplotypecaller.vcf.gz"
```

CBIcall validates the selected YAML, pins the registered release, writes a
params file in the run directory, and runs `nextflow run <source> -r <release>`. The
workflow keeps its upstream output layout; CBIcall does not force it into
`01_*`, `02_*`, or `03_*` directories.

Use `canonical_outputs` for external workflows that produce final deliverables
in upstream workflow directories. CBIcall resolves these patterns under the declared output directory, records matches in `run-report.json`, and computes normalized
VCF hashes for `compare-runs`.

When adding another nf-core workflow, add a registry entry with `provider`,
`source`, `release`, `default_outdir`, and any canonical final outputs, then add
a compatible resource catalog entry if users should be able to select it through
the `resource` key.
For nf-core workflows where references and containers are managed by the
nf-core profile, use a lightweight `type: nextflow-managed` resource entry.
See [Adding Resources](adding-resources) for the resource catalog contract.

## 2. Register the Pipeline

Edit:

```text
workflows/registry/cbicall-workflow-registry.yaml
```

Add the pipeline under the correct backend and software stack.

The registry hierarchy is:

```text
workflows
  <backend>
    base_dir
    software_stacks
      <software_stack>
        helpers
        profiles
        pipelines
          <pipeline>
            <mode>
              default_registry_version
              registry_versions
                <registry_version>
                  script | provider/source/release
```

```yaml
workflows:
  bash:
    base_dir: "workflows/bash"
    software_stacks:
      gatk-4.6:
        helpers:
          env: "env.sh"
          coverage: "coverage.sh"
          jaccard: "jaccard.sh"
          vcf2sex: "vcf2sex.sh"
          vcf2hash: "vcf2hash.sh"
        profiles:
          cnag-hpc:
            env: "cnag-hpc-env.sh"
        pipelines:
          mypipe:
            single:
              default_registry_version: "v1"
              registry_versions:
                v1:
                  script: "mypipe_single.sh"
```

For Bash, the helper set is required at the software-stack level. If you add a
pipeline under an existing stack such as `gatk-4.6`, keep the existing helpers
and add only the new `pipelines.<pipeline>` block.

If cohort mode is supported:

```yaml
pipelines:
  mypipe:
    single:
      default_registry_version: "v1"
      registry_versions:
        v1:
          script: "mypipe_single.sh"
    cohort:
      default_registry_version: "v1"
      registry_versions:
        v1:
          script: "mypipe_cohort.sh"
```

:::info[Registry paths]
Registry paths are relative to the backend/software-stack directory. For Bash GATK 4.6, `mypipe_single.sh` resolves below `workflows/bash/gatk-4.6/`.
:::

:::tip[Registry versions]
The software-stack level above identifies the execution stack, for example
`gatk-4.6` or `nf-core`. The nested `v1` is the CBIcall registry
version. Normal run YAML files do not need to set it because the registry
`default_registry_version` is used. If a workflow script or external launch
contract changes in a way that may affect outputs, add a new registry version
such as `v2`, point it to the new script or release, and choose whether to move
`default_registry_version` to that version.
:::

## 3. Make the Pipeline Selectable

The registry makes the script discoverable, but Python validation still controls which `pipeline`, `mode`, and `software_stack` combinations are allowed.

If `mypipe` is a new pipeline name, update the validation sets in `src/cbicall/config.py`:

```python
PIPELINE_VALUES = {"wes", "wgs", "mit", "mypipe"}

_ALLOWED_COMBOS = {
    "gatk-4.6": {
        "wes": ["single", "cohort"],
        "wgs": ["single", "cohort"],
        "mypipe": ["single"],
    },
}
```

If you are only adding a missing mode or script for an existing pipeline and software stack, the registry entry may be enough.

## 4. Add a User-Facing YAML Example

Create a minimal parameters file that exercises the new workflow:

```yaml
mode:            single
pipeline:        mypipe
workflow_backend: bash
software_stack:    gatk-4.6
input_dir:       SAMPLE01
genome:          b37
```

Run it:

```bash
bin/cbicall run -p mypipe.yaml -t 4
```

CBIcall should create a run directory similar to:

```text
SAMPLE01/cbicall_bash_gatk-4.6_mypipe_single_b37_<run-id>/
```

## 5. Validate the Addition

Check the pipeline at three levels.

| Level | What to verify |
| --- | --- |
| Registry | The workflow appears under the correct backend/software-stack/pipeline/mode/registry-version path. |
| Resource catalog | Compatible bundle entries point to real registry workflow keys. |
| Files | Referenced scripts exist; Bash scripts are executable. |
| Runtime | CBIcall creates a run directory, writes `log.json`, `cbicall-execution-contract.json`, `run-report.json`, and produces expected outputs. |

Good first checks:

```bash
bin/cbicall validate-registry
bin/cbicall validate-resources
bin/cbicall validate-parameters -p mypipe.yaml
bin/cbicall run -p mypipe.yaml -t 2
```

Then inspect:

```text
<run-dir>/log.json
<run-dir>/cbicall-execution-contract.json
<run-dir>/logs/
```


:::tip[Validate the workflow language too]
`validate-registry` checks that CBIcall can resolve the entrypoint and compatible
resources. It does not replace backend-native workflow checks. Use `bash -n` for
Bash scripts, Snakemake lint or dry-run checks for Snakefiles, Nextflow's own
validation/dry-run tools for Nextflow workflows, and `womtool validate` for WDL
files. When `WOMTOOL_JAR` is set, CBIcall runs the WDL syntax check as part of
`validate-registry`.
:::


## When Python Changes Are Needed

Most workflow additions only need validation changes in `config.py` if the pipeline name is new. Broader Python changes are needed when the execution model changes.

| Change | Likely file |
| --- | --- |
| New YAML key or default | `src/cbicall/config.py` |
| New value in the typed runtime model | `src/cbicall/models.py` |
| New registry resolution behavior | `src/cbicall/workflow_registry.py` |
| New execution backend | `src/cbicall/execution.py` |
| Different command-line launch behavior | `src/cbicall/execution.py` |

:::warning[New backend]
Adding a backend is different from adding a pipeline. Prefer adding a new runner class rather than expanding conditional logic inside an existing runner.
:::

## Contributor Checklist

- [ ] Pick backend, software stack, pipeline name, and mode.
- [ ] Inspect the closest existing workflow in `workflows/{bash|snakemake}/{software-stack}/`.
- [ ] Add workflow entrypoint scripts.
- [ ] Make Bash scripts executable.
- [ ] Register scripts in `workflows/registry/cbicall-workflow-registry.yaml`.
- [ ] Update Python validation if the pipeline name or compatibility matrix changes.
- [ ] Add a minimal YAML example.
- [ ] Run CBIcall and inspect `log.json`, `cbicall-execution-contract.json`, `run-report.json`, logs, and expected outputs.
- [ ] Add or update tests when validation or execution behavior changes.

## Next Steps

- Review the execution model in [Architecture](architecture).
- Check YAML behavior in [Configuration Reference](../help/configuration-reference).
- Document generated files in [Outputs](../help/outputs) if the new pipeline produces user-facing outputs.
