# Adding a new pipeline

!!! abstract "What you will learn"

    This page explains how to integrate a new workflow into CBIcall by registering it in the workflow registry and adding the corresponding workflow scripts.

CBIcall is designed to be **extensible**. New pipelines are added by providing workflow scripts and registering them in the workflow registry.

The components involved are:

- `workflows/config/cbicall.workflows.yaml` — workflow registry describing available pipelines
- `workflows/schema/workflows.schema.json` — schema used to validate the registry

The execution driver reads the user configuration, resolves the requested workflow from the registry, creates a runtime directory, and executes the workflow script.

---

## 1. How CBIcall runs workflows

At runtime, CBIcall:

1. Reads the user configuration file (e.g. `params.yaml`)
2. Validates parameters and compatibility rules
3. Loads the workflow registry
4. Resolves the workflow script corresponding to the selected pipeline
5. Creates a runtime directory and launches the workflow

Because the workflow registry defines the available pipelines, adding a new workflow typically involves only **a script and a registry entry**.

### Execution layers

The runtime path is split into two responsibilities:

- `src/cbicall/config.py` validates parameters, applies semantic rules, and assembles resolved runtime metadata
- `src/cbicall/workflow_registry.py` loads the workflow registry, resolves the workflow entrypoint, and validates referenced files
- `src/cbicall/dnaseq.py` dispatches execution to an engine-specific runner

The current execution runners are:

- `BashRunner`
- `SnakemakeRunner`

This separation is intentional: adding a new workflow usually means adding scripts plus a registry entry, while adding a new **engine** should be done by introducing a new runner class rather than growing a single dispatcher.

---

## 2. Create workflow scripts

Workflow entrypoints are stored under the engine directory and tool version.

Example structure:

```
workflows/bash/gatk-4.6/
  mypipe_single.sh
  mypipe_cohort.sh
  env.sh
```

or

```
workflows/snakemake/gatk-4.6/
  mypipe_single.smk
  mypipe_cohort.smk
  config.yaml
```

### Engine-specific expectations

For Bash workflows:

- the registry must define the common helper scripts required by the selected version
- the pipeline entrypoint is executed directly as a command
- `GENOME` is exported in the runtime environment

For Snakemake workflows:

- the registry must define a version-specific `config.yaml`
- CBIcall launches `snakemake` with the resolved Snakefile via `-s` and the shared config file via `--configfile`
- `genome` is always passed to Snakemake via `--config`
- `pipeline`, `sample_map`, and `workspace` are added when required by the selected mode/version
- if `workflow_rule` is set, CBIcall uses that rule as the Snakemake target instead of `all`

---

## 3. Naming conventions

Workflow filenames follow the pattern:

```
{pipeline}_{mode}.{sh|smk}
```

Where:

- `pipeline` is the pipeline name (e.g. `wes`, `wgs`, `mit`)
- `mode` is either `single` or `cohort`

Examples:

```
wes_single.sh
wes_cohort.sh
mypipe_single.smk
```

---

## 4. Register the pipeline

Edit the workflow registry:

```
workflows/config/cbicall.workflows.yaml
```

Example entry:

```yaml
workflows:
  bash:
    base_dir: "workflows/bash"
    versions:
      gatk-4.6:
        common:
          env: "env.sh"
        pipelines:
          mypipe:
            single: "mypipe_single.sh"
```

If cohort mode is supported:

```yaml
pipelines:
  mypipe:
    single: "mypipe_single.sh"
    cohort: "mypipe_cohort.sh"
```

---

## 5. Minimal pipeline example

The following script illustrates a minimal workflow.

!!! note "Execution context"

    Workflow scripts are executed **inside the runtime directory** created by CBIcall.  
    Input FASTQ files are typically located in the **parent sample directory**.

Example script:

```
workflows/bash/gatk-4.6/mypipe_single.sh
```

```bash
#!/usr/bin/env bash
set -eu

BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNDIR=$(pwd)

mkdir -p logs

echo "Running example pipeline in $RUNDIR"

for R1 in ../*_R1_*fastq.gz; do
    echo "Found FASTQ: $R1"
done

echo "Pipeline finished"
```

Make it executable:

```
chmod +x workflows/bash/gatk-4.6/mypipe_single.sh
```

---

## 6. Select the pipeline in params.yaml

Example configuration:

```yaml
pipeline: mypipe
mode: single
workflow_engine: bash
gatk_version: gatk-4.6
input_dir: SAMPLE01
```

Run:

```
cbicall -p params.yaml -t 4
```

CBIcall validates the configuration, creates the runtime directory, and executes the workflow script.

If `input_dir` is defined, the runtime directory is created below that input directory with a generated run identifier appended to the configured `project_dir` prefix.

---

## 7. Validation

CBIcall checks that:

- the pipeline is registered in the workflow registry
- the selected execution backend exists
- the selected backend/version/pipeline/mode combination is allowed by the Python validation layer
- referenced workflow scripts are present
- Bash scripts are executable

This ensures that workflows are always executed from **registered and validated configurations**.

---

## 8. Optional Python configuration

Most pipelines only require workflow scripts and registry entries.

If additional parameters or configuration options are needed, the Python changes may span several modules depending on the type of change:

- `src/cbicall/config.py` for parameter defaults, semantic validation, and resolved runtime metadata
- `src/cbicall/workflow_registry.py` for workflow resolution and registry-related validation
- `src/cbicall/models.py` if a new option must become part of the typed internal runtime/config model
- `src/cbicall/dnaseq.py` if the change affects execution behavior or engine-specific command building

For most workflow additions, no Python code changes are required beyond scripts and registry entries.

If you are unsure where a change belongs, inspect the existing implementations in the repository first and follow the closest matching pattern.

---

# Quick checklist

- [ ] Inspect existing pipelines in `workflows/{bash|snakemake}/{gatk-version}/`
- [ ] Check the existing Python implementation paths in `src/cbicall/` when adding new parameters or execution behavior
- [ ] Add workflow script(s)
- [ ] Register the pipeline in `workflows/config/cbicall.workflows.yaml`
- [ ] Validate registry against the JSON schema
- [ ] Make Bash scripts executable (`chmod +x`)
- [ ] Provide example `params.yaml`
