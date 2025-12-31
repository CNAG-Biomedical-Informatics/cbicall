# Adding a new pipeline

CBIcall is designed to be extensible: you can add new analysis pipelines without changing the core framework *as long as the workflow wiring is declared in the workflow registry*.

In CBIcall, **the source of truth for available workflows is**:

- `workflows/config/cbicall.workflows.yaml` (registry; *what exists*)
- `workflows/schema/workflows.schema.json` (schema; *what a valid registry looks like*)

The Python code validates user parameters (enums, defaults, compatibility rules) and then resolves the correct workflow script **from the YAML registry**.

---

## 1. How CBIcall discovers workflows

At runtime, CBIcall:

1. Reads the user parameter YAML (e.g. `params.yaml`)
2. Applies defaults and validates semantic rules (pipeline/mode/genome/engine compatibility)
3. Loads and schema-validates the workflow registry:
   - `workflows/config/cbicall.workflows.yaml`
   - `workflows/schema/workflows.schema.json`
4. Resolves workflow scripts from the registry and checks:
   - referenced files exist
   - Bash scripts are executable (`+x`)
5. Dispatches the workflow through the selected engine (Bash or Snakemake)

**Key point:** adding a pipeline is mostly **YAML + scripts**, not Python.

---

## 2. Create workflow scripts

Create one or more workflow entrypoints under the engine directory and GATK version:

### Bash

```
workflows/bash/gatk-4.6/
  mypipe_single.sh
  mypipe_cohort.sh
  parameters.sh
```

### Snakemake

```
workflows/snakemake/gatk-4.6/
  mypipe_single.smk
  mypipe_cohort.smk
  config.yaml
```

### Naming conventions

CBIcall follows the pattern:

- `{pipeline}_{mode}.{sh|smk}`
- `pipeline`: e.g. `wes`, `wgs`, `mit`, `mypipe`
- `mode`: `single` or `cohort`

---

## 3. Register the pipeline in the workflow registry

Edit `workflows/config/cbicall.workflows.yaml` and add your pipeline under the appropriate engine + version.

### Minimal example (Bash, single mode)

```yaml
workflows:
  bash:
    base_dir: "workflows/bash"
    versions:
      gatk-4.6:
        common:
          parameters: "parameters.sh"
          coverage: "coverage.sh"
          jaccard: "jaccard.sh"
          vcf2sex: "vcf2sex.sh"
        pipelines:
          mypipe:
            single: "mypipe_single.sh"
```

### Adding cohort mode

```yaml
        pipelines:
          mypipe:
            single: "mypipe_single.sh"
            cohort: "mypipe_cohort.sh"
```

### Snakemake example

```yaml
workflows:
  bash:
    base_dir: "workflows/bash"
    versions:
      gatk-4.6:
        common:
          parameters: "parameters.sh"
          coverage: "coverage.sh"
          jaccard: "jaccard.sh"
          vcf2sex: "vcf2sex.sh"
        pipelines:
          wes:
            single: "wes_single.sh"

  snakemake:
    base_dir: "workflows/snakemake"
    versions:
      gatk-4.6:
        common:
          config: "config.yaml"
        pipelines:
          mypipe:
            single: "mypipe_single.smk"
```

> The schema currently requires `workflows.bash` to exist, even if you primarily use Snakemake. Keep a minimal bash block if needed.

---

## 4. Expose user parameters via the parameter YAML

Users select the workflow by providing:

```yaml
pipeline: mypipe
mode: single
workflow_engine: bash
gatk_version: gatk-4.6
projectdir: cbicall
```

Then add pipeline-specific fields as needed (reference, targets, etc.). The recommended approach is:

- keep “policy” options in YAML (things the user should control)
- keep script-specific internal wiring inside the workflow scripts

---

## 5. Single vs cohort mode

If you support both modes, register both in the YAML and provide both entrypoints.

If cohort mode does not make sense, **omit it** from the registry. CBIcall will raise a clear error if a user requests an unavailable mode.

---

## 6. Validation and guardrails

CBIcall will fail fast when:

- `pipeline`, `mode`, `workflow_engine`, `gatk_version` are invalid
- the pipeline/mode combination is not allowed for the selected GATK version
- incompatible combinations are selected (example: `snakemake` with `gatk-3.5` if restricted)
- a registry entry points to a missing file
- a Bash workflow script is not executable (`+x`)

This keeps “what is allowed” stable and “what is wired” configurable.

---

## 7. Recommended: Add a minimal example dataset

Create a tiny, deterministic example to support CI and user onboarding:

```
examples/mypipe_test/
  params.yaml
  input/
  expected/
```

Goals:

- fast execution (ideally < 1 minute)
- deterministic outputs
- validates one or two key artifacts

---

## 8. Document the pipeline page

Add a new page:

```
docs/pipelines/mypipe.md
```

Include:

- requirements (tools, references)
- supported modes
- required YAML parameters
- example command(s)
- output layout

Then add it to `mkdocs.yml`:

```yaml
nav:
  - Pipelines:
      - Adding a pipeline: pipelines/adding-a-pipeline.md
      - MyPipe: pipelines/mypipe.md
```

---

## 9. Optional: Container support

If the pipeline requires extra software, pin versions for reproducibility:

- extend your Dockerfile, or
- publish a new container image

Make sure the workflow scripts remain runnable in the container environment.

---

## Quick checklist

- [ ] Add workflow script(s) under `workflows/{bash|snakemake}/{gatk-version}/`
- [ ] Register pipeline in `workflows/config/cbicall.workflows.yaml`
- [ ] Ensure registry passes `workflows/schema/workflows.schema.json`
- [ ] Make Bash scripts executable (`chmod +x`)
- [ ] Add example `params.yaml`
- [ ] Add docs page and link it in `mkdocs.yml`
