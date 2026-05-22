# Architecture

CBIcall is a **thin orchestration layer** around concrete variant-calling workflows.
Its main responsibilities are:

- Reading a parameters YAML file
- Validating required parameters, compatibility rules, and runtime resources
- Resolving the selected **workflow backend**, **software stack**, **pipeline**, **mode**, and **registry version**
- Preparing a deterministic run directory
- Dispatching the appropriate Bash, Snakemake, Nextflow, Cromwell, or nf-core workflow
- Recording logs, workflow fingerprints, resource provenance, and compact run reports

The actual bioinformatics work, such as alignment, variant calling, mtDNA
analysis, and QC, is implemented in modular workflow files that can be extended
or replaced.

---

## Architecture diagram

![CBIcall architecture diagram](/img/architecture-diagram.png)

_CBIcall resolves a validated parameters YAML into one registered workflow implementation and records the resolved runtime context._

<div className="cbicallNotePanel">
  <p><strong>Mental model:</strong> users describe the run once in YAML; CBIcall resolves that request against the workflow registry and resource catalog, then records exactly what was executed.</p>
</div>

---

## Main components

| Component | Role | Main files or directories |
| --- | --- | --- |
| **CLI and validation layer** | Reads the parameters YAML, applies defaults, checks schema-level and compatibility rules, resolves runtime profiles, and starts the selected command. | `src/cbicall/cli.py`, `src/cbicall/config.py` |
| **Workflow registry** | Developer-facing routing table that maps `workflow_backend`, `software_stack`, `pipeline`, `mode`, and `registry_version` to one implementation. Validate it with `bin/cbicall validate-registry` after editing. | `workflows/registry/cbicall-workflow-registry.yaml`, `src/cbicall/workflow_registry.py` |
| **Resource catalog** | Declares external dependency sets, resource identifiers, compatible workflow keys, and checksum metadata for downloadable CBIcall bundles. | `resources/cbicall-resource-catalog.json`, `src/cbicall/resources.py` |
| **Workflow runners** | Execute the resolved implementation. Bash runs local scripts directly; Snakemake, native Nextflow, and Cromwell run bundled workflow files; external nf-core entries are launched through Nextflow. | `src/cbicall/dnaseq.py` |
| **Workflow implementations** | Contain the analysis logic for native WES, WGS, and mtDNA pipelines, plus registered external nf-core workflows. | `workflows/bash/`, `workflows/snakemake/`, `workflows/nextflow/`, `workflows/cromwell/` |
| **Run audit layer** | Writes `log.json`, `run-report.json`, workflow fingerprints, selected resource identity, output fingerprints, and comparison reports. | `src/cbicall/cli.py`, `bin/cbicall compare-runs` |
| **Contract tests** | Run small examples and validate expected output contracts without keeping full `ref_*` run directories in the repository. | `src/cbicall/integration_tests.py`, `tests/fixtures/integration/` |

---

## Directory structure
```
<sample_or_project_dir>/
  cbicall_<backend>_<software-stack>_<pipeline>_<mode>_<genome>_<run-id>/
    01_bam/
    02_varcall/
    03_stats/
    logs/
    log.json
    run-report.json
```

External nf-core runs use a related layout:

```text
<project_dir>/
  cbicall_nextflow_nf-core_<pipeline>_<mode>_<display-genome>_<run-id>/
    <nf-core-output-dir>/
    work/
    pipeline_info/
    log.json
    run-report.json
```

Typical native usage:

- Intermediate alignment and BAM files are stored under `01_bam/`.
- Variant-calling outputs, gVCFs, VCFs and related files are stored under `02_varcall/`.
- Summary statistics and QC metrics are collected under `03_stats/`.
- Log files for workflow steps are stored under `logs/`.

For nf-core workflows, native nf-core outputs remain under the pipeline output
directory, while CBIcall adds the same audit files at the run-directory root.

---

## Resolution model

The parameters YAML selects a workflow in layers:

```yaml
workflow_backend: bash
software_stack: gatk-4.6
pipeline: wes
mode: single
```

CBIcall then resolves `default_registry_version` from the registry unless an
advanced run pins a specific `registry_version`:

```yaml
registry_version: v1
```

For native CBIcall workflows, `software_stack` identifies the implementation
family, for example `gatk-4.6`. For external nf-core workflows, users select
`workflow_provider: nf-core`; CBIcall then resolves the internal software stack
to `nf-core`, and the external workflow release is declared by the registry
entry.

The selected run is therefore identified by:

```text
workflow_backend + software_stack + pipeline + mode + registry_version
```

Resource compatibility is resolved separately through the resource catalog. This
keeps workflow routing and external dependency identity explicit but independent.

<div className="cbicallNotePanel">
  <p><strong>Run identity:</strong> the workflow registry identifies the code path; the resource catalog identifies the external dependency set. Both are written to run metadata.</p>
</div>

---

## Execution model

CBIcall supports two analysis modes:

- **Single mode**  
  One sample is processed independently.

- **Cohort mode**  
  Joint or cohort-level analysis, depending on the selected pipeline.

The workflow backend is selected in the YAML:

- `workflow_backend: bash`
- `workflow_backend: snakemake`
- `workflow_backend: nextflow`
- `workflow_backend: cromwell`

::::tip[Backend choice]
Bash workflows are direct and transparent. Snakemake, native Nextflow, and
Cromwell workflows provide backend-managed orchestration for bundled CBIcall workflows.
External nf-core workflows are also launched through the Nextflow backend, but
their test data, containers, and references are managed by nf-core/Nextflow
rather than by the CBIcall resource bundle.
::::

---

## Runtime profiles

Native CBIcall workflows can expose more than one environment mapping in the
workflow registry. The default profile is `local`; an operator can select another
profile from the CLI:

```bash
bin/cbicall run -p parameters.yaml --runtime-profile cnag-hpc
```

Profiles are intentionally not part of the user YAML. They are deployment
choices, not biological analysis parameters. The resolved profile and selected
environment file are recorded in `log.json` for real runs.

---

## Run outputs and audit files

Every run records two complementary metadata files:

| File | Purpose |
| --- | --- |
| **`log.json`** | Detailed run log with resolved configuration, selected workflow, runtime profile, resource bundle identity, and command context. |
| **`run-report.json`** | Compact audit report for comparing runs, including workflow fingerprints, selected versions, output inventory hashes, canonical output fingerprints, status, and elapsed time. |

Use `bin/cbicall compare-runs` to compare two `run-report.json` files and, by
default, generate both terminal and HTML summaries.

---

## Supported workflows

### Native CBIcall workflows

| Pipeline | Mode | Genome | Software stack | Bash | Snakemake | Nextflow | Cromwell |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **WES** | `single` | `b37` | `gatk-3.5`, `gatk-4.6` | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | `gatk-4.6` | `gatk-4.6` | `gatk-4.6` |
| **WES** | `cohort` | `b37` | `gatk-3.5`, `gatk-4.6` | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | `gatk-4.6` | `gatk-4.6` | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> |
| **WGS** | `single` | `b37`, `hg38` | `gatk-4.6` | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> |
| **WGS** | `cohort` | `b37`, `hg38` | `gatk-4.6` | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> |
| **mtDNA** | `single` | `rsrs` | `gatk-3.5` | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> |
| **mtDNA** | `cohort` | `rsrs` | `gatk-3.5` | <span className="cbicallTestBadge cbicallTestBadgeYes">V</span> | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> |

::::warning[mtDNA platform limitation]
The bundled mtDNA workflow is not supported on ARM / aarch64 because of legacy
third-party dependencies.
::::

### Registered nf-core workflows

| Pipeline | Mode | Registry source | Release | Resource model |
| --- | --- | --- | --- | --- |
| `demo` | `single` | `nf-core/demo` | `1.1.0` | Nextflow/nf-core managed |
| `sarek` | `cohort` | `nf-core/sarek` | `3.8.1` | Nextflow/nf-core managed |

---

## Validation commands

| Command | Checks |
| --- | --- |
| **`bin/cbicall validate-registry`** | Workflow registry structure, schema, referenced files, and compatible resource keys. |
| **`bin/cbicall validate-parameters -p parameters.yaml`** | One concrete run setup, including parameter compatibility, selected workflow, runtime profile, and installed resource identity when applicable. |
| **`bin/cbicall test --wes-bash`** and related flags | Contract integration examples that run selected workflows and check expected output fingerprints. |

---

## Extensibility

New native pipelines can be added without modifying the core execution contract:

- Register the workflow under `workflows/registry/cbicall-workflow-registry.yaml`.
- Place implementation files under the backend and software-stack directory.
- Add a `registry_versions` entry so changes can be pinned and audited.
- Add or update resource catalog compatibility when the workflow needs a CBIcall bundle.
- Add an integration-test fixture when a small runnable example exists.

External workflows can also be registered through the Nextflow backend when they
can be launched through a stable source and release, such as `nf-core/sarek`.

See:

[Adding a pipeline](adding-a-pipeline)
