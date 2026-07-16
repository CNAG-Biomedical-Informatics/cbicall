import WorkflowCompatibilityMatrix from '@site/src/components/WorkflowCompatibilityMatrix.mdx';

# Architecture

CBIcall is a **configuration-driven execution framework** for reproducible
variant-calling pipelines. It does not replace Bash, Snakemake, Nextflow, or
Cromwell. Instead, it validates one user request, resolves it against approved
workflow-registry and resource definitions, launches the selected backend, and records the
evidence needed to audit the run afterwards.

The architecture is organized around a three-layer execution contract:

| Layer | Purpose |
| --- | --- |
| **User parameters YAML** | Describes the analysis intent: pipeline, mode, genome, input data, workflow backend, workflow provider, software stack, and run-specific parameters. |
| **Workflow registry YAML** | Maps the requested analysis to an approved workflow implementation, entrypoint, backend, provider, version, expected files, and canonical outputs. |
| **Resource catalog JSON** | Defines external resources: reference bundles, workflow compatibility, resource identity, availability, and integrity metadata. |

After these layers are validated and resolved, CBIcall creates a deterministic
run directory, dispatches the selected pipeline implementation, and writes structured audit
artifacts such as `log.json`, `cbicall-execution-contract.json`,
`run-report.json`, optional HTML reports, and comparison reports.

The actual bioinformatics work, such as alignment, variant calling, mtDNA
analysis, and QC, remains implemented in backend workflow files or upstream external
pipelines.

---

## Architecture diagram

![CBIcall execution contract diagram](/img/cbicall-execution-landscape.svg)

_CBIcall validates the user parameters YAML, resolves it against the workflow registry and resource catalog, launches the selected backend, and records audit outputs._

<div className="cbicallNotePanel">
  <p><strong>Mental model:</strong> the user parameters file describes what should be run; the workflow registry decides how it is run; the resource catalog defines which external dependencies are valid.</p>
</div>

---

## Main components

| Component | Role | Main files or directories |
| --- | --- | --- |
| **Execution controller** | Reads the user parameters YAML, applies defaults, validates schema-level and compatibility rules, resolves runtime profiles, and starts the selected backend command. | `src/cbicall/cli.py`, `src/cbicall/config.py` |
| **Validation schemas** | Define the structure and required fields for workflow registry entries and resource catalog entries; user parameters are validated by the configuration layer and compatibility checks. | `workflows/schema/`, `resources/cbicall-resource-catalog.schema.json`, `src/cbicall/config.py` |
| **Workflow registry** | Developer-facing routing table that maps `workflow_provider`, `workflow_backend`, `software_stack`, `pipeline`, `mode`, and `registry_version` to one approved implementation. | `workflows/registry/cbicall-workflow-registry.yaml`, `src/cbicall/workflow_registry.py` |
| **Resource catalog** | Declares external dependency sets, resource identifiers, compatible workflow keys, availability rules, and checksum metadata for downloadable CBIcall bundles. | `resources/cbicall-resource-catalog.json`, `src/cbicall/resources.py` |
| **Workflow runners** | Execute the resolved implementation. Bash runs local scripts directly; Snakemake, native Nextflow, and Cromwell run bundled workflow files; external nf-core entries are launched through Nextflow. | `src/cbicall/execution.py` |
| **Workflow implementations** | Contain native CBIcall workflows and registered external workflow entries. Native workflows follow the CBIcall output contract; external workflows keep their upstream output layout. | `workflows/bash/`, `workflows/snakemake/`, `workflows/nextflow/`, `workflows/cromwell/` |
| **Run audit layer** | Writes `log.json` and `run-report.json`; records runtime versions, workflow and resource identity, output inventories, execution contracts, and VCF fingerprints. | `src/cbicall/run_audit.py`, `src/cbicall/runtime_info.py` |
| **Report layer** | Loads and refreshes completed audits, renders text/HTML/MultiQC summaries, and compares two or more runs. | `src/cbicall/report_commands.py`, `src/cbicall/comparison_commands.py`, `src/cbicall/html_reports.py`, `src/cbicall/multiqc.py` |
| **Contract tests** | Run small examples and validate expected output contracts without keeping full `ref_*` run directories in the repository. | `src/cbicall/integration_tests.py`, `tests/fixtures/integration/` |

:::note[Source and installed layouts]
Source checkouts keep workflows and catalogs in the repository directories
shown above. PyPI wheels package the same files under CBIcall's private runtime
directory. `src/cbicall/paths.py` resolves either layout, so user YAML and
workflow identities do not change with the installation method.

`CBICALL_ROOT` is an advanced developer/test override for selecting another
complete runtime tree. Normal installations do not need it.
:::

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
    cbicall-execution-contract.json
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
    cbicall-execution-contract.json
    run-report.json
```

Typical native usage:

- Intermediate alignment and BAM files are stored under `01_bam/`.
- Variant-calling outputs, gVCFs, VCFs and related files are stored under `02_varcall/`.
- Summary statistics and QC metrics are collected under `03_stats/`.
- Log files for workflow steps are stored under `logs/`.

For external nf-core workflows, upstream nf-core outputs remain under the pipeline
output directory, while CBIcall adds audit files at the run-directory root.

---

## Native and External Implementations

CBIcall uses **native** to describe the output contract, not the workflow
language. A pipeline implementation is CBIcall-native when it writes the expected CBIcall run
layout inside the generated `cbicall_*` directory. An implementation is external when
CBIcall launches and audits it but the upstream pipeline keeps its own output
layout. External provider entries can still participate in run comparison when the
workflow registry declares canonical final outputs.

## Resolution model

The parameters YAML selects an analysis in layers:

```yaml
workflow_provider: cbicall
workflow_backend: bash
software_stack: gatk-4.6
pipeline: wes
mode: single
genome: b37
```

For native pipeline implementations, `workflow_provider: cbicall` is the default. For external
nf-core provider entries, users select `workflow_provider: nf-core` and
`workflow_backend: nextflow`; the external release is declared by the registry
entry.

CBIcall resolves `default_registry_version` from the registry unless an advanced
run pins a specific `registry_version`. The resolved workflow identity therefore
combines the provider, backend, software stack, pipeline, mode, genome, registry
version, and workflow fingerprints. Resource compatibility is resolved
separately through the resource catalog, so workflow routing and external
dependency identity remain explicit but independent.

<div className="cbicallNotePanel">
  <p><strong>Run identity:</strong> the workflow registry identifies the code path and provider; the resource catalog identifies the external dependency set. Both are written to run metadata and comparison reports.</p>
</div>

---

## Execution model

CBIcall supports two analysis modes:

| Mode | Meaning |
| --- | --- |
| `single` | One sample is processed independently. |
| `cohort` | Joint or cohort-level analysis, depending on the selected pipeline. |

Workflow backends execute workflow definitions:

| Workflow backend | Role |
| --- | --- |
| `bash` | Direct native shell workflow execution. |
| `snakemake` | Native Snakemake workflow execution. |
| `nextflow` | Native Nextflow workflow execution or external nf-core dispatch. |
| `cromwell` | Native WDL/Cromwell workflow execution. |

Workflow providers identify where the implementation comes from:

| Workflow provider | Meaning |
| --- | --- |
| `cbicall` | Pipeline implementation is maintained as part of CBIcall and follows the CBIcall output contract. |
| `nf-core` | Pipeline implementation comes from nf-core and keeps its upstream output layout. |

::::tip[Backend choice]
Bash workflows are direct and transparent. Snakemake, Nextflow, and Cromwell
can provide backend-managed orchestration for CBIcall-native implementations when
those implementations produce the CBIcall output contract.
External nf-core provider entries are also launched through the Nextflow backend, but
their test data, containers, and references are managed by nf-core/Nextflow
rather than by the CBIcall resource bundle.
::::

---

## Runtime profiles

Native CBIcall workflows can expose more than one environment mapping in the
workflow registry. The default profile is `local`; an operator can select another
profile from the CLI:

```bash
cbicall run -p parameters.yaml --runtime-profile cnag-hpc
```

Profiles are intentionally not part of the user YAML. They are deployment
choices, not biological analysis parameters. The resolved profile and selected
environment file are recorded in `log.json` for real runs.

---

## Run outputs and audit files

Each run directory created for workflow launch records three complementary
metadata files. Validation failures detected before launch do not create a run
directory; execution failures after launch retain these audit artifacts.

| File | Purpose |
| --- | --- |
| **`log.json`** | Detailed run log with resolved configuration, selected implementation, runtime profile, resource bundle identity, and command context. |
| **`cbicall-execution-contract.json`** | Backend-ready execution plan created after validation and registry/resource resolution, including the command, generated backend files, and controlled environment overrides. |
| **`run-report.json`** | Compact audit report for comparing runs, including execution-contract, workflow, resource, software, output inventory, and canonical output fingerprints, status, and elapsed time. |

Use `cbicall compare-runs` to compare two `run-report.json` files and, by
default, generate both terminal and HTML summaries.

---

## Supported Pipelines and Backends

<WorkflowCompatibilityMatrix />

---

## Validation commands

| Command | Checks |
| --- | --- |
| **`cbicall validate-registry`** | Workflow registry structure, schema, referenced files, and compatible resource keys. |
| **`cbicall validate-parameters -p parameters.yaml`** | One concrete run setup, including parameter compatibility, selected implementation, runtime profile, and installed resource identity when applicable. |
| **`cbicall test --wes-bash`**, `--backend-equivalence`, and related flags | Contract integration examples that run selected pipeline implementations; backend-equivalence mode compares normalized WES VCF output across native backends. |

---

## Extensibility

New native pipelines can usually be added without modifying the core execution layer:

- Register the workflow under `workflows/registry/cbicall-workflow-registry.yaml`.
- Place implementation files under the backend and software-stack directory.
- Add a `registry_versions` entry so changes can be pinned and audited.
- Add or update resource catalog compatibility when the workflow needs a CBIcall bundle.
- Add an integration-test fixture when a small runnable example exists.

External workflows can also be registered through the Nextflow backend when they
can be launched through a stable source and release, such as `nf-core/sarek`.

See:

[Adding a pipeline](adding-a-pipeline)
