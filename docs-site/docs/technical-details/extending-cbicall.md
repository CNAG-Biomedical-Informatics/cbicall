# Developer Overview

Use this section when you are maintaining CBIcall itself: adding workflows,
editing the workflow registry, or declaring resource entries. The **workflow
registry** maps YAML workflow choices to concrete Bash, Snakemake, or Nextflow entrypoints;
the **resource catalog** lists external resource entries, their versions, and
their compatibility metadata. Normal users usually only need the Run, Pipelines, and
Reproducibility sections.

CBIcall can be extended in two main ways:

<div className="cbicallCardGrid">
  <div className="cbicallCard">
    <span className="cbicallCardLabel">Workflows</span>
    <h3>Add a Pipeline</h3>
    <p>Add a Bash, Snakemake, or Nextflow entrypoint and register it in the workflow registry, <code>workflows/registry/cbicall-workflow-registry.yaml</code>. Users can then select it from the parameters YAML.</p>
  </div>

  <div className="cbicallCard">
    <span className="cbicallCardLabel">Resources</span>
    <h3>Add a Resource</h3>
    <p>Add a resource catalog entry that declares its resource version, which workflow keys it supports, and how installed resource identity can be checked.</p>
  </div>
</div>

## What Extension Means

Extending CBIcall does not usually mean changing the Python execution driver.
For most additions, the stable contract is:

| Extension point | Main file | Purpose |
| --- | --- | --- |
| Workflow registry | `workflows/registry/cbicall-workflow-registry.yaml` | Developer-facing file that declares available workflow implementations. |
| Workflow scripts | `workflows/bash/...`, `workflows/snakemake/...`, or `workflows/nextflow/...` | Implement the actual analysis steps. |
| Resource catalog | `resources/cbicall-resource-catalog.json` | Declares external resource versions, compatibility, and identity metadata. |
| Parameters YAML | user-provided `*.yaml` | Selects one registered workflow and resource for a run. |

Python changes are usually needed only when introducing a new execution model,
a new top-level parameter, or a new compatibility rule.

## Where To Start

- To add a workflow implementation, start with [Adding a Pipeline](adding-a-pipeline).
- To add a reference/tool resource set, start with [Adding Resources](adding-resources).
- To understand how the pieces connect, start with [Architecture](architecture).

After editing the workflow registry or resource files, run the relevant checks:

```bash
bin/cbicall validate-registry
bin/cbicall validate-resources
bin/cbicall validate-param -p parameters.yaml
```

`validate-registry` is mainly for pipeline developers. It checks that CBIcall's
workflow registry is structurally valid; it does not run a workflow.
