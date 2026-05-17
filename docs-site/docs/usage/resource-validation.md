# Resource Validation

Resource validation checks whether the declared external dependency layer is compatible with the selected workflow. It covers catalog structure, compatible workflow keys, and runtime identity checks when local metadata is present.

## Validate the Catalog

Validate the default catalog:

```bash
bin/cbicall validate-resources
```

Validate one resource key:

```bash
bin/cbicall validate-resources --bundle cbicall-germline-resources-v1
```

Validate a custom catalog:

```bash
bin/cbicall validate-resources \
  --catalog /path/to/cbicall-resource-catalog.json \
  --bundle my-center-germline-v1
```

This confirms that the catalog entry is well formed and that its `compatible_workflows` keys exist in `workflows/registry/workflows.yaml`.

## Validate One Run

Catalog validation is not enough to prove a concrete run is wired correctly. Use `doctor` with the YAML that will be launched:

```bash
bin/cbicall doctor -p my-center-wes.yaml
```

`doctor` resolves the workflow implementation, checks that the selected `resource` is compatible with it, resolves the backend resource location, and verifies installed metadata when `cbicall-resource-id.json` or `cbicall-resource-installation.json` is present in `DATADIR`.

## Runtime Provenance

Completed runs record resource provenance in `log.json` and `run-report.json`. The compact run report can then be compared across runs:

```bash
bin/cbicall compare-runs run_a/ run_b/ --output compare-report.txt --html compare-report.html
```

For adding a new resource entry, see [Adding Resources](../technical-details/adding-resources).
