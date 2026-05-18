# Resource Validation

In CBIcall, a **resource** is the external dependency layer used by workflows:
third-party tools, reference genomes, known-sites files, interval lists, and
auxiliary databases. A run selects one resource key in the parameters YAML:

```yaml
resource: cbicall-germline-resources-v1
```

The **resource catalog** is the JSON inventory of those resource entries:
`resources/cbicall-resource-catalog.json`. It records **resource type**,
**resource version**, **workflow compatibility**, and optional identity metadata.

<div className="cbicallNotePanel">
  <p><strong>Use this page when the question is:</strong> does this resource catalog, selected resource key, or installed resource directory match the workflow I am about to run?</p>
</div>

Use [Integration Tests](integration-tests) when you want to run the shipped
example workflows and compare output VCFs.

For the resource model and JSON examples, see [Adding Resources](../technical-details/adding-resources).
For the current CBIcall-provided bundle, see [Bundle v1](../technical-details/resource-bundle-v1).

## Validate the Catalog

Validate the default catalog:

```bash
bin/cbicall validate-resources
```

Validate one resource key:

```bash
bin/cbicall validate-resources --resource cbicall-germline-resources-v1
```

For external nf-core adapters:

```bash
bin/cbicall validate-resources --resource nf-core-demo-managed-resources-v1
bin/cbicall validate-resources --resource nf-core-sarek-managed-resources-v1
```

Validate a custom catalog:

```bash
bin/cbicall validate-resources \
  --catalog /path/to/cbicall-resource-catalog.json \
  --resource my-center-germline-v1
```

This checks the **catalog shape** and confirms that declared workflow
compatibility keys exist.

## Validate One Run

Use `validate-param` with the parameters YAML that will be launched:

```bash
bin/cbicall validate-param -p my-center-wes.yaml
```

This checks the **selected resource** against the resolved workflow and
installed resource metadata when present. Bundle resources can be checked
against `DATADIR` metadata; externally managed resources, such as nf-core/Sarek,
record catalog identity and compatibility but do not use CBIcall bundle
installation checks.

## Runtime Provenance

Completed runs record resource provenance in `log.json` and `run-report.json`.
Compare repeated runs with:

```bash
bin/cbicall compare-runs run_a/ run_b/ --output compare-report.txt --html compare-report.html
```

For adding a new resource entry, see [Adding Resources](../technical-details/adding-resources).
