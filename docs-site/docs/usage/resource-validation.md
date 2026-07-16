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

Use [Integration Tests](../validation/integration-tests) when you want to run the shipped
example workflows and compare output VCFs.

For the resource model and JSON examples, see [Adding Resources](../technical-details/adding-resources).
For the current CBIcall-provided bundle, see [Bundle v1](../technical-details/resource-bundle-v1).

## Check the Installation

Run the installation-level check without a parameters YAML:

```bash
cbicall doctor
```

`doctor` validates the packaged workflow registry and resource catalog, then
checks the resource metadata at the `CBICALL_DATA` root. Current installations
are verified from `cbicall-resource-installation.json`, including the catalog
fingerprint, recorded archive-checksum result, and expected top-level layout.
Installations carrying the catalog-pinned `cbicall-resource-id.json` are also
recognized. The command does not rehash the large resource archive. It reports
whether the four supported execution backends are available; missing optional
backends are warnings.

This command verifies the installation as a whole. It does not select a resource
for an analysis or determine whether a particular YAML contract can run.

## Set the Installed Bundle Location

Native workflows read the external bundle location from `CBICALL_DATA`:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
```

The directory should contain `Databases/`, `NGSutils/`, and the installation
metadata created by `cbicall install-resources`. The Python driver applies this
single value to all native backends; users should not edit packaged `env.sh` or
`config.yaml` files inside `site-packages`.

## Validate the Catalog

Validate the default catalog:

```bash
cbicall validate-resources
```

This checks the **catalog shape** and confirms that declared workflow
compatibility keys exist.

<details>
<summary>Validate one resource key or a custom catalog</summary>

| Case | Command |
| --- | --- |
| Native bundle | `cbicall validate-resources --resource cbicall-germline-resources-v1` |
| nf-core/demo | `cbicall validate-resources --resource nf-core-demo-managed-resources-v1` |
| nf-core/Sarek | `cbicall validate-resources --resource nf-core-sarek-managed-resources-v1` |
| Custom catalog | `cbicall validate-resources --catalog /path/to/cbicall-resource-catalog.json --resource my-center-germline-v1` |

</details>

## Validate One Run

Use `validate-parameters` with the parameters YAML that will be launched:

```bash
cbicall validate-parameters -p my-center-wes.yaml
```

This checks the **selected resource** against the resolved workflow and
installed resource metadata when present. Bundle resources can be checked
against `DATADIR` metadata; externally managed resources, such as nf-core/Sarek,
record catalog identity and compatibility but do not use CBIcall bundle
installation checks.

## Runtime Provenance

Runs that reach workflow launch record resource provenance in `log.json` and
`run-report.json`, including failed executions.
Compare repeated runs with:

```bash
cbicall compare-runs run_a/ run_b/ --output compare-report.txt
```

For adding a new resource entry, see [Adding Resources](../technical-details/adding-resources).
