# Adding Resources

CBIcall resources are the external dependency layer used by workflows:
third-party tools, reference genomes, known-sites files, interval lists, and
auxiliary databases. The resource catalog can describe different resource
types, for example `bundle` today and container-oriented entries in the future.

Most resource additions today therefore mean adding a bundle. That needs two
pieces:

1. A bundle entry in `resources/cbicall-resource-catalog.json`
2. Backend-specific runtime configuration that points workflows to the installed files

:::info[Scope]
The CBIcall bundle downloader is intentionally scoped to CBIcall-provided
bundle entries. Custom local, HPC module, or institutional resource layouts
can still be declared as bundle entries, but they do not need to use
`scripts/download_cbicall_bundle.py`.
:::

## Resource Model

The current resource catalog structure is:

```text
resources/cbicall-resource-catalog.json
  resources
    <resource-key>
      type: bundle
```

A resource key is the selectable value users put in `resource`. The shipped
CBIcall resource uses `type: bundle`; the schema keeps `type` open so additional
resource models can be introduced without renaming the catalog.

:::note[Resource type templates]
Non-bundle types can be used as catalog templates for future resource models.
For example, `type: docker` may describe an image name, tag, and digest for
provenance or future validation. Current workflows still resolve executable and
reference paths through their backend configuration, such as Bash `env.sh` or
Snakemake `config.yaml`.
:::

The catalog is validated against:

```text
resources/cbicall-resource-catalog.schema.json
```

## What You Are Adding

Before editing files, define the bundle entry.

| Decision | Example | Why it matters |
| --- | --- | --- |
| Resource key | `my-center-germline-v1` | The value users set as `resource`. |
| Compatible workflows | `bash/wes/single/gatk-4.6/v1` | Prevents incompatible workflow/resource combinations. |
| Runtime location | `DATADIR` or Snakemake `datadir` | Tells the selected backend where files are installed. |
| Integrity metadata | archive checksum or resource identifier | Lets CBIcall verify downloaded or installed bundle identity. |

For a custom bundle, users select it in the run YAML:

```yaml
resource: my-center-germline-v1
```

## 1. Add the Bundle Entry

Edit:

```text
resources/cbicall-resource-catalog.json
```

Add the resource under `resources` and declare its type:

```json
{
  "schema_version": 1,
  "resources": {
    "my-center-germline-v1": {
      "type": "bundle",
      "description": "Institutional germline bundle for b37 and hg38 workflows.",
      "compatible_workflows": [
        "bash/wes/single/gatk-4.6/v1",
        "bash/wes/cohort/gatk-4.6/v1"
      ],
      "layout": {
        "expected_top_level": [
          "Databases",
          "NGSutils"
        ]
      }
    }
  }
}
```

`compatible_workflows` uses:

```text
engine/pipeline/mode/gatk_version/pipeline_version
```

The workflow keys must exist in `workflows/registry/workflows.yaml`.

The essential contract for a bundle entry is intentionally small: the resource
key, `type`, and `compatible_workflows`. Other fields are optional metadata used
for human description, CBIcall-provided downloads, or lightweight runtime
identity checks.

:::tip[Keep the catalog backend-agnostic]
The catalog should describe resource identity, compatibility, broad layout, and
optional integrity metadata. Do not put Bash-only environment variable bindings
in the catalog. Backend-specific paths belong in Bash `env.sh` or Snakemake
`config.yaml`.
:::

## 2. Add Integrity Metadata When Available

For CBIcall-provided downloadable bundles, include the download and checksum
metadata needed by `scripts/download_cbicall_bundle.py`.

```json
{
  "remote_identifier": {
    "filename": "cbicall-resource-id.json",
    "sha256": "64-character-sha256-hex",
    "expected": {
      "resource_key": "my-center-germline-v1"
    }
  },
  "archive": {
    "source_name": "data.tar.gz",
    "canonical_name": "my-center-germline-v1.tar.gz",
    "checksum_file": "data.tar.gz.md5",
    "checksum_algorithm": "md5"
  }
}
```

`archive.canonical_name` is only needed by the downloader/installer. It gives
the verified archive a stable local name after `data.tar.gz` has been assembled.
CBIcall does not need the tarball name when launching workflows.

For local or HPC-managed bundles, a full archive description may not exist.
In that case, keep the catalog entry focused on identity and compatibility, and
place an identifier file beside the installed resources when possible:

```json
{"resource_key": "my-center-germline-v1"}
```

Save it as:

```text
DATADIR/cbicall-resource-id.json
```

When the catalog pins `remote_identifier.sha256`, CBIcall checks that file
against the expected SHA-256 at runtime.

## 3. Point Workflows to the Installed Files

The catalog does not launch tools directly. The selected backend resolves
concrete paths.

For Bash workflows, update the relevant profile `env.sh`:

```bash
DATADIR=/path/to/my-center-resources
DBDIR=$DATADIR/Databases
NGSUTILS=$DATADIR/NGSutils
```

For Snakemake workflows, update the relevant `config.yaml`:

```yaml
datadir: /path/to/my-center-resources
```

The run provenance records the selected resource key and fingerprint. Workflow
logs record the concrete executable and reference paths used during execution.

## 4. Validate the Bundle Entry

Validate the whole catalog:

```bash
bin/cbicall validate-resources
```

Validate one non-default bundle entry:

```bash
bin/cbicall validate-resources --bundle my-center-germline-v1
```

Validate a bundle from a custom catalog file:

```bash
bin/cbicall validate-resources \
  --catalog /path/to/cbicall-resource-catalog.json \
  --bundle my-center-germline-v1
```

This checks that the catalog entry is well formed and that its
`compatible_workflows` keys exist in the workflow registry.

## 5. Validate One Real Run

Resource catalog validation does not know which profile or runtime directory a user will
select. Use `validate-param` with a concrete parameters file to validate the selected
workflow and bundle together:

```yaml
pipeline: wes
mode: single
workflow_engine: bash
gatk_version: gatk-4.6
resource: my-center-germline-v1
genome: b37
input_dir: SAMPLE01
```

Then run:

```bash
bin/cbicall validate-param -p my-center-wes.yaml
```

`validate-param` resolves the workflow implementation, checks that the selected
`resource` is compatible with it, resolves the backend resource location,
and verifies installed bundle metadata when `cbicall-resource-id.json` or
`cbicall-resource-installation.json` is present in `DATADIR`.

## Contributor Checklist

- [ ] Choose a resource key that users can place in `resource`.
- [ ] Add the bundle entry to `resources/cbicall-resource-catalog.json`.
- [ ] Declare compatible workflow implementation keys.
- [ ] Add archive and identifier checksums for CBIcall-provided downloads.
- [ ] Keep backend-specific paths in `env.sh` or Snakemake `config.yaml`.
- [ ] Run `bin/cbicall validate-resources --bundle <resource-key>`.
- [ ] Run `bin/cbicall validate-param -p <parameters.yaml>` with that resource selected.
