# Reviewer Response Notes: Reproducibility and Framework Validation

These notes summarize technical changes and proposed response language for
reviewer comments about reproducibility, external resources, and framework-level
validation. They are working notes for the manuscript revision, not public user
documentation.

## Reviewer Concern

The reviewers raised three related technical points:

1. Third-party tools and resources such as GATK, BWA-MEM, Samtools, reference genomes, known-sites files, and MToolBox resources are downloaded separately and placed in a designated resource directory.
2. Manual resource preparation can become the main source of workflow divergence, even if the CBIcall wrapper and workflow scripts are versioned.
3. CBIcall's reproducibility claims should be demonstrated at the framework level, rather than inferred only from plausible variant-calling results.

The core issue is not only software provenance. It is resource-level reproducibility: the exact external executables, reference files, interval files, known-sites databases, and directory layout used by the workflows.

The installation utility verifies the downloaded archive after the split archive parts are reassembled. That protects against a corrupted or incomplete download. Because Google Drive can stall during large downloads, the same utility also supports a manual-download route: users can download the listed files in a browser and then rerun the setup with `--skip-download` to assemble, verify, extract, and record the installation.

The remaining reproducibility gap is not transfer integrity. It is semantic installation identity: the bundle must be unpacked, placed under the expected directory layout, and connected to workflow-specific `env.sh` or Snakemake configuration files. CBIcall now addresses this with a resource catalog (`resources/cbicall-resource-catalog.json`) and runtime checks that connect a selected workflow implementation to a declared bundle.

## Current Design

CBIcall separates three layers.

| Layer | Current behavior |
| --- | --- |
| Execution layer | Python wrapper, versioned workflow registry, Bash/Snakemake workflow scripts, and optional containers. |
| External resource layer | Catalogued bundle with third-party tools, reference genomes, known-sites VCFs, interval lists, and MToolBox resources installed under a configured data directory. |
| Runtime provenance | `log.json` and `run-report.json` record CLI arguments, resolved configuration, selected workflow implementation version, bundle identity, paths, and runtime metadata. Workflow logs capture executed commands and tool paths. Nuclear WES/WGS single-sample Bash workflows also write a VCF SHA-256 fingerprint report under `03_stats`. |

The connection between CBIcall and the external resource layer is made through workflow-specific configuration files, especially `env.sh`.

For example, `env.sh` maps the external directory layout to variables used by the workflow:

| Variable | Meaning |
| --- | --- |
| `DATADIR` | Root data/resource directory. |
| `DBDIR` | Reference database directory below `DATADIR`. |
| `NGSUTILS` | Tool installation directory below `DATADIR`. |
| `GATK4_BIN` | GATK 4 executable path. |
| `BWA` | BWA executable path. |
| `SAM` | Samtools executable path. |
| `REF` / `REFGZ` | Reference genome FASTA paths. |
| `dbSNP` | dbSNP VCF path. |
| `MILLS_INDELS` / `KG_INDELS` | Known INDEL resource paths. |
| `INTERVAL_LIST` | WES target interval list. |
| `MTOOLBOXDIR` / `MTOOLBOXDB` | MToolBox executable and database locations. |

This means the workflow is reproducible when the code version, user parameter file, workflow registry entry, pipeline implementation version, selected bundle, resource layout, and resource contents are all preserved. The current implementation records the selected resource key and catalog fingerprint, and checks local bundle metadata when `cbicall-resource-installation.json` or `cbicall-resource-id.json` is present in `DATADIR`.

## Minimal Reproducibility Check

The documentation now exposes the bundled integration test as a minimal framework-level reproducibility check:

```bash
bin/cbicall validate-registry
bin/cbicall validate-resources
bin/cbicall doctor -p examples/input/param.yaml
bin/cbicall test --wes -t 1
```

This sequence demonstrates:

1. The workflow registry conforms to its JSON Schema.
2. The resource catalog is well formed and all compatible workflow keys exist in the workflow registry.
3. The user-facing YAML resolves to a declared workflow, profile, pipeline implementation version, and selected bundle without launching GATK.
4. The normal CLI can execute the bundled WES workflow and reproduce the shipped reference VCF under deterministic comparison rules.
5. The integration test reports a normalized SHA-256 fingerprint after removing VCF headers and sorting variant records.

This is deliberately framed as framework-level integration validation. It does not claim to validate the biological accuracy of GATK, BWA-MEM, or MToolBox.

## Cloud Portability Check

The documentation now includes a Google Cloud Compute Engine recipe using Docker. The recipe runs the same minimal reproducibility check on a fresh cloud VM:

```bash
bin/cbicall validate-registry
bin/cbicall validate-resources
bin/cbicall doctor -p examples/input/param.yaml
bin/cbicall test --wes -t 1
```

The purpose is to provide direct portability evidence for one sample, not to benchmark cloud-scale throughput. The recommended evidence to keep for the rebuttal is:

- CBIcall commit hash.
- Docker image digest.
- Output of `validate-registry`, `validate-resources`, `doctor`, and `test --wes`.
- `run-report.json`.
- Normalized variant-record SHA-256 for the generated VCF.

Suggested response language:

> To directly assess portability, we added and executed a documented Google Cloud Compute Engine recipe using Docker. The run uses the same CBIcall checkout, user YAML, workflow registry, Docker image, and selected bundle as the local test. CBIcall resolves the same workflow implementation and bundle identity, and the generated VCF is checked against the shipped reference output using the deterministic comparison implemented in the integration test.

## Output Fingerprints

The Bash WES/WGS single-sample workflows now call `vcf2hash.sh` after the final QC VCF is written. The report is written to `03_stats/<sample>.vcf.sha256.txt` and includes:

| Field | Meaning |
| --- | --- |
| `RAW_SHA256` | Byte-level file hash. Useful for exact artifact identity, but too strict for cross-run VCF reproducibility. |
| `NORMALIZED_PATTERN` | Header pattern removed before normalized hashing, currently `^#`. |
| `NORMALIZED_SORT` | Sort rule used before hashing, currently `LC_ALL=C`. |
| `NORMALIZED_RECORDS` | Number of non-header VCF records included in the normalized hash. |
| `NORMALIZED_SHA256` | Compact fingerprint of the header-stripped, sorted variant-record stream. |

The normalized hash is intended as a reporting fingerprint, not a replacement for the integration-test diff. It avoids unstable VCF headers and gzip metadata while still detecting changes in the variant records.

## Limitation of Versioning Only Tools

Tool versions alone are not sufficient.

For example, two runs can both use GATK 4.6 and BWA 0.7.18 but still diverge if they use:

- different reference FASTA builds
- different dbSNP releases
- different known-sites VCFs
- different interval lists
- modified local resource files
- different MToolBox database contents

Therefore, a stronger reproducibility model should describe the full resource set, not just executable versions.

## Implemented Resource Catalog

CBIcall now includes a resource catalog at `resources/cbicall-resource-catalog.json`. Its supported resource type today is `bundle`.

The catalog describes supported resource sets and their compatibility with workflow engine, pipeline, mode, GATK version, and CBIcall pipeline implementation version. It also describes how the downloaded archive maps to an installed resource tree.

Example sketch:

```yaml
catalog_version: 1
resource_key: cbicall-gatk46-b37-wes-2026.05
type: bundle
description: GATK 4.6 WES/WGS resources for b37 workflows
workflow_engine: bash
gatk_version: gatk-4.6
genome: b37
pipelines:
  - wes
  - wgs

layout:
  datadir: /cbicall-data
  dbdir: Databases
  ngsutils: NGSutils

archive:
  filename: data.tar.gz
  canonical_name: cbicall-germline-resources-v1.tar.gz
  checksum_file: data.tar.gz.md5
  integrity_check: assemble parts, compute MD5, compare with data.tar.gz.md5
  install_action: safe extract into DATADIR
  manual_download_supported: true

remote_identifier:
  filename: cbicall-resource-id.json
  expected:
    resource_key: "cbicall-germline-resources-v1"
  sha256: optional registry-pinned SHA-256 of the small identifier file

tools:
  gatk:
    version: 4.6.2.0
    path: NGSutils/gatk/gatk-4.6.2.0/gatk
  bwa:
    version: 0.7.18
    path: NGSutils/bwa-0.7.18/bwa
  samtools:
    version: 0.1.19
    path: NGSutils/samtools-0.1.19/samtools

resources:
  reference_fasta:
    label: b37 reference FASTA
    path: Databases/GATK_bundle/b37/references_b37_Homo_sapiens_assembly19.fasta
  dbsnp:
    label: dbSNP b144 GRCh37p13
    path: Databases/dbSNP/human_9606_b144_GRCh37p13/All_20160408.vcf.gz
  mills_indels:
    label: Mills and 1000G gold-standard INDELs
    path: Databases/GATK_bundle/b37/b37_Mills_and_1000G_gold_standard.indels.b37.vcf.gz
  exome_intervals:
    label: Broad b37 human exome intervals
    path: Databases/GATK_bundle/b37/b37_Broad.human.exome.b37.interval_list
```

Future command-line extensions could expose this catalog directly, for example:

```bash
cbicall resources list
cbicall resources install --bundle cbicall-gatk46-b37-wes-2026.05
cbicall resources verify --bundle cbicall-gatk46-b37-wes-2026.05
cbicall resources doctor --param wes_single.yaml
```

In this model, the catalog is not just a list of checksums. It is the bridge between:

- the local CBIcall compatibility rules and resource descriptions
- the small remote resource identifier file
- the downloaded archive
- the unpacked directory layout
- the variables expected by `env.sh`
- the resources used by each workflow
- the CBIcall and workflow implementation versions known to be compatible

## Archive Integrity vs Bundle Versioning

Checksums are used at two levels. The local CBIcall registry can optionally pin the SHA-256 of the small `cbicall-resource-id.json` file, which verifies that the remote bundle declares the expected resource key. The installer then verifies the assembled archive using the distributed archive checksum. This confirms that the `data.tar.gz` assembled from the split parts matches the expected bundle payload.

That is different from versioning the selected resource entry. After the archive is unpacked, CBIcall still needs a machine-readable way to know which bundle was installed and whether the local `env.sh` or Snakemake configuration points to the expected files.

A practical registry should therefore make the installed bundle identifiable and verifiable without requiring a checksum for every file in every resource tree.

Reasons:

- Some third-party bundles contain many files, indexes, caches, and generated side files.
- Public resources do not always come with stable semantic bundle versions.
- Some files may be rebuilt locally, for example FASTA indexes or sequence dictionaries.
- Strict file-by-file checksums can make harmless local regeneration look like incompatibility.

A practical registry could support selective integrity checking.

| Resource type | Suggested checksum policy |
| --- | --- |
| Primary reference FASTA | Strongly recommended |
| Known-sites VCFs | Strongly recommended |
| Interval/BED files | Strongly recommended |
| Tool binaries or release archives | Recommended |
| Index files generated from primary resources | Optional or derived |
| Large third-party directories | Optional directory manifest or release identifier |
| Locally generated cache files | Not required |

The registry could therefore record the archive-level checksum that is already used during download, plus optional selected checksums or release labels where they are most meaningful, without requiring every file to be pinned.

Example:

```yaml
resources:
  archive:
    path: data.tar.gz
    checksum_file: data.tar.gz.md5
    checksum_policy: install-time
  reference_fasta:
    path: Databases/GATK_bundle/b37/references_b37_Homo_sapiens_assembly19.fasta
    sha256: "..."
    checksum_policy: required
  reference_fai:
    path: Databases/GATK_bundle/b37/references_b37_Homo_sapiens_assembly19.fasta.fai
    checksum_policy: derived
  mtoolbox_database:
    path: Databases/mtDNA
    release_label: local-mtoolbox-db-2026.05
    checksum_policy: manifest-optional
```

## Compatibility Matrix

The resource catalog is linked to a compatibility matrix.

This matrix describes which bundles are valid for each workflow family and pipeline implementation version.

| CBIcall version | Workflow key | Bundle |
| --- | --- | --- |
| `1.x` | `bash/wes/single/gatk-4.6/v1` | `cbicall-germline-resources-v1` |
| `1.x` | `bash/wgs/single/gatk-4.6/v1` | `cbicall-germline-resources-v1` |
| `1.x` | `snakemake/wes/single/gatk-4.6/v1` | `cbicall-germline-resources-v1` |
| `1.x` | `snakemake/wgs/single/gatk-4.6/v1` | `cbicall-germline-resources-v1` |
| `1.x` | `bash/mit/single/gatk-3.5/v1` | `cbicall-germline-resources-v1` |

This is preferable to hard-coding one global bundle in the Python wrapper, because different sequencing types and analysis modes may need different resources.

## Possible Manuscript / Rebuttal Language

> We agree that long-term reproducibility depends not only on workflow code but also on the external resource layer. We have therefore extended CBIcall to record the selected bundle, resource catalog fingerprint, resolved workflow implementation version, runtime profile, and workflow paths in the run provenance. At runtime, CBIcall resolves the selected `DATADIR` from the active Bash `env.sh` or Snakemake configuration and verifies local bundle metadata when `cbicall-resource-installation.json` or `cbicall-resource-id.json` is present.

> We also agree that manual resource preparation can be a source of workflow divergence. CBIcall now includes a resource catalog with bundle entries describing the archive, expected installation layout, tool versions, reference file hints, and compatible workflow implementation versions. Backend-specific runtime paths are resolved from the active workflow configuration, for example Bash `env.sh` or Snakemake `config.yaml`. The catalog is used during configuration resolution to reject incompatible workflow/resource combinations before execution.

> To make framework-level reproducibility explicit, we added a minimal reproducibility check to the documentation: `bin/cbicall validate-registry`, `bin/cbicall validate-resources`, `bin/cbicall doctor -p examples/input/param.yaml`, and `bin/cbicall test --wes -t 1`. This validates the workflow registry, validates the resource catalog and its workflow compatibility keys, resolves the user-facing YAML without launching GATK, and then executes the bundled WES example against shipped reference outputs. We present this as framework-level integration validation, distinct from biological validation of the underlying variant-calling algorithms.

> We further added a documented Google Cloud Docker recipe for the same one-sample reproducibility check. This provides a concrete portability test outside the local institutional environment. For output-level reporting, the nuclear WES/WGS single-sample Bash workflows now write a VCF hash report containing both the raw file SHA-256 and a normalized SHA-256 over header-stripped, sorted variant records.

## Implementation Status

Implemented in the current revision:

- Documented the exact expected resource layout.
- Documented the relationship between `DATADIR`, `env.sh`, Snakemake configuration, and workflow variables.
- Added a minimal reproducibility command block to the Quickstart.
- Added explicit pipeline implementation versions in the workflow registry.
- Recorded selected bundle provenance in `log.json` and `run-report.json`.
- Added runtime bundle identity checks using local metadata files and lightweight hash validation.
- Validated bundle compatibility against resolved workflow implementation versions.
- Added `doctor`, `validate-registry`, `validate-resources`, and `test` commands.
- Added a Google Cloud Docker recipe for a one-sample portability check.
- Added `vcf2hash.sh` to the Bash workflow registry and WES/WGS single-sample workflows, writing VCF SHA-256 fingerprint reports under `03_stats`.

Still worth considering:

- Add a resource installation/initialization command that writes or validates the expected configuration.
- Add `cbicall resources verify` for optional runtime checks of selected files.
- Add optional resource fetching or initialization for supported bundles.
- Support multiple catalogs for different sequencing assays or institutional deployments.
