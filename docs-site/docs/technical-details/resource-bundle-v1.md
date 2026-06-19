# Bundle v1

CBIcall production workflows require a CBIcall-provided bundle containing third-party executables, reference genomes, known-sites files, interval lists, and auxiliary databases.

The current bundle is selected in the run YAML with:

```yaml
resource: "cbicall-germline-resources-v1"
```

CBIcall resolves this key against the resource catalog, the JSON inventory of
resource entries and their compatibility metadata:

```text
resources/cbicall-resource-catalog.json
```

The resource catalog can contain different resource types. This page documents
the current `bundle` type and the CBIcall-provided bundle used by the packaged
workflows.

The `scripts/download_cbicall_bundle.py` utility is intentionally scoped to CBIcall-provided bundle entries in the resource catalog. It is not a general-purpose installer for arbitrary local, HPC module, or third-party layouts. The downloader uses the local catalog when it is present. If the script is used standalone and the local catalog is absent, it fetches the canonical catalog URL before selecting the bundle entry.

## Bundle Identity

| Field | Value |
| --- | --- |
| Resource key | `cbicall-germline-resources-v1` |
| Version | `v1` |

:::info[Why this is explicit]
The resource key is the identifier users select with `resource`. `log.json`
records this key, the resource version, and a catalog fingerprint, so two runs
can be checked for the same declared external dependency set.
:::

## Compatible Native Workflows

This bundle is compatible with the packaged CBIcall-native workflows below.

| Workflow backend | WES | WGS | mtDNA |
| --- | --- | --- | --- |
| `bash` | `single`, `cohort` with `gatk-3.5` or `gatk-4.6` | `single`, `cohort` with `gatk-4.6` | `single`, `cohort` with `gatk-3.5` |
| `snakemake` | `single`, `cohort` with `gatk-4.6` | `single`, `cohort` with `gatk-4.6` | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> |
| `nextflow` | `single`, `cohort` with `gatk-4.6` | `single`, `cohort` with `gatk-4.6` | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> |
| `cromwell` | `single`, `cohort` with `gatk-4.6` | `single`, `cohort` with `gatk-4.6` | <span className="cbicallTestBadge cbicallTestBadgeNo">X</span> |

CBIcall checks the exact workflow key internally, for example
`nextflow/wgs/cohort/gatk-4.6/v1`.

## Downloaded Files

The production bundle is distributed as a small identifier JSON, split archive parts, and a checksum file.

| File | Purpose |
| --- | --- |
| `cbicall-resource-id.json` | Declares the resource key and is pinned by SHA-256 in the catalog. |
| `data.tar.gz.md5` | MD5 checksum file. The current bundle records the split archive parts. |
| `data.tar.gz.part-00` | Split archive part. |
| `data.tar.gz.part-01` | Split archive part. |
| `data.tar.gz.part-02` | Split archive part. |
| `data.tar.gz.part-03` | Split archive part. |
| `data.tar.gz.part-04` | Split archive part. |
| `data.tar.gz.part-05` | Split archive part. |

The setup utility verifies the files covered by `data.tar.gz.md5`, reassembles the parts into an archive, and verifies the downloaded CBIcall-provided bundle payload before extraction.

An optional small remote identifier file can also be used:

```json
{"resource_key": "cbicall-germline-resources-v1"}
```

When available, this file is named `cbicall-resource-id.json`. Its SHA-256 can be pinned in the local catalog to confirm that the remote bundle declares the expected resource key.

## Expected Layout

After extraction, `DATADIR` should contain:

```text
DATADIR/
  Databases/
  NGSutils/
```

The bundle layout uses these conventional top-level names:

| Variable | Meaning |
| --- | --- |
| `DATADIR` | Root of the installed CBIcall-provided bundle. |
| `DBDIR` | `DATADIR/Databases` |
| `NGSUTILS` | `DATADIR/NGSutils` |

Workflow-specific files such as Bash `env.sh` and Snakemake/Nextflow/Cromwell `config.yaml` resolve the installed bundle layout into concrete executable and reference paths.

## Tools

| Tool | Version | Path hint |
| --- | --- | --- |
| GATK 4 | `4.6.2.0` | `NGSutils/gatk/gatk-4.6.2.0/gatk` |
| BWA | `0.7.18` | `NGSutils/bwa-0.7.18/bwa` |
| Samtools | `0.1.19` | `NGSutils/samtools-0.1.19/samtools` |

:::note
Some workflow branches may use architecture-specific executable paths or legacy tool paths. The catalog records the intended bundle identity; the workflow logs and `log.json` record the concrete paths resolved during a run.
:::

## Reference Resources

### b37

| Resource | Path hint |
| --- | --- |
| Reference FASTA | `Databases/GATK_bundle/b37/references_b37_Homo_sapiens_assembly19.fasta` |
| dbSNP | `Databases/dbSNP/human_9606_b144_GRCh37p13/All_20160408.vcf.gz` |
| Mills / 1000G INDELs | `Databases/GATK_bundle/b37/b37_Mills_and_1000G_gold_standard.indels.b37.vcf.gz` |

### hg38

| Resource | Path hint |
| --- | --- |
| Reference FASTA | `Databases/GATK_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta` |

### rsrs / mtDNA

The mtDNA workflows use MToolBox-related mitochondrial files from the bundle.

## Provenance in Runs

The selected resource bundle is stored in `log.json` under `config.resources.bundle`.

Example:

```json
{
  "key": "cbicall-germline-resources-v1",
  "compatible": true,
  "fingerprint": "..."
}
```

Two runs used the same declared external dependency set when their `config.resources.bundle.fingerprint` values match.

## Installation Manifest

The setup utility also writes:

```text
cbicall-resource-installation.json
```

This local manifest records the installed resource key, archive checksum result, source files, extraction status, and optional remote identifier provenance.

## Runtime Check

Before launching a workflow, CBIcall resolves `DATADIR` from the selected Bash `env.sh` or Snakemake/Nextflow/Cromwell `config.yaml`.

If bundle metadata exists beside `DATADIR`, CBIcall validates it:

| File | Runtime check |
| --- | --- |
| `cbicall-resource-id.json` | Resource key must match the selected `resource`; SHA-256 must match the catalog when pinned. |
| `cbicall-resource-installation.json` | Installed resource key must match the selected `resource`; the manifest catalog entry must match the local catalog fingerprint. |

This check is intentionally small. It validates the installed bundle identity without hashing the full resource archive on every run.
