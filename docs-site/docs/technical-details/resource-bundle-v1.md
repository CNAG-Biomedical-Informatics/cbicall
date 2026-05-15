# Resource Bundle v1

CBIcall production workflows require an external resource bundle containing third-party executables, reference genomes, known-sites files, interval lists, and auxiliary databases.

The current bundle is selected in the run YAML with:

```yaml
resource_bundle: "cbicall-germline-resources-v1"
```

CBIcall resolves this key against the resource catalog:

```text
resources/cbicall-resource-catalog.json
```

The downloader uses the local catalog when it is present. If the script is used standalone and the local catalog is absent, it fetches the canonical catalog URL before selecting the bundle entry.

## Bundle Identity

| Field | Value |
| --- | --- |
| Catalog key | `cbicall-germline-resources-v1` |
| Bundle ID | `cbicall-germline-resources` |
| Bundle version | `1` |
| Status | `current` |
| CBIcall compatibility | `>=1.0,<2.0` |
| Archive name after setup | `cbicall-germline-resources-v1.tar.gz` |

:::info[Why this is explicit]
The bundle key is human-readable, while the catalog entry is machine-readable. `log.json` records both the selected key and a catalog fingerprint, so two runs can be checked for the same declared external dependency set.
:::

## Supported Workflows

| Engine | Pipeline | Mode | GATK version |
| --- | --- | --- | --- |
| `bash` | `wes` | `single` | `gatk-3.5` |
| `bash` | `wes` | `cohort` | `gatk-3.5` |
| `bash` | `wes` | `single` | `gatk-4.6` |
| `bash` | `wes` | `cohort` | `gatk-4.6` |
| `bash` | `wgs` | `single` | `gatk-4.6` |
| `bash` | `wgs` | `cohort` | `gatk-4.6` |
| `bash` | `mit` | `single` | `gatk-3.5` |
| `bash` | `mit` | `cohort` | `gatk-3.5` |
| `snakemake` | `wes` | `single` | `gatk-4.6` |
| `snakemake` | `wes` | `cohort` | `gatk-4.6` |

## Downloaded Files

The production bundle is distributed as a small identifier JSON, split archive parts, and a checksum file.

| File | Purpose |
| --- | --- |
| `cbicall-bundle-id.json` | Declares the bundle catalog key and is pinned by SHA-256 in the catalog. |
| `data.tar.gz.md5` | MD5 checksum file. The current bundle records the split archive parts. |
| `data.tar.gz.part-00` | Split archive part. |
| `data.tar.gz.part-01` | Split archive part. |
| `data.tar.gz.part-02` | Split archive part. |
| `data.tar.gz.part-03` | Split archive part. |
| `data.tar.gz.part-04` | Split archive part. |
| `data.tar.gz.part-05` | Split archive part. |

The setup utility verifies the files covered by `data.tar.gz.md5`, reassembles the parts into `data.tar.gz`, and then renames the verified archive to `cbicall-germline-resources-v1.tar.gz`.

An optional small remote identifier file can also be used:

```json
{"bundle": "cbicall-germline-resources-v1"}
```

When available, this file is named `cbicall-bundle-id.json`. Its SHA-256 can be pinned in the local catalog to confirm that the remote bundle declares the expected catalog key.

## Expected Layout

After extraction, `DATADIR` should contain:

```text
DATADIR/
  Databases/
  NGSutils/
```

The catalog maps these names to workflow variables:

| Variable | Meaning |
| --- | --- |
| `DATADIR` | Root of the installed external resource bundle. |
| `DBDIR` | `DATADIR/Databases` |
| `NGSUTILS` | `DATADIR/NGSutils` |

Workflow-specific files such as `env.sh` and Snakemake `config.yaml` expand these variables into concrete executable and reference paths.

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

The mtDNA workflows use MToolBox-related mitochondrial resources from the external bundle.

## Provenance in Runs

The selected bundle is stored in `log.json` under `config.resource_bundle`.

Example:

```json
{
  "key": "cbicall-germline-resources-v1",
  "bundle_id": "cbicall-germline-resources",
  "bundle_version": "1",
  "compatible": true,
  "fingerprint": "..."
}
```

Two runs used the same declared external dependency set when their `config.resource_bundle.fingerprint` values match.

## Installation Manifest

The setup utility also writes:

```text
cbicall-resource-installation.json
```

This local manifest records the installed bundle key, archive checksum result, source files, extraction status, and optional remote identifier provenance.
