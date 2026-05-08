# Naming Conventions

CBIcall supports both Illumina-style FASTQ names and run accession-based FASTQ names used in public datasets such as 1000 Genomes, ENA, and SRA.

:::tip[Main rule]
For active WES/WGS single-sample workflows, paired FASTQs must be matchable as `_R1_` and `_R2_`, and all FASTQs for the sample should live in the same sample directory.
:::

The FASTQ filename does not need to contain the sample ID, `_ex`, or `_wg` as long as mate pairing is unambiguous.

## What CBIcall Uses

| Item | Required? | Why it matters |
| --- | --- | --- |
| `_R1_` / `_R2_` in FASTQ names | Yes | Used to pair reads. |
| One sample per sample directory | Yes | The sample directory defines the analysis unit. |
| Sample ID inside FASTQ filename | No | Useful for humans, but not required by current single-sample workflows. |
| `_ex` / `_wg` suffix | Legacy-dependent | Useful for older CNAG layouts, GATK 3.5, and mtDNA workflows. |

## Supported FASTQ Styles

### Illumina-Style Names

```text
NA12878_S1_L001_R1_001.fastq.gz
NA12878_S1_L001_R2_001.fastq.gz
```

This is the standard pattern for many local sequencing runs.

### Run Accession-Based Names

```text
HG00096/
├── ERR3239276_S1_L001_R1_001.fastq.gz
├── ERR3239276_S1_L001_R2_001.fastq.gz
├── ERR3239277_S1_L002_R1_001.fastq.gz
└── ERR3239277_S1_L002_R2_001.fastq.gz
```

Here, the FASTQ prefix is the run accession. The sample is defined by the directory, not by the FASTQ basename.

## Legacy `_ex` / `_wg` Layout

Keep the historical `_ex` / `_wg` convention when a project may use:

- GATK 3.5 workflows
- MToolBox / mtDNA workflows
- older CNAG run layouts

:::info[Directory convention]
In legacy-compatible layouts, the sample directory name is more important than the FASTQ filename itself.
:::

Project directory:

```text
CN99999_exome
```

Sample directory:

```text
CN9999901P_ex
```

| Part | Meaning |
| --- | --- |
| `CN99999` | Family or project identifier. |
| `01` | Individual/sample index within the project convention. |
| `P` | Role code, for example proband/patient. |
| `_ex` | Exome sample directory suffix. |

Common role codes:

| Code | Meaning |
| --- | --- |
| `P` | Proband / patient |
| `F` | Father |
| `M` | Mother |

Legacy-compatible layout:

```text
CN99999_exome/
└── CN9999901P_ex/
    ├── CN9999901P_ex_S1_L001_R1_001.fastq.gz
    └── CN9999901P_ex_S1_L001_R2_001.fastq.gz
```

The same layout also works with run-based FASTQ names:

```text
CN99999_exome/
└── CN9999901P_ex/
    ├── ERR3239276_S1_L001_R1_001.fastq.gz
    ├── ERR3239276_S1_L001_R2_001.fastq.gz
    ├── ERR3239277_S1_L002_R1_001.fastq.gz
    └── ERR3239277_S1_L002_R2_001.fastq.gz
```

## Equivalent Examples

For current single-sample WES/WGS workflows, both examples are valid:

```text
CNAG99901P_ex_S2_L001_R1_001.fastq.gz
CNAG99901P_S2_L001_R1_001.fastq.gz
```

What matters is that the matching mate exists as `_R2_` and that all FASTQs for the sample are grouped together.

## Next Steps

- Build a YAML file with [Configuration Reference](configuration-reference).
- Run a test with [Quickstart](../usage/quickstart).
- Inspect generated files with [Outputs](outputs).
