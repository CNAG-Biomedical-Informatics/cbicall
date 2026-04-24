## Naming Conventions

CBIcall supports more than one FASTQ naming style.

!!! info "Short version"
    CBIcall accepts both **Illumina-style** FASTQ names and **run accession-based** FASTQ names such as those commonly seen in **1000 Genomes / ENA / SRA** datasets.

    The main workflow requirement is that paired-end files can be matched as `_R1_` and `_R2_`.

    The historical `_ex` / `_wg` convention is mainly a **legacy sample-directory naming** convention used for compatibility with older workflows such as **GATK 3.5** and **MTOOLBox**.

---

## What Is Required

!!! success "Required by active workflows"
    For current single-sample workflows, the important parts are:

    - paired FASTQs must follow the `_R1_` / `_R2_` convention
    - FASTQs for the same sample should be grouped under the same sample directory
    - FASTQ prefixes may differ within that directory

This means the FASTQ filename does **not** need to contain:

- the sample ID
- `_ex`
- `_wg`

as long as read pairing is still unambiguous.

---

## Supported FASTQ Styles

### Illumina-style naming

Typical format:

```text
SampleName_SampleNumber_Lane_Read_001.fastq.gz
```

Example:

```text
NA12878_S1_L001_R1_001.fastq.gz
NA12878_S1_L001_R2_001.fastq.gz
```

This format is supported.

### 1000 Genomes / ENA / SRA-style naming

Run accession-based prefixes are also supported, as long as mate pairing still uses `_R1_` / `_R2_`.

Example:

```text
HG00096/
├── ERR3239276_S1_L001_R1_001.fastq.gz
├── ERR3239276_S1_L001_R2_001.fastq.gz
├── ERR3239277_S1_L002_R1_001.fastq.gz
└── ERR3239277_S1_L002_R2_001.fastq.gz
```

In this layout, the FASTQ prefix is the run accession, and the sample is defined by the sample directory rather than by the FASTQ filename.

---

## Legacy `_ex` / `_wg` Convention

!!! warning "Use this for legacy compatibility"
    Keep the historical `_ex` / `_wg` convention if your project may use:

    - **GATK 3.5**
    - **MTOOLBox / mtDNA workflows**
    - older runs that already follow the original CNAG layout

In those cases, the important part is usually the **sample directory name**, not the FASTQ filename itself.

### Internal CNAG example

Project directory:

```text
CN99999_exome
```

- `CN99999` is the family / project identifier
- `_exome` indicates the analysis type

Sample directory:

```text
CN9999901P_ex
```

- `CN9999901` is the sample ID
- `P` is the role code
- `_ex` indicates exome

Common role codes:

- `P` = proband / patient
- `F` = father
- `M` = mother

Legacy-compatible directory layout:

```text
CN99999_exome/
└── CN9999901P_ex/
    ├── CN9999901P_ex_S1_L001_R1_001.fastq.gz
    └── CN9999901P_ex_S1_L001_R2_001.fastq.gz
```

Legacy-compatible directory layout with run-based FASTQ names:

```text
CN99999_exome/
└── CN9999901P_ex/
    ├── ERR3239276_S1_L001_R1_001.fastq.gz
    ├── ERR3239276_S1_L001_R2_001.fastq.gz
    ├── ERR3239277_S1_L002_R1_001.fastq.gz
    └── ERR3239277_S1_L002_R2_001.fastq.gz
```

So `_ex` / `_wg` should be thought of mainly as a directory convention for legacy compatibility.

---

## Equivalent FASTQ Examples

For active single-sample workflows, these filenames are equally valid:

- `CNAG99901P_ex_S2_L001_R1_001.fastq.gz`
- `CNAG99901P_S2_L001_R1_001.fastq.gz`

What matters is that the matching mate exists as `_R2_` and that all FASTQs for the sample are kept in the same sample directory.

---

## Notes

- `_ex` means exome sequencing
- `_wg` means whole-genome sequencing
- Illumina-style FASTQ naming follows the standard described here:
  <https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm>
