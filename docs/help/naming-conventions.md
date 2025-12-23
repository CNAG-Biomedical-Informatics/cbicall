## Scope and Tool Compatibility

This document defines **naming conventions required to run analyses in `cohort` mode with legacy tools**:

- **GATK 3.5**
- **MTOOLBox**

!!! warning "Legacy tools only"
    The directory and subdirectory naming conventions described in this document are **mandatory for GATK 3.5 and MTOOLBox**.
    They are **not required for GATK ≥ 4.6**.

    With **GATK 4.6**, sample directory names are **ignored by the pipeline** and may be **any arbitrary string**
    (e.g. `NEUMA640ZBR`).

---

## Directory Naming

### Format

```
[ProjectCode]_[SampleType]
```

- `ProjectCode`  
  Exactly **7 characters**: `[a-zA-Z0-9]`  
  Example: `CN99999`

- `SampleType`  
  One of:
    - `exome`
    - `wgs`

### Examples

```
CN99999_exome
CN99999_wgs
```

---

## Subdirectory Naming

### Format

```
[ProjectCode][SampleID][Role]_[SampleTypeShort]
```

- `ProjectCode` — 7 characters (`[a-zA-Z0-9]`)
- `SampleID` — 2 characters (e.g. `01`)
- `Role` — 1 character (`P`, `F`, `M`)
- `SampleTypeShort` — must be `ex`

### Example

```
CN9999901F_ex
```

---

## FASTQ Naming Convention

This FASTQ naming convention applies **to all tools**, including **GATK 3.5** and **GATK 4.6**.

It is based on the Illumina specification:
<https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm>

### Sequencing-type suffix

- `_ex` — exome sequencing
- `_wg` — whole genome sequencing

### Example

```
CN9999901P_ex_S5_L001_R1_001.fastq.gz
```

---

## Expected Directory Layout

!!! warning "Required for cohort mode (legacy)"
    This directory layout is **required** when running **GATK 3.5** or **MTOOLBox** in `cohort` mode.

```
CN99999_exome/
└── CN9999901P_ex/
    ├── CN9999901P_ex_S1_L001_R1_001.fastq.gz
    └── CN9999901P_ex_S1_L001_R2_001.fastq.gz
```
