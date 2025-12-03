## Directory Naming

- Format: `[ProjectCode]_[SampleType]`

    - - `ProjectCode`: Exactly 7 characters \[a-zA-Z0-9\] (e.g., `CN99999`)
    - - `SampleType`: Must be `exome` (5 characters) or `wgs` (3 characters)

    Example (exome):

        CN99999_exome

    Total: 7 + 1 + 5 = 13 characters.

    Example (wgs):

        CN99999_wgs

    Total: 7 + 1 + 3 = 11 characters.


## Subdirectory Naming

- Format: `[ProjectCode][SampleID][Role]_[SampleTypeShort]`

    - - `ProjectCode`: 7 characters (\[a-zA-Z0-9\] e.g., `CN99999`)
    - - `SampleID`: 2 characters (e.g., `01`)
    - - `Role`: 1 character (e.g., `P` for Proband, `F` for Father, `M` for Mother)
    - - `SampleTypeShort`: Must be `ex` (2 characters)

    Example:

        CN9999901F_ex

    Total: 7 + 2 + 1 + 1 + 2 = 13 characters (excluding any file extension).

## FASTQ Naming Convention

This convention is adapted from the following document:

[Document](https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)

We have added a custom suffix to indicate sequencing type:

    - _ex for exome sequencing
    - _wg for whole genome sequencing

Example:

    CN9999901P_ex_S5_L001_R1_001.fastq.gz

In summary, you need to have something like this:

    CN99999_exome/CN9999901P_ex/CN9999901P_ex_S1_L001_R?_001.fastq.gz
