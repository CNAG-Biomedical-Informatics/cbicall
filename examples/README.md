# Examples README

This directory contains sample datasets used to validate the
variant-calling pipeline.

------------------------------------------------------------------------

## Input

Pre-generated input files are provided in the `input` directory:

### Chromosome 22 (chr22)

-   **50,000 paired-end reads** simulated with
    **[wgsim](https://github.com/lh3/wgsim)**
-   Reads originate only from chromosome 22
> **Note:** Sex determination is not possible because X/Y reads are
    not included.

------------------------------------------------------------------------

### mtDNA Reads -- HG00152

Mitochondrial reads from the **[1000 Genomes Project (HG00152)](https://www.internationalgenome.org/data-portal/sample/HG00152)**: 
- **R1:** ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR769/SRR769545/SRR769545_1.fastq.gz
- **R2:** ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR769/SRR769545/SRR769545_2.fastq.gz

------------------------------------------------------------------------

## Nomenclature Requirements

Directory and file names must follow a fixed-length, standardized
structure. Example:

``` bash
CNAG999_exome/
├── CNAG99901P_ex/
```

FASTQ files should follow this pattern (`01P` = Proband 1):

    CNAG99901P_ex_S2_L001_R1_001.fastq.gz
    CNAG99901P_ex_S2_L001_R2_001.fastq.gz
    ...

------------------------------------------------------------------------

## Creating Additional Input Examples

Adjust paths based on your local setup. Reference data is provided (see
the [installation instructions](../README.md#INSTALLATION)).\
You must install **[wgsim](https://github.com/lh3/wgsim)** separately.

``` bash
DBDIR=/media/mrueda/2TBS/Databases
BUNDLE=$DBDIR/GATK_bundle/b37
REF=$BUNDLE/references_b37_Homo_sapiens_assembly19.fasta
WGSIM=/media/mrueda/2TBS/NGSutils/wgsim

$WGSIM -S 42 -N 50000 -1 150 -2 150 -r 0.001 -e 0.02 -R 0.001 -X 0.001   $REF sim_R1.fastq sim_R2.fastq > wgsim.log
```

### wgsim Parameters

| Parameter                      | Meaning                       | Example Value        |
|--------------------------------|-------------------------------|-----------------------|
| `-S`                           | Random seed                   | `42`                  |
| `-N`                           | Number of read pairs          | `50000`               |
| `-1`                           | Read length (R1)              | `150`                 |
| `-2`                           | Read length (R2)              | `150`                 |
| `-r`                           | Mutation rate                 | `0.001`               |
| `-e`                           | Sequencing error rate         | `0.02`                |
| `-R`                           | Fraction of reads with indels | `0.001`               |
| `-X`                           | Indel rate per base           | `0.001`               |
| `ref.fasta`                    | Reference genome              | `b37` or your own     |
| `sim_R1.fastq`, `sim_R2.fastq` | Output FASTQ files            | R1/R2                 |

To apply the standard nomenclature:

``` bash
mkdir -p CNAG999_exome/CNAG99901P_ex
gzip -c sim_R1.fastq > CNAG999_exome/CNAG99901P_ex/CNAG99901P_ex_S2_L001_R1_001.fastq.gz
gzip -c sim_R2.fastq > CNAG999_exome/CNAG99901P_ex/CNAG99901P_ex_S2_L001_R2_001.fastq.gz
```

------------------------------------------------------------------------

## Scripts

The `scripts` directory contains examples for running `cbicall` in bulk.
