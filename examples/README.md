# CBIcall example data

The compact paired-end fixture under
`input/CNAG999_exome/CNAG99901P_ex/` is used by the CBIcall integration
tests. It combines a simulated nuclear component with real mitochondrial
read sequences so that the same input can exercise the native WES and mtDNA
workflows. It is test data, not an analytical benchmark sample.

## Fixture composition

Each mate file contains 55,053 reads:

| Component | Read pairs | Origin |
| --- | ---: | --- |
| Nuclear | 50,000 | Chromosome 22 reads simulated from the GATK b37 reference with `wgsim` |
| Mitochondrial | 5,053 | HG00152 reads from ENA run SRR769545, selected by an earlier CBIcall/MToolBox execution |

### Nuclear reads

Chromosome 22 was extracted from
`references_b37_Homo_sapiens_assembly19.fasta`. The paired 150-bp reads were
then generated with the `wgsim` distributed with SAMtools 0.1.19:

```bash
samtools faidx references_b37_Homo_sapiens_assembly19.fasta 22 > chr22.fa

wgsim \
  -S 42 \
  -N 50000 \
  -1 150 \
  -2 150 \
  -r 0.001 \
  -e 0.02 \
  -R 0.001 \
  -X 0.001 \
  chr22.fa sim_R1.fastq sim_R2.fastq
```

The simulated read identifiers begin with `@22_`.

### Mitochondrial reads

The source sample is
[HG00152](https://www.internationalgenome.org/data-portal/sample/HG00152),
sequencing run SRR769545:

- [R1 FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR769/SRR769545/SRR769545_1.fastq.gz)
- [R2 FASTQ](https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR769/SRR769545/SRR769545_2.fastq.gz)

The full SRR769545 FASTQs were processed by the native CBIcall mtDNA
workflow. The historical construction record identifies run
`cbicall_bash_mit_single_gatk-3.5_648472991873040`; its MToolBox
`outmt1.fastq` and `outmt2.fastq` outputs supplied 5,053 mtDNA-enriched read
pairs to the fixture.

All 5,053 SRR769545 identifiers and nucleotide sequences in the fixture match
the corresponding source reads. The fixture records are not byte-identical to
the downloaded ENA FASTQs because they were taken from the MToolBox FASTQ
outputs: their headers are shortened and their quality strings reflect that
intermediate representation.

## Fixture assembly

The two components were concatenated by mate and compressed:

```bash
mkdir -p CNAG999_exome/CNAG99901P_ex

cat sim_R1.fastq outmt1.fastq | gzip -c \
  > CNAG999_exome/CNAG99901P_ex/CNAG99901P_ex_S2_L001_R1_001.fastq.gz

cat sim_R2.fastq outmt2.fastq | gzip -c \
  > CNAG999_exome/CNAG99901P_ex/CNAG99901P_ex_S2_L001_R2_001.fastq.gz
```

The reconstructed uncompressed streams match the repository fixtures exactly.
The current checksums are:

| File | SHA-256 of `.fastq.gz` | SHA-256 after decompression |
| --- | --- | --- |
| `CNAG99901P_ex_S2_L001_R1_001.fastq.gz` | `b1332b602dfe1e4e06c819c09c9316dcfd37bb8a56bb0fc0b04046ee10014a5e` | `018c30c973d6b9426efc56b7a5af03df09ebedba376fd1f08edc01aedc5b26ae` |
| `CNAG99901P_ex_S2_L001_R2_001.fastq.gz` | `169f90aa84d3f6c4d46e64a94fbf6a02780950c25fd96af69b96bc4c1811b4af` | `da28e07152f256eea185209eff9fc2509ebe833ac218e7c0e59440699bf68187` |

These fixture files first entered the repository in commit `37d1b6d` and have
the same content as the recovered construction copies.

## Layout and naming

```text
CNAG999_exome/
`-- CNAG99901P_ex/
    |-- CNAG99901P_ex_S2_L001_R1_001.fastq.gz
    `-- CNAG99901P_ex_S2_L001_R2_001.fastq.gz
```

The filenames intentionally follow the CBIcall paired-end naming convention.
Because the nuclear component contains chromosome 22 only, it is suitable for
compact workflow integration tests but not for genome-wide coverage or sex
inference validation.

The `scripts/` directory contains examples for running CBIcall in bulk.
