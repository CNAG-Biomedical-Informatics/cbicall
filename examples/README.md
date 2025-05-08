# Examples README

This directory contains data to verify that the variant calling pipeline functions correctly.

---

## Input:

The input is pre-generated for you in the `input` directory.

### About Nomenclature

It is extremely important that directories and files follow a standardized nomenclature with the same number of characters. For example:

```bash
CNAG999_exome/  
├── CNAG99901P_ex  
```

In particular, apart from the number of characters, the `_exome` and `_ex`

Then, the `fastq` files should look like this:

```
CNAG99901F_ex_S2_L001_R1_001.fastq.gz  
CNAG99901F_ex_S2_L001_R2_001.fastq.gz  
...  
```

### Creating additional input examples

Note: Adjust paths to your local setup. All necessary files are provided with the [external data](../README.md#INSTALLATION).

```bash
DBDIR=/media/mrueda/2TBS/Databases
BUNDLE=$DBDIR/GATK_bundle/b37
REF=$BUNDLE/references_b37_Homo_sapiens_assembly19.fasta
WGSIM=/media/mrueda/2TBS/NGSutils/wgsim

$WGSIM -S 42 -N 50000 -1 150 -2 150 -r 0.001 -e 0.02 -R 0.001 -X 0.001 \
  $REF sim_R1.fastq sim_R2.fastq > wgsim.log
```

| Parameter | Description |
|-----------|-------------|
| `-S 42`   |  Seed       |
| `-N 1000000` | Generates **1,000,000 read pairs** |
| `-1 150` | Read length of **150bp** for read 1 |
| `-2 150` | Read length of **150bp** for read 2 |
| `-r 0.001` | Mutation rate (**0.1% of bases mutated**) |
| `-e 0.02` | Sequencing error rate (**2% of bases will have errors**) |
| `-R 0.001` | Fraction of reads with an **indel** |
| `-X 0.001` | Indel rate per base |
| `ref.fasta` | Your reference genome (input) |
| `sim_R1.fastq` | Output file for **forward reads (R1)** |
| `sim_R2.fastq` | Output file for **reverse reads (R2)** |

To follow the suggested nomenclature:

```bash
mkdir -p CNAG999_exome/CNAG99901P_ex
gzip -c sim_R1.fastq > CNAG999_exome/CNAG99901P_ex/CNAG99901F_ex_S2_L001_R1_001.fastq.gz
gzip -c sim_R2.fastq > CNAG999_exome/CNAG99901P_ex/CNAG99901F_ex_S2_L001_R2_001.fastq.gz
```

## Run `cbicall`

#### Wes Single

To test your installation, execute:

```bash
cd input
../../bin/cbicall -p param.yaml -t 4
```

This should complete in ~15 minutes.

---

##### Test Results

Check results using:

```bash
wc -l CNAG999_exome/CNAG99901P_ex/cbicall_/02_varcall/*QC*vcf
```

Expected result: ~128 lines, containing exactly 1 variant.

```
CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CNAG99901rPF
3	30686407	.	G	T	15.85	LowCoverage;LowQual;VeryLowQual	AC=2;AF=1.00;AN=2;DP=2;Dels=0.00;ExcessHet=3.0103;FS=0.000;HaplotypeScore=0.0000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=7.93;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:43,6,0
```

---

> **Note:** We will not be able to determine the sex of the sample as it does not contain enogh reads

---

#### WES Cohort

Generate larger FASTQ files using `wgsim` for multiple samples:

- `CNAG999_exome/CNAG99901F_ex`
- `CNAG999_exome/CNAG99901M_ex`
- (additional samples as necessary)

---

## Scripts

In the directory `scripts` you will find examples of how to run `cbicall` in bulk.
