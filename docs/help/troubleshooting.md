# Common errors and troubleshooting

## Installation

### External data not found

**Error message (example from logs)**
```
/usr/bin/bash: line 9: /media/mrueda/2TBS/NGSutils/gatk/gatk-4.6.2.0/gatk: No such file or directory
```

**Cause**

The `DATADIR` variable is not correctly configured.

**Solution**

Update the `DATADIR` variable in the following files:

```bash
workflows/bash/gatk-4.6/env.sh
workflows/snakemake/gatk-4.6/config.yaml
```

Set it to your data directory (e.g., `/cbicall-data` or your mounted path).

---

## GATK and Picard errors  
*(wes_single.sh or wes_cohort.sh)*

### NaN LOD value during recalibration

**Error message**
```
NaN LOD value assigned
```

**Cause**

Too few INDEL variants (e.g., <8000) are available to train a reliable VQSR model.

**Solution**

Modify the pipeline to skip VQSR when the INDEL count is below the threshold.

- Increase or enforce the minimum INDEL threshold
- Rerun only the affected samples

This prevents model training on insufficient data.

---

### Not enough columns in dbSNP line

**Error message**
```
there aren't enough columns for line ... dbsnp_137.hg19.vcf
```

**Cause**

Malformed or truncated lines in the dbSNP VCF file.

**Solution**

- Locate the problematic line in the dbSNP VCF
- Remove or correct it
- Document the modification in a local README or changelog

---

### Error parsing text SAM file

**Error message**
```
Error parsing text SAM file. Not enough fields; File /dev/stdin; Line ...
```

**Cause**

Secondary or supplementary alignments (common in SRA/dbGaP data) can introduce invalid or conflicting records when piped directly into Picard/GATK.

**Solution**

Filter out secondary (0x100) and supplementary (0x800) alignments before adding read groups.

In `wes_single.sh`, ensure the alignment pipe includes filtering:

```bash
bwa mem -M -t "$THREADS" "$REFGZ" "$R1" "$R2" \
  | samtools view -bSh -F 0x900 - \
  | gatk AddOrReplaceReadGroups ...
```

---

## MToolBox errors and mtDNA-specific issues

### Unsupported `N` CIGAR operations

**Symptom**

MToolBox fails due to unsupported `N` operations in CIGAR strings.

**Solution**

Add the following flag in `MToolBox.sh` (around the main alignment or SAM processing step):

```
--filter_reads_with_N_cigar
```

This removes reads with `N` CIGAR operations before downstream processing.

---

### Low coverage and unreliable heteroplasmy fractions

**Symptom**

- Very low mtDNA coverage
- Noisy or unstable heteroplasmy fraction (HF) estimates

**Guideline**

- Below ~10× mtDNA coverage, HF estimates are unreliable
- At higher coverage, HF values are generally robust

**Recommended actions**

- Flag samples with <10× median mtDNA coverage
- Interpret HF cautiously in low-coverage samples
- Exclude such samples from HF-based analyses if necessary
- Consider resequencing if mtDNA analysis is critical

---

## Variant interpretation

### Expected de novo rates in trios

!!! note
    Typical de novo variant rates:

    - Probands: ~1%
    - Parents: ~10%

    Significant deviations may indicate technical or pipeline issues and should be investigated.
