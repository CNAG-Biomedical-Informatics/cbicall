# Common errors and troubleshooting

!!! note "Expected de novo rates in trios"
    For trio analyses, approximate de novo variant rates are:

    - Probands: ~1%
    - Parents: ~10%

    Large deviations from these ranges may indicate technical or pipeline issues that warrant investigation.

---

## GATK and Picard errors  
*(wes_single.sh or wes_cohort.sh)*

### NaN LOD value assigned during recalibration

**Error message example**

`NaN LOD value assigned` during VariantRecalibrator or ApplyVQSR.

**Cause**

This typically occurs when there are too few INDEL variants (for example, fewer than about 8000) to train a robust negative model. The default minimum INDEL count threshold is 8000 in the VQSR step.

**Solution**

Increase the minimum INDEL count threshold in the relevant pipeline script so that VQSR is skipped for samples with low INDEL counts. Only rerun the affected samples.  
This prevents VariantRecalibrator from trying to build a model on too few variants.

---

### Not enough columns in dbSNP line

**Error message example**

`there aren't enough columns for line ... dbsnp_137.hg19.vcf`

**Cause**

One or more lines in the dbSNP VCF file do not conform to the expected VCF column structure (for example, a truncated or malformed record).

**Solution**

- Identify the problematic line in the dbSNP VCF.
- Remove or fix that line.
- Document the change in a local README or change log so the modification is traceable.

---

### Error parsing text SAM file

**Error message example**

`Error parsing text SAM file. Not enough fields; File /dev/stdin; Line 105120626...`

**Cause**

Some SRA or dbGaP datasets include duplicate or problematic reads. When piping BWA output directly into AddOrReplaceReadGroups, secondary and supplementary alignments can cause issues and lead to collisions or invalid lines as seen by Picard or GATK.

**Solution**

Remove secondary (0x100) and supplementary (0x800) alignments from the BWA stream before adding read groups.

In `wes_single.sh`, uncomment the filtering step in the alignment pipe, for example:

```bash
bwa mem -M -t "$THREADS" "$REFGZ" "$R1" "$R2"   | samtools view -bSh -F 0x900 -   | gatk AddOrReplaceReadGroups ...
```

This filtering prevents problematic alignments from reaching Picard or GATK and avoids the parsing error.

---

## MToolBox errors and mtDNA specific issues

### Unsupported N CIGAR operations

**Symptom**

MToolBox fails with an error related to unsupported N operations in CIGAR strings.

**Solution**

Add the `--filter_reads_with_N_cigar` flag in `MToolBox.sh` (around the main `bwa mem` or SAM tools invocations, typically near line 386 in your local copy).

This discards reads with N operations in the CIGAR string before downstream processing, avoiding MToolBox failures.

---

### Low coverage and unreliable heteroplasmy fractions

**Symptom**

- Very low coverage mtDNA samples.
- Heteroplasmic fraction (HF) estimates appear noisy or unreliable.

**Guideline**

- Below about 10x mtDNA coverage, HF estimates become unreliable and may be biologically meaningless.
- At higher coverage, HF values tend to be robust, even when coverage varies across samples.

**Recommended actions**

- Flag samples with less than 10x median mtDNA coverage for review.
- Interpret HF with caution in low coverage samples, or exclude them from HF based analyses.
- Consider resequencing or deeper coverage if mtDNA heteroplasmy is critical to the study.
