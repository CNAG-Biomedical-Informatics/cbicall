# End-to-end example (WES single-sample run)

!!! note "Prerequisites"
    Installation, reference bundles, and all dependencies must be completed beforehand.  

    [➡️ Installation](../installation/non-containerized.md){ .md-button .md-button--primary }

This example demonstrates how to run CBIcall on a real WES sample from FASTQ files through final VCF and QC outputs.

---

## 1. Prepare your FASTQ files


CBIcall expects paired-end FASTQ files with a shared prefix, for example:

```
# Project    / Sample (Proband WES)
CN00001_exome/CN0000101P_ex/
  CN0000101P_ex_S1_L001_R1_001.fastq.gz
  CN0000101P_ex_S1_L001_R2_001.fastq.gz
```

!!! Note "Note on nomenclature"

    Please see [this page](../help/naming-conventions.md).
---

## 2. Create a parameters file

Create a YAML file, e.g. `wes_single.yaml`:

```yaml
mode:            single
pipeline:        wes
workflow_engine: bash
gatk_version:    gatk4.6
sample:          CN00001_exome/CN0000101P_ex
cleanup_bam:     false
```

Notes:

- `mode` selects single-sample or cohort (joint genotyping).  
- `pipeline` switches between WES, WGS or mtDNA.  
- `workflow_engine` chooses the backend (bash or snakemake).  

---

## 3. Run CBIcall

```bash
bin/cbicall -p wes_single.yaml -t 4
```

- `-p` selects the YAML parameters file  
- `-t` sets the number of threads

---

## 4. Inspect outputs

After completion, you will find:

```
CN00001_exome/CN0000101P_ex/cbicall_bash_wes_single_gatk-4.6_*/
  01_bam/
  02_varcall/
  03_stats/
  logs/
```

Where:

- VCF files are stored in `02_varcalls/`  
- QC metrics (coverage, sample stats, sex prediction) are in `03_stats`  
- Logs for all pipeline steps are under `logs/`  

These files are ready for downstream analysis, annotation or integration with cohort-level studies.

---

For advanced parameters, multi-sample analyses, mtDNA workflows and troubleshooting, see the **Usage** and **FAQ** sections.
