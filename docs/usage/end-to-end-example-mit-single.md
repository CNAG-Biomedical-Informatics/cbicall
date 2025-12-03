# End-to-end example (MIT single-sample run)

!!! note "Prerequisites"
    Installation, reference bundles, and all dependencies must be completed beforehand.  

    [➡️ Installation](../installation/non-containerized.md){ .md-button .md-button--primary }

!!! warning "Architecture"
    MToolBox supports **x86_64 only**. ARM-based systems, including Apple Silicon (M1/M2/M3), are **not supported**.

---

## 1. Before running mtDNA calling you must have a `bam` file coming from wes/wgs

!!! example "Does it matter if I ran WES/WGS with GATK 3.5 or GATK 4.6?"
    No. CBICall will detect and use the `bam` files produced by either version.  
    Just make sure that `bam` files are available — **FASTQ input is not supported**.

CBIcall expects BAM file infrom a previous run:

```
CN00001_exome
└── CN0000101P_ex  <--- ID taken from here
    └── cbicall_bash_mit_single_gatk-3.5_* <- The script expects that you are submitting the job from inside this directory
```

!!! Note "Note on nomenclature"

    Please see [this page](../help/naming-conventions.md).
---

## 2. Create a parameters file

Create a YAML file, e.g. `mit.yaml`:

```yaml
mode:            single
pipeline:        mit
workflow_engine: bash
gatk_version:    gatk3.5
sample:          CN00001_exome/CN0000101P_ex
```

Notes:

- `mode` selects single-sample or cohort (joint genotyping).  
- `pipeline` switches between WES, WGS or mtDNA.  
- `workflow_engine` chooses the backend (bash or snakemake).  

---

## 3. Run CBIcall

```bash
bin/cbicall -p wes_mit.yaml -t 4
```

- `-p` selects the YAML parameters file  
- `-t` sets the number of threads

---

## 4. Inspect outputs

After completion, you will find:

```
CN00001_exome/CN0000101P_ex/cbicall_bash_wes_mit_gatk-3.5_*/
  01_mtoolbox/
  02_browser/
```

Where:

- VCF files are stored in `01_mtoolbox/`  

## 5. Visualize variants with web browser

Foo Bar


---

For advanced parameters, multi-sample analyses, mtDNA workflows and troubleshooting, see the **Usage** and **FAQ** sections.
