import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# End-to-end examples (MToolBox)

> **Prerequisites**
    Installation, reference bundles, and all dependencies must be completed beforehand.  

    [Installation](../installation/non-containerized)

> **Architecture**
    MToolBox supports **x86_64 only**. ARM-based systems, including Apple Silicon (M1/M2/M3), are **not supported**.

---

<Tabs groupId="workflow-mode">
<TabItem value="mit-single-sample-run" label="MIT single-sample run" default>

## 1. Before running mtDNA calling you must have a BAM file from WES/WGS

> **Does it matter if I ran WES/WGS with GATK 3.5 or GATK 4.6?**
    No. CBIcall will detect and use the `bam` files produced by either version.  
    Just make sure that `bam` files are available — **FASTQ input is not supported**.

CBIcall expects a BAM file from a previous WES/WGS run:
```
CNAG999_exome
└── CNAG99901P_ex  <--- ID taken from here
    └── *cbicall_bash_w?s_single_gatk-* <- The script expects that you have a BAM file inside this directory
```

> **Note on nomenclature**
    Please see [this page](../help/naming-conventions).
---

## 2. Create a parameters file

Create a YAML file, e.g. `mit_single.yaml`. 

> **Important**
    Please make sure you use the same value for the key `sample` that you used for WES/WGS.

Example:

```yaml
mode:            single
pipeline:        mit
workflow_engine: bash
input_dir:       CNAG999_exome/CNAG99901P_ex
```

See [Configuration Reference](../help/configuration-reference) for all YAML keys and supported combinations.

---

## 3. Run CBIcall

```bash
bin/cbicall run -p mit_single.yaml -t 4
```

- `-p` selects the YAML parameters file  
- `-t` sets the number of threads

---

## 4. Inspect outputs

After completion, you will find:
```
CNAG999_exome/CNAG99901P_ex/cbicall_bash_mit_single_rsrs_gatk-3.5_*/
  01_mtoolbox/
  02_browser/
```

:::tip[What you get]
- Final mtDNA report: `01_mtoolbox/mit_prioritized_variants.txt`
- mtDNA VCF: `01_mtoolbox/VCF_file.vcf`
- Browser report: `02_browser/<run-id>.html`

See [Outputs](../help/outputs) for the full file reference.
:::

## 5.  Visualize variants in the browser

Please see:
```
02_browser/README.txt
```

The results are reported both as an HTML table and as downloadable files.

> **See snapshot**
    ![browser](/img/browser.png)

### Downloadable files:

* **mtDNA JSON** A JSON file with the results from `mit_prioritized_variants.txt`.
* **Report:** A tsv file including all the annotations for each variant. Name of the file `mit_prioritized_variants.txt`.
* **Haplog:** A tsv file including the predicted [haplogroup](https://en.wikipedia.org/wiki/Human_mitochondrial_DNA_haplogroup) for each sample. Name of the file `mt_classification_best_results.csv`.
* **VCF:** A text file consisting of all the variants in the VCF format. Name of the file `VCF_file.vcf`.

### HTML table:

The CBIcall mtDNA variation browser displays a browsable table consisting of the most relevant fields relative to the variant annotation:

* **Sample**: The full name of each sample.
* **Locus**: The location on the mitochondrial chromosome.
* **Variant_Allele**: The position in the mitochondrial chromosome + the alternative allele format.
* **Ref**: The reference allele (mitochondrial reference genome: RSRS).
* **Alt**: The alternative allele(s).
* **Aa_change**: The amino acid change if the variant falls in a coding region.
* **GT**: Genotype. 0:Ref, ≥1:Alt(s).
* **Depth**: The number of times this position is covered by reads.
* **Heterop_Frac**: The heteroplasmic fraction. Note that the confidence interval can be retrieved from the downloadable VCF file.
* **Other**: For other fields please consult [MToolBox's manual](https://github.com/mitoNGS/MToolBox/wiki/Output-files).

> **Filtered variants**
    The table shows **pre-filtered variants**. Variants were excluded if:

    - **HF ≤ 0.30** (maximum HF observed in **any** sample)
    - **1000 Genomes frequency ≥ 0.01**
    - **Not present in the input VCF**

    By default, variants with **missing HF values (`NA`,`N/A`,`.`) are excluded**.
    Use the `--keep-missing-hf` option to retain them.

---

For advanced parameters, multi-sample analyses, mtDNA workflows and troubleshooting, see the **Usage** and **FAQ** sections.

</TabItem>
<TabItem value="mit-cohort-run" label="MIT cohort run">

## 1. Before running mtDNA calling you must have BAM files from WES/WGS

> **Does it matter if I ran WES/WGS with GATK 3.5 or GATK 4.6?**
    No. CBIcall will detect and use the `bam` files produced by either version.  
    Just make sure that `bam` files are available — **FASTQ input is not supported**.

CBIcall expects BAM files from previous WES/WGS runs:
```
CNAG999_exome
└── CNAG99901P_ex  <--- ID taken from here
    └── *cbicall_bash_w?s_single_gatk-* <- The script expects that you have a BAM file inside this directory
    CNAG99902M_ex  <--- ID taken from here
    └── *cbicall_bash_w?s_single_gatk-* <- The script expects that you have a BAM file inside this directory
```

> **Note on nomenclature**
    Please see [this page](../help/naming-conventions).
---

## 2. Create a parameters file

Create a YAML file, e.g. `mit_cohort.yaml`:

```yaml
mode:            cohort
pipeline:        mit
workflow_engine: bash
gatk_version:    gatk-3.5
input_dir:       CNAG999_exome
```

See [Configuration Reference](../help/configuration-reference) for all YAML keys and supported combinations.

---

## 3. Run CBIcall

```bash
bin/cbicall run -p mit_cohort.yaml -t 4
```

- `-p` selects the YAML parameters file  
- `-t` sets the number of threads

---

## 4. Inspect outputs

After completion, you will find:
```
CNAG999_exome/cbicall_bash_mit_cohort_rsrs_gatk-3.5*
  01_mtoolbox/
  02_browser/
```

:::tip[What you get]
- Final joint mtDNA report: `01_mtoolbox/mit_prioritized_variants.txt`
- Joint mtDNA VCF: `01_mtoolbox/VCF_file.vcf`
- Browser report: `02_browser/<run-id>.html`

See [Outputs](../help/outputs) for the full file reference.
:::

## 5.  Visualize variants in the browser

Please see:
```
02_browser/README.txt
```

The results are reported both as an HTML table and as downloadable files.

> **See snapshot**
    ![browser](/img/browser-cohort.png)

### Downloadable files:

* **mtDNA JSON** A JSON file with the results from `mit_prioritized_variants.txt`.
* **Report:** A tsv file including all the annotations for each variant. Name of the file `mit_prioritized_variants.txt`.
* **Haplog:** A tsv file including the predicted [haplogroup](https://en.wikipedia.org/wiki/Human_mitochondrial_DNA_haplogroup) for each sample. Name of the file `mt_classification_best_results.csv`.
* **VCF:** A text file consisting of all the variants in the VCF format. Name of the file `VCF_file.vcf`.

### HTML table:

The CBIcall mtDNA variation browser displays a browsable table consisting of the most relevant fields relative to the variant annotation:

* **Sample**: The full name of each sample.
* **Locus**: The location on the mitochondrial chromosome.
* **Variant_Allele**: The position in the mitochondrial chromosome + the alternative allele format.
* **Ref**: The reference allele (mitochondrial reference genome: RSRS).
* **Alt**: The alternative allele(s).
* **Aa_change**: The amino acid change if the variant falls in a coding region.
* **GT**: Genotype. 0:Ref, ≥1:Alt(s).
* **Depth**: The number of times this position is covered by reads.
* **Heterop_Frac**: The heteroplasmic fraction. Note that the confidence interval can be retrieved from the downloadable VCF file.
* **Other**: For other fields please consult [MToolBox's manual](https://github.com/mitoNGS/MToolBox/wiki/Output-files).

> **Filtered variants**
    The table shows **pre-filtered variants**. Variants were excluded if:

    - **HF ≤ 0.30** (maximum HF observed in **any** sample)
    - **1000 Genomes frequency ≥ 0.01**
    - **Not present in the input VCF**

    By default, variants with **missing HF values (`NA`,`N/A`,`.`) are excluded**.
    Use the `--keep-missing-hf` option to retain them.

</TabItem>
</Tabs>
> Genetic data interpretation disclaimer: review the project disclaimer before clinical or diagnostic interpretation.
    ---
