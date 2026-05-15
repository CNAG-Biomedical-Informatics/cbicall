import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# End-to-end examples (GATK 4.6)

> **Prerequisites**
    Installation, reference bundles, and all dependencies must be completed beforehand.  

    [Installation](../installation/non-containerized)

---

<Tabs groupId="workflow-mode">
<TabItem value="wes-single-sample-run" label="WES single-sample run" default>

This example demonstrates how to run CBIcall on a real WES sample from FASTQ files through final VCF and QC outputs.

## 1. Prepare your FASTQ files


CBIcall expects paired-end FASTQ files with a shared prefix, for example:
```
# Project    / Sample (Proband WES)
CNAG999_exome/CNAG99901P_ex/
  CNAG99901P_ex_S1_L001_R1_001.fastq.gz
  CNAG99901P_ex_S1_L001_R2_001.fastq.gz
```

> **Note on nomenclature**
    Please see [this page](../help/naming-conventions).
---

## 2. Create a parameters file

Create a YAML file, e.g. `wes_single.yaml`:

```yaml
mode:            single
pipeline:        wes
workflow_engine: bash
gatk_version:    gatk-4.6
input_dir:       CNAG999_exome/CNAG99901P_ex
genome:          b37
cleanup_bam:     false
```

Notes:

- `mode` selects single-sample or cohort (joint genotyping).  
- `pipeline` switches between WES, WGS or mtDNA.  
- `workflow_engine` chooses the backend (bash or snakemake).  
- See [Configuration Reference](../help/configuration-reference) for all YAML keys and supported combinations.

> **How can I perform WGS?**
    Simply change the parameter `pipeline` to `wgs`. Like this:

    ```yaml
    mode:            single
    pipeline:        wgs
    workflow_engine: bash
    gatk_version:    gatk-4.6
    input_dir:       CNAG999_exome/CNAG99901P_ex
    genome:          b37
    cleanup_bam:     false
```

---

## 3. Run CBIcall

```bash
bin/cbicall run -p wes_single.yaml -t 4
```

- `-p` selects the YAML parameters file  
- `-t` sets the number of threads


You should see something like this on the screen:

```bash
CBIcall 1.0.0
  Executable   => .../cbicall/bin/cbicall
  Workflow     => bash -> wes -> single
  Genome       => b37
  Threads      => 4
  Project      => .../CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_b37_gatk-4.6_177447031761843
  Run ID       => 177447031761843

Inputs
  Param file   => wes_single.yaml
  Input dir    => .../input/CNAG999_exome/CNAG99901P_ex
  Sample map   => (undef)
  GATK         => gatk-4.6

Resolved
  Entrypoint   => .../bash/gatk-4.6/wes_single.sh
  Log          => /media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/examples/input/CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_b37_gatk-4.6_177447031761843/bash_wes_single_b37_gatk-4.6.log

Running
  Workflow     => bash -> wes -> single
  This workflow may take a while depending on input size and pipeline.

Completed
  Status       => Finished successfully
  Elapsed      => 1m 30s
  Log          => /media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/examples/input/CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_b37_gatk-4.6_177447031761843/bash_wes_single_b37_gatk-4.6.log
Do Widzenia
```
---

## 4. Inspect outputs

After completion, you will find:
```
CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_b37_gatk-4.6_*/
  01_bam/
  02_varcall/
  03_stats/
  logs/
```

Where:

- VCF files are stored in `02_varcall/`  
- QC metrics (coverage, sample stats, sex prediction) are in `03_stats`  
- Logs for all pipeline steps are under `logs/`  

:::tip[What you get]
- Final VCF for interpretation: `02_varcall/<id>.hc.QC.vcf.gz`
- gVCF for cohort joint genotyping: `02_varcall/<id>.hc.g.vcf.gz`
- Run metadata: `log.json`

See [Outputs](../help/outputs) for the full file reference.
:::

---

For advanced parameters, multi-sample analyses, mtDNA workflows and troubleshooting, see the **Usage** and **FAQ** sections.

</TabItem>
<TabItem value="wes-cohort-run" label="WES cohort run">

> **Important**
    In order to run a `cohort` based calculation you first have to create `GVCF` for each sample. This is being done by running `wes` mode `single`.


## 1. Create a sample map file like the one we display below:


[View source](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/examples/input/sample_map.tsv)

**GATK** needs **absolute paths** for the files.

---

> **Scaling joint genotyping for large cohorts**
    For very large cohorts (hundreds to thousands of samples), the joint
    genotyping step can become computationally demanding when executed as a
    single job.

    A common strategy is to split the analysis by chromosome and run one job
    per chromosome through the HPC scheduler. This reduces memory usage per
    job and allows parallel execution across compute nodes.

    After all chromosomes have finished, the resulting VCF files can be merged (if needed)
    into a single cohort callset.

    Example:

    ```bash
    bcftools concat -Oz -o cohort_merged.vcf.gz chr*.vcf.gz
    bcftools index cohort_merged.vcf.gz
```

    In this example, `chr*.vcf.gz` simply represents a set of per-chromosome
    VCF files (e.g., `chr1.vcf.gz`, `chr2.vcf.gz`, …). The naming pattern is
    arbitrary and should be adapted to the filenames generated by your
    workflow.

## 2. Create a parameters file

Create a YAML file, e.g. `wes_cohort.yaml`:

```yaml
mode:            cohort
pipeline:        wes
workflow_engine: bash
gatk_version:    gatk-4.6
genome:          b37
sample_map:      ./sample_map.tsv
```

---

## 3. Run CBIcall

```bash
bin/cbicall run -p wes_cohort.yaml -t 4
```

- `-p` selects the YAML parameters file  
- `-t` sets the number of threads

---

## 4. Inspect outputs

After completion, you will find:
```
cbicall_bash_wes_cohort_b37_gatk-4.6_*/
  02_varcall/
  logs/
```

Where:

- Final VCF files are stored in `02_varcall/`  
- Logs for all pipeline steps are under `logs/`  

:::tip[What you get]
- Final joint VCF: `02_varcall/cohort.gv.QC.vcf.gz`
- GenomicsDB workspace and raw cohort VCF in `02_varcall/`
- Run metadata: `log.json`

See [Outputs](../help/outputs) for the full file reference.
:::

</TabItem>
</Tabs>
> **Do you have examples in how to run CBIcall programatically?**
    Yes, you can find examples at [https://github.com/CNAG-Biomedical-Informatics/cbicall/tree/main/examples/scripts](https://github.com/CNAG-Biomedical-Informatics/cbicall/tree/main/examples/scripts).

> **Any suggestions for performing annotation?**
    We recommend using [beacon2-cbi-tools](https://github.com/CNAG-Biomedical-Informatics/beacon2-cbi-tools). This tool allows you not only to annotate data, but also to convert it into a data exchange format compatible with the [Beacon v2 API](https://genomebeacons.org/).
