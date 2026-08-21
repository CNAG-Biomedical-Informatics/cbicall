import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# mtDNA pipelines

The bundled `cbicall-core` mtDNA workflows use mtDNA-only BAMs exported by
native GATK 4.6 WES/WGS single-sample runs. MToolBox then calls, annotates, and
prioritizes mtDNA variants.

| Mode | Use |
| --- | --- |
| `single` | Analyze one individual. |
| `cohort` | Analyze all eligible samples in one project together. |

:::warning[Architecture]
The bundled MToolBox workflow supports x86_64 Linux only.
:::

:::warning[Required input]
Run WES/WGS single-sample processing with `export_mtdna_bam: true` first.
mtDNA workflows require `04_mtdna_input/<sample>-DNA_MIT.bam` and its index;
they do not extract reads from a full BAM.
:::

The exported BAM inherits mitochondrial read placement from the preceding
BWA-MEM alignment. CBIcall uses a fixed BWA input batch size so changing the
workflow thread count does not change which reads are passed to MToolBox. See
the [WES/WGS single-sample pipeline](./wes-wgs-single) for details.

<Tabs groupId="workflow-mode">
<TabItem value="single" label="Single sample" default>

## Single-sample workflow

![mtDNA single-sample workflow](/img/diagram-mtdna-single.svg)

Set `input_dir` to a sample directory containing a completed native WES/WGS
single-sample run with an exported mtDNA BAM. CBIcall copies that BAM into the
MToolBox workspace and appends genotype (`GT`), depth (`DP`), and heteroplasmy
fraction values to the prioritized report.

```yaml
mode:             single
pipeline:         mit
workflow_backend: bash
software_stack:   gatk-3.5
input_dir:        CNAG999_exome/CNAG99901P_ex
```

</TabItem>
<TabItem value="cohort" label="Cohort">

## Cohort workflow

![mtDNA cohort workflow](/img/diagram-mtdna-cohort.svg)

Set `input_dir` to a project directory containing the sample directories and
their completed native WES/WGS single-sample runs. Each run must contain its
exported mtDNA BAM and index. CBIcall collects these files and runs MToolBox
jointly.

```yaml
mode:             cohort
pipeline:         mit
workflow_backend: bash
software_stack:   gatk-3.5
input_dir:        CNAG999_exome
```

Use cohort mode for a family, maternal-lineage analysis, or a project-level
mtDNA table. A WES/WGS cohort run is not required beforehand.

</TabItem>
</Tabs>

## Outputs

Both modes produce the same principal artifacts:

| File | Description |
| --- | --- |
| `01_mtoolbox/mit_prioritized_variants.txt` | Annotated variants with CBIcall-added `GT`, `DP`, and heteroplasmy values. |
| `01_mtoolbox/VCF_file.vcf` | mtDNA VCF from MToolBox. |
| `01_mtoolbox/mt_classification_best_results.csv` | Predicted mitochondrial haplogroups. |
| `01_mtoolbox/mit.filtered.json` | Canonical filtered records used to generate the browser. |
| `02_browser/<run-id>.html` | Standalone interactive mtDNA report. |

:::note[Haplogroup assignments]
When the available mtDNA sequence does not distinguish a single haplogroup,
MToolBox may report multiple candidate assignments. CBIcall preserves these
assignments and does not select among them.
:::

See the [mtDNA end-to-end example](../usage/end-to-end-example-mit) for the
complete run procedure and browser screenshots.

<details>
<summary>Implementation details and source files</summary>

The workflow builds on
[MToolBox v1.0](https://github.com/mitoNGS/MToolBox). It converts and prepares
the input BAM, aligns mitochondrial reads to RSRS, calls variants, predicts
haplogroups, and performs functional annotation and prioritization. CBIcall then
creates the canonical filtered JSON and standalone browser report.

- [Single-sample Bash workflow](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-3.5/mit_single.sh)
- [Cohort Bash workflow](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-3.5/mit_cohort.sh)
- Calabrese C. et al. *MToolBox: a highly automated pipeline for heteroplasmy annotation and prioritization analysis of human mitochondrial variants in high-throughput sequencing.* Bioinformatics (2014). [Article](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btu483)

</details>
