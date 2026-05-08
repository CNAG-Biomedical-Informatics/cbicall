# Performance

## Runtime Behavior

CBIcall itself uses very little memory for orchestration. The Python wrapper typically remains below **2% of a 16 GB system**. Most memory and CPU usage comes from external tools:

- **BWA-MEM**
  Memory usage increases with thread count and reference size.
  BWA does not provide an internal memory cap, so limiting RAM requires external mechanisms such as `ulimit`.

- **GATK and Picard**
  These tools default to using **8 GB** of memory.
  This value can be adjusted through the CBIcall GATK 4.6 [environment file](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-4.6/env.sh#L19) or the Snakemake [configuration file](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/snakemake/gatk-4.6/config.yaml#L32).

:::info[Cohort mode with GATK 4.6]
Joint genotyping defaults to **64 GB** of RAM for `GenomicsDBImport` and `GenotypeGVCFs`.
The value is controlled by `MEM_GENOTYPE` in the GATK 4.6 [environment file](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-4.6/env.sh#L20) and by `mem_genotype` in the Snakemake workflow.
:::

### Parallelization

Parallel execution is supported, but performance does not scale linearly with additional threads. In practice, optimal throughput is usually achieved with **4-6 threads per task**.

For example, on a 12-core workstation:

- Running **3 tasks with 4 threads each** is typically preferable to
- Running **1 task with all 12 threads**

The benchmark below shows the shape of this scaling for one WES single-sample run. The biggest gain comes from moving from 2 to 4 threads; after 6 threads, the improvement is small.

![Run time versus number of threads for WES single-mode](/img/run-time.svg)

| Threads | Runtime (minutes) |
| ---: | ---: |
| 2 | 28.5 |
| 4 | 23.4 |
| 6 | 22.1 |
| 8 | 22.0 |
| 10 | 21.9 |
| 12 | 21.4 |

:::tip[Practical default]
For batch processing, start with **4 threads per task** and scale by running more tasks in parallel when the machine or scheduler has available cores.
:::
