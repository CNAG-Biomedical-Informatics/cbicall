import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Performance

## Runtime Behavior

CBIcall adds negligible orchestration overhead. The Python wrapper typically remains below **2% of a 16 GB system**, does not process reads or variants, and does not create Python worker threads. The `-t/--threads` value is passed to the selected workflow backend.

Most memory and CPU usage comes from external tools:

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

:::note[Python driver overhead]
The Python driver is expected to require one CPU core only during short setup phases. For long-running variant-calling jobs, scheduler CPU and memory requests should be sized for the selected external tools and the requested workflow threads, not for the CBIcall Python process itself.
:::

### Parallelization

<Tabs groupId="performance-benchmark">
<TabItem value="wes" label="WES" default>

Parallel execution is supported, but performance does not scale linearly with additional threads. In practice, optimal throughput is usually achieved with **4-6 threads per task**.

For example, on a 12-core workstation:

- Running **3 tasks with 4 threads each** is typically preferable to
- Running **1 task with all 12 threads**

The benchmark below shows the shape of this scaling for **WES single-sample calling** on the 1000 Genomes WES sample **SRR1596639**. The paired FASTQ inputs were **1.9 GB** (`R1`) and **2.0 GB** (`R2`). Runs used an **HP Z2 G8 Tower Workstation** (`x86_64`) with an **Intel Xeon W-1350P @ 4.00 GHz**, **6 physical cores / 12 hardware threads**, and **31 GiB RAM**. The GATK/Picard memory setting was fixed at `MEM=8G` in `env.sh` for all thread counts.

The biggest gain comes from moving from 2 to 4 threads; after 6 threads, the improvement is small.

![Run time versus number of threads for WES single-mode](/img/run-time.svg)

| Threads | Elapsed seconds | Runtime (minutes) |
| ---: | ---: | ---: |
| 2 | 7856.302 | 130.9 |
| 4 | 6342.747 | 105.7 |
| 6 | 5812.155 | 96.9 |
| 8 | 5695.981 | 94.9 |
| 10 | 5525.616 | 92.1 |
| 12 | 5452.732 | 90.9 |

</TabItem>
<TabItem value="wgs" label="WGS">

WGS cluster benchmark data will be added here.

</TabItem>
</Tabs>

:::tip[Practical default]
For batch processing, start with **4 threads per task** and scale by running more tasks in parallel when the machine or scheduler has available cores.
:::
