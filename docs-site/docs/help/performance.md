import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Performance

## Runtime Behavior

Throughout this page, **GiB** is used for file sizes reported by Linux tools such as `ls -lh` and for Java heap settings such as `-Xmx8G`. **GB** is reserved for decimal storage specifications when the source explicitly reports decimal gigabytes.

The `-t/--threads` value is passed to the selected workflow backend. Most memory and CPU usage comes from external tools:

- **BWA-MEM**
  Memory usage increases with thread count and reference size.
  BWA does not provide an internal memory cap, so limiting RAM requires external mechanisms such as `ulimit`.

- **GATK and Picard**
  These tools default to using **8 GiB** of memory (`-Xmx8G`).
  This value can be adjusted through the CBIcall GATK 4.6 [environment file](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-4.6/env.sh#L19) or the Snakemake [configuration file](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/snakemake/gatk-4.6/config.yaml#L32).

:::info[Cohort mode with GATK 4.6]
Joint genotyping defaults to **64 GiB** of RAM for `GenomicsDBImport` and `GenotypeGVCFs`.
The value is controlled by `MEM_GENOTYPE` in the GATK 4.6 [environment file](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/workflows/bash/gatk-4.6/env.sh#L20) and by `mem_genotype` in the Snakemake workflow.
:::

:::note[Python driver overhead]
CBIcall adds negligible orchestration overhead. The Python wrapper typically remains below **2% of a 16 GiB system**, does not process reads or variants, and does not create Python worker threads. It is expected to require one CPU core only during short setup phases. For long-running variant-calling jobs, scheduler CPU and memory requests should be sized for the selected external tools and workflow threads, not for the CBIcall Python process itself.
:::

### Parallelization

:::note[Total CPU time]
Some workflow steps can be split further, for example across FASTQ chunks or genomic intervals. This may shorten one job, but it does not necessarily reduce total CPU time. On a Slurm cluster with a fixed CPU allocation, **1 job with 24 threads** and **6 jobs with 4 threads each** use the same simultaneous CPU budget. CBIcall deliberately favors moderate per-job thread counts because the intended production use case is running thousands of jobs on Slurm. This matters most for WGS, where one job can run for days.
:::

### Benchmark Setup

The WES and WGS benchmarks below were run on the same **HP Z2 G8 Tower Workstation** running **Linux Mint 20.3 (Una)** (`x86_64`) with an **Intel Xeon W-1350P @ 4.00 GHz**, **6 physical cores / 12 hardware threads**, and **32 GB RAM**. Runs wrote to an **HDD** rather than an SSD. The GATK/Picard Java heap setting was fixed at `MEM=8G` in `env.sh` for all thread counts, corresponding to an 8 GiB heap request.

<Tabs groupId="performance-benchmark">
<TabItem value="wes" label="WES" default>

Parallel execution is supported, but performance does not scale linearly with additional threads. In practice, optimal throughput is usually achieved with **4-6 threads per task**.

For example, on a 12-core workstation:

- Running **3 tasks with 4 threads each** is typically preferable to
- Running **1 task with all 12 threads**

The benchmark below shows the shape of this scaling for **WES single-sample calling** on the 1000 Genomes sample **HG00103** using run accession **SRR1596639**. The paired FASTQ inputs were **1.9 GiB** (`R1`) and **2.0 GiB** (`R2`) by `ls -lh`, with mean coverage of **97.2x**.

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

WGS single-sample runs are substantially longer than WES runs. Measurements for the 1000 Genomes sample **HG00097** used paired FASTQ inputs totaling approximately **155 GiB** compressed by `ls -lh`, with mean coverage of **31.9x**, and show modest improvement when increasing from 8 to 12 threads.

![Run time versus number of threads for WGS single-mode](/img/run-time-wgs.svg)

| Threads | Runtime (hours) |
| ---: | ---: |
| 4 | 70.0 |
| 8 | 63.45 |
| 12 | 60.0 |

</TabItem>
</Tabs>

:::tip[Practical default]
For batch processing, start with **4 threads per task** and scale by running more tasks in parallel when the machine or scheduler has available cores.
:::
