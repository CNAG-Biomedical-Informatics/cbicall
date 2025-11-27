### Runtime behavior

CBICall itself uses very little memory for orchestration. The Python wrapper typically remains below **2% of a 16 GB system**. Most memory and CPU usage comes from external tools:

- **BWA-MEM**
  Memory usage increases with thread count and reference size.
  BWA does not provide an internal memory cap, so limiting RAM requires external mechanisms such as `ulimit`.

- **GATK and Picard**
  These tools default to using **8 GB** of memory.
  This value can be adjusted through the CBICall configuration file.

### Parallelization

Parallel execution is supported, but performance does not scale linearly with additional threads. In practice, optimal throughput is usually achieved with **about 4 threads per task**.

For example, on a 12-core workstation:

- Running **3 tasks with 4 threads each** is typically faster than
- Running a **single task with all 12 threads**

A runtime example is shown in the figure below.

!!! Example "Performance"

    ![Time](https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/docs/img/run-time.png)
