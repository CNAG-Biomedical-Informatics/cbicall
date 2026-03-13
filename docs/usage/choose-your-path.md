# Choose Your Path

Use this page to decide which installation route and workflow to run before diving into the detailed documentation.

---

## 1. Choose how you will run CBIcall

### Use Docker if:

- You are on a local Linux workstation or server.
- You want the fastest setup.
- You do not need an HPC scheduler integration.

Go to [Docker installation](../installation/docker.md).

### Use Apptainer if:

- You are on an HPC cluster.
- Docker is not available.
- You need a container workflow that works with batch schedulers.

Go to [Apptainer installation](../installation/apptainer.md).

### Use the non-containerized setup if:

- You want to inspect or modify the source tree directly.
- You already manage Python and Java on the host.
- You are developing or debugging workflows.

Go to [Non-containerized installation](../installation/non-containerized.md).

---

## 2. Choose the analysis you want to run

### Run `wes` if:

- Your input data is whole-exome sequencing.
- You want standard germline calling for one sample or a cohort.

See [WES/WGS single-sample pipeline](../pipelines/wes-wgs-single.md) and [WES/WGS cohort pipeline](../pipelines/wes-wgs-cohort.md).

### Run `wgs` if:

- Your input data is whole-genome sequencing.
- You want the same general calling flow as WES, but genome-wide.

See [WES/WGS single-sample pipeline](../pipelines/wes-wgs-single.md) and [WES/WGS cohort pipeline](../pipelines/wes-wgs-cohort.md).

### Run `mit` if:

- You want mtDNA calling with MToolBox.
- You already have BAM output from a previous WES/WGS single-sample run.

Important: mtDNA mode does not start from FASTQ in this project documentation flow.

See [mtDNA pipelines](../pipelines/mtdna.md) and [End-to-end example (mtDNA)](end-to-end-example-mit.md).

---

## 3. Choose `single` or `cohort`

### Use `single` if:

- You are processing one individual.
- You need a per-sample VCF or gVCF.
- You are preparing inputs for a later cohort run.

### Use `cohort` if:

- You want joint genotyping across multiple samples.
- You already have the required per-sample gVCFs for WES/WGS.
- For mtDNA, you already have the required WES/WGS BAM outputs for all samples.

---

## 4. Choose the workflow engine

### Use `bash` if:

- You want the most broadly supported option.
- You are running mtDNA.
- You want behavior that matches the legacy and example scripts most closely.

### Use `snakemake` if:

- You want Snakemake orchestration for supported WES/WGS workflows.
- You are using `gatk-4.6`.

Important: `snakemake` is not supported for `mit`, and it is not supported with `gatk-3.5`.

---

## Recommended Starting Points

- If you want the shortest validated smoke test, go to [Quickstart](quickstart.md).
- If you want a realistic WES example from inputs to outputs, go to [End-to-end example (WES)](end-to-end-example-wes.md).
- If you want a realistic mtDNA example, go to [End-to-end example (mtDNA)](end-to-end-example-mit.md).
