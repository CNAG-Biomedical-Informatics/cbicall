# Adding Your Own Pipeline to CBIcall

CBICall is designed to be extensible, allowing developers to integrate new analysis pipelines with minimal changes to the core framework.  
This guide explains the recommended structure, required components, and best practices for adding a custom pipeline.

---

## 1. Overview

A CBICall pipeline consists of:

- A directory under `pipelines/<pipeline_name>/`
- One or more workflow scripts (Bash or Snakemake)
- Registration inside the Python wrapper
- Configuration options exposed through the YAML parameters
- (Optional) container definitions and test data

Following the existing structure ensures that your pipeline integrates cleanly with the logging, directory layout, QC modules, and project management mechanisms.

---

## 2. Create the Pipeline Directory

Inside the `pipelines/` folder, create a directory for your pipeline:

```
pipelines/
  mypipeline/
```

Within this folder, you can choose one or both workflow engines:

### **Bash workflow**

```
pipelines/mypipeline/
  mypipeline_single.sh
  mypipeline_cohort.sh
```

### **Snakemake workflow**

```
pipelines/mypipeline/
  Snakefile
  rules/
    step1.smk
    step2.smk
```

The naming convention follows the existing WES/WGS/mit modules.

---

## 3. Register the Pipeline in the Python Wrapper

Edit the main driver script (e.g., `cbicall.py`) so that it recognizes the new pipeline:

```python
if cfg.pipeline == "mypipeline":
    run_mypipeline(cfg)
```

If your pipeline supports both engines:

```python
if cfg.workflow_engine == "bash":
    run_mypipeline_bash(cfg)

elif cfg.workflow_engine == "snakemake":
    run_mypipeline_snakemake(cfg)
```

A minimal wrapper function typically prepares paths, validates inputs, and dispatches the workflow script.

---

## 4. Define YAML Configuration Options

Every pipeline requires a documented set of parameters.  
At minimum, support:

```yaml
pipeline: mypipeline
mode: single
reference: /path/to/ref
workflow_engine: bash
sample: /path/to/sample_prefix
projectdir: /path/to/output
```

If your pipeline requires extra options (e.g., target BED, annotation files, specialized tools), add them under:

```yaml
mypipeline_options:
    param1: value
    param2: value
```

Keep the structure consistent with the existing pipelines.

---

## 5. Follow CBICall Directory Conventions

CBICall generates a standardized output structure:

```
cbicall_run/
  logs/
  work/
  results/
    <pipeline_name>/
```

To ensure compatibility:

- Write temporary or intermediate files into `work/`
- Place final outputs under `results/<pipeline_name>/`
- Do not overwrite data from other pipelines
- Use logging consistent with the rest of CBICall (stdout + file)

---

## 6. Support Single and/or Cohort Mode

If your pipeline supports both:

- `mypipeline_single.sh` handles standalone samples
- `mypipeline_cohort.sh` merges or jointly processes samples

If not applicable, document the limitation:

```yaml
mode: single   # cohort not supported
```

CBICall will raise an informative error if a user selects an unsupported mode.

---

## 7. Optional: Snakemake Integration

If using Snakemake:

- Include a `Snakefile`
- Place reusable rules in a `rules/` subdirectory
- Use CBICall's existing pattern for cluster submission, logging, and checkpoints
- Expose thread/memory settings using Snakemake resources

Example rule snippet:

```python
rule align:
    input:
        r1 = "{sample}_R1.fastq.gz",
        r2 = "{sample}_R2.fastq.gz"
    output:
        bam = "work/{sample}.bam"
    threads: 8
    resources:
        mem_mb = 16000
    shell:
        "bwa mem -t {threads} {input.r1} {input.r2} | samtools sort -o {output.bam}"
```

---

## 8. Add Test Data (Recommended)

Create a minimal test dataset:

```
examples/mypipeline_test/
  input.fastq.gz
  test.yaml
```

The test should:

- run quickly (<1 minute ideally)
- validate key steps
- generate deterministic output

This helps users and CI workflows verify that the pipeline works.

---

## 9. Document Your Pipeline

Add a new page under:

```
docs/pipelines/mypipeline.md
```

Document:

- Requirements
- YAML parameters
- Supported modes (single/cohort)
- Example commands
- Output structure

Link the page in `mkdocs.yml` under Technical Details.

---

## 10. (Optional) Add Container Support

If your pipeline requires additional software:

- Extend the existing Dockerfile, or
- Create a new image under `docker/`

Ensure reproducibility by pinning versions of required tools.

---

## Summary

To add a new pipeline:

1. Create a pipeline directory  
2. Add Bash or Snakemake workflow scripts  
3. Register the pipeline in the Python wrapper  
4. Define YAML configuration options  
5. Follow output + logging conventions  
6. Support single/cohort modes if needed  
7. Provide optional test data + docs  

Following this structure ensures your pipeline integrates cleanly with the CBICall framework and remains compatible as the project evolves.

