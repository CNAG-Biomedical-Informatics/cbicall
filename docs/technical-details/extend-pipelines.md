CBIcall is designed to be extensible, allowing developers to integrate new analysis pipelines with minimal changes to the core framework.  
This guide explains the recommended structure, required components, and best practices for adding a custom pipeline.

---

## 1. Overview

A CBIcall pipeline (or workflow) consists of:

- One workflow scripts (Bash or Snakemake) under `workflows/{bash,snakemake}/gatk_{3.5,4.6}/`
- Registration inside the Python wrapper
- Configuration options exposed through the YAML parameters

Following the existing structure ensures that your pipeline integrates cleanly with the logging, directory layout, QC modules, and project management mechanisms.

---

## 2. Create the Pipeline

Within the `workflows/` folder, you can choose one or both workflow engines:

### **Bash workflow**

```
workflows/bash/gatk_4.6/
  mypipeline_single.sh
```

### **Snakemake workflow**

```
workflows/snakemake/gatk_4.6/
    mypipeline_single.smk
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

If your pipeline requires extra options (e.g., target BED, annotation files, specialized tools), add them to `parameters.sh`.

---

## 5. Follow CBIcall Directory Conventions

CBIcall generates a standardized output structure:

```
cbicall_run_<pipeline_name>/
  logs/
  01_step/
  02_step/
  . . .
```

---

## 6. Support Single and/or Cohort Mode

If your pipeline supports both:

- `mypipeline_single.sh` handles standalone samples
- `mypipeline_cohort.sh` merges or jointly processes samples

If not applicable, document the limitation:

```yaml
mode: single   # cohort not supported
```

CBIcall will raise an informative error if a user selects an unsupported mode.

---

## 7. Add Test Data (Recommended)

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

## 8. Document Your Pipeline

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
