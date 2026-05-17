# Run Comparison

`cbicall compare-runs` compares completed CBIcall run directories or `run-report.json` files. It is intended for reproducibility checks across repeated local runs, HPC runs, container runs, or cloud runs.

```bash
bin/cbicall compare-runs run_a/ run_b/ run_c/ \
  --output compare-report.txt \
  --html compare-report.html
```

With two runs, CBIcall prints a direct pairwise comparison. With three or more runs, the first run is used as the baseline and the remaining runs are compared against it.

## What Is Compared

| Layer | Fields |
| --- | --- |
| Framework | CBIcall version. |
| Pipeline | Workflow key, pipeline implementation version, entrypoint, workflow fingerprint. |
| Workflow files | Entrypoint and helper/config file paths and SHA-256 values. |
| Resources | Resource key and resource fingerprint. |
| Outputs | Normalized VCF fingerprints when `03_stats/*.vcf.sha256.txt` is present. |

The status vocabulary is intentionally small:

| Status | Meaning |
| --- | --- |
| `same` | Values or fingerprints match. |
| `different` | Values or fingerprints exist in all compared runs but differ. |
| `missing` | Evidence is present in only some runs. |
| `not available` | Evidence is not recorded in any compared run. |

## Reports

The text report is the canonical audit artifact because it is easy to diff, archive, and attach to review material. The HTML report is a static rendering of the same information for browsing.

![Screenshot of the CBIcall HTML run comparison report showing summary counters, baseline run metadata, and status pills for framework, pipeline, and workflow file checks.](/img/compare-runs-html-report.png)

## Interpretation

A changed workflow fingerprint means the resolved workflow files are not byte-identical. That is audit evidence, not automatically a failed analysis. For example, editing a comment in a Bash workflow changes the workflow hash but may leave the normalized VCF fingerprint unchanged.

For output reproducibility, prioritize the normalized VCF fingerprint. For implementation provenance, inspect the workflow and helper file fingerprints.
