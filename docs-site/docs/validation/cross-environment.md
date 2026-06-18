# Cross-Environment Reproducibility

Cross-environment reproducibility checks whether CBIcall can run the same
validated analysis contract on different machines and recover the same
workflow/resource identity and canonical output fingerprints. It is separate
from truth-set benchmarking, which is covered in [GIAB Benchmarking](giab).

The reviewer-facing check uses one real WES sample. The bundled WES test can be
used first as an optional smoke test, but it should not be presented as the main
cross-environment reproducibility result.

| Check | Dataset | Purpose |
| --- | --- | --- |
| Optional smoke test | Bundled WES test data | Confirms installation, resource resolution, workflow dispatch, and report generation before running the real sample. |
| Reproducibility evidence | One public WES sample | Tests whether the same workflow/resource contract and normalized final VCF fingerprint are recovered across environments. |

## Environments

Run the checks in each environment being claimed. For the manuscript response,
the intended matrix is:

| Environment | Runtime style | What to run |
| --- | --- | --- |
| Local workstation | Native Bash or Docker-supported local run | One real WES sample |
| Google Cloud VM | Fresh source install on Docker-capable VM | Same real WES sample |
| SLURM/CentOS HPC | Native Bash with site runtime profile | Same real WES sample |

The Google Cloud setup recipe is documented separately in
[Google Cloud](../installation/google-cloud-docker).

## Optional Smoke Test

Use the bundled test only to check that a new environment is ready before
spending time on real data. This is an integration/smoke test, not the main
cross-environment reproducibility result.

```bash
bin/cbicall test --wes-bash -t 1
```

Keep the generated `cbicall_*` run directory if you want to compare smoke-test
runs later. If storage is tight, keep at least:

```text
run-report.json
log.json
cbicall-execution-contract.json
```

## 1. Run One Real WES Sample

Run one real WES sample with the same parameters and resource bundle in each
environment. This is the practical check that the reproducibility machinery also
works on realistic input.

Example parameter file:

```yaml
mode: single
pipeline: wes
workflow_backend: bash
software_stack: gatk-4.6
genome: hg38
input_dir: /path/to/public_wes_sample_fastqs
cleanup_bam: true
```

Run it locally or on the cloud VM:

```bash
bin/cbicall run -p real-wes-single.yaml -t 4
```

On HPC, use the site runtime profile if that is how production runs are launched:

```bash
bin/cbicall run -p real-wes-single.yaml -t 4 --runtime-profile cnag-hpc
```

Use the same FASTQs, same CBIcall version, same workflow backend, and same
registered resource bundle across environments whenever possible.

## 2. Compare Completed Runs

After the runs finish, collect the run directories or `run-report.json` files in
one place and compare them.

Real WES comparison:

```bash
bin/cbicall compare-runs \
  local_wes_run/ cloud_wes_run/ hpc_wes_run/ \
  --alias local-wes cloud-wes hpc-wes \
  --output cross-env-real-wes-compare.txt
```

CBIcall also writes an HTML report next to the text report:

```text
cross-env-real-wes-compare.html
```

## 3. Interpret The Report

For the reviewer response, focus on the layers that define reproducible CBIcall
execution and final biological output.

| Report layer | Expected result |
| --- | --- |
| CBIcall version | Same for strict reproducibility claims. |
| Workflow key/version/hash | Same when the same registered workflow was executed. |
| Resource key/version/hash | Same when the same declared resource bundle was used. |
| Execution contract | Same after normalization when the same backend-ready plan was launched. |
| Normalized final VCF fingerprint | Same final-output fingerprint for reproducible variant output. |
| Output inventory | May differ because logs, timestamps, backend metadata, and scheduler files can be environment-specific. |
| Hostname, paths, thread count | Runtime context; useful audit evidence but not by itself a failed reproducibility check. |

For the real WES sample, the normalized final VCF fingerprint is the main output
claim. Differences in logs or auxiliary execution files should be described as
environment metadata unless they change the canonical variant output.

## Reporting Template

Use a compact table in the manuscript or response letter:

| Check | Environments compared | Workflow/resource result | Final-output result |
| --- | --- | --- | --- |
| One public WES sample | local, Google Cloud, SLURM/CentOS | pending | pending |

Replace `pending` with the actual `compare-runs` result labels after the runs
finish.

<details className="cbicallCodeDetails">
<summary>Evidence files to keep</summary>
<div>

| File | Why it matters |
| --- | --- |
| `cross-env-real-wes-compare.txt` and `.html` | Cross-environment comparison for the real WES sample. |
| `run-report.json` | Compact per-run audit report used by `compare-runs`. |
| `log.json` | Resolved parameters, workflow, resource, and runtime evidence. |
| `cbicall-execution-contract.json` | Backend-ready launch plan and normalized execution fingerprint. |
| Workflow log | Backend-specific execution evidence. |

</div>
</details>
