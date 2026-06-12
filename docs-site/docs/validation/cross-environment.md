# Cross-Environment Reproducibility

Use this page to document CBIcall reproducibility across machines, operating
systems, cloud instances, and HPC environments. The goal is to audit the
framework layer: the same parameters, workflow registry entry, resource bundle,
and final-output fingerprint should be recovered across environments.

:::important[Scope]
This is a **framework reproducibility** check. It asks whether CBIcall can launch
the same validated workflow contract across environments and recover comparable
audit evidence. It is separate from analytical benchmarking against truth sets,
which is covered in [GIAB Benchmarking](giab).
:::

## Validation Matrix

A compact validation matrix can use two WES datasets across four execution
environments.

| Dataset | Local workstation | ARM workstation | Google Cloud | SLURM/CentOS HPC |
| --- | --- | --- | --- | --- |
| Bundled WES test dataset | Run native WES test | Run native WES test | Run native WES test | Run native WES test |
| Public 1000G WES sample | Run same WES YAML | Run same WES YAML | Run same WES YAML | Run same WES YAML |

The bundled test dataset is useful because the expected normalized VCF hash is
known. The public 1000G WES sample is useful because it exercises a less trivial
real-data run outside the tiny integration fixture.

## Run The Same Contract

Use the same parameters YAML and resource bundle in each environment. For the
bundled test dataset:

```bash
bin/cbicall test --wes-bash -t 1
```

For a real WES sample, run the same YAML everywhere:

```bash
bin/cbicall run -p cloud-wes-real.yaml -t 4
```

On HPC, use the same runtime profile used for production runs:

```bash
bin/cbicall run -p cloud-wes-real.yaml -t 4 --runtime-profile cnag-hpc
```

:::tip[Keep the run directory]
Each command prints the generated `cbicall_*` run directory. Keep that directory
or at least its `run-report.json`; `compare-runs` can use either.
:::

## Compare Environments

After all runs finish, compare the completed run directories from one machine
that can access the reports:

```bash
bin/cbicall compare-runs run_local/ run_arm/ run_cloud/ run_hpc/ \
  --alias local arm cloud hpc \
  --output cross-env-compare.txt
```

For three or more runs, CBIcall automatically writes both baseline and all-to-all
evidence. The HTML report is written next to the text report by default:

```text
cross-env-compare.txt
cross-env-compare.html
```

For the detailed comparison UI and screenshots, see [Run Comparison](run-comparison).

## What To Check

| Layer | Expected interpretation |
| --- | --- |
| CBIcall version | Should match for strict release reproducibility. |
| Python, Java, backend version | Runtime evidence; may differ when testing portability. |
| Workflow fingerprint | Should match when the same workflow files were executed. |
| Resource fingerprint | Should match when the same declared resource bundle was used. |
| Execution contract fingerprint | Should match after normalization when the same backend-ready plan was launched. |
| File inventory | May differ because logs and backend metadata can differ across environments. |
| Normalized final VCF hash | Primary final-output reproducibility check. |

:::note[Interpretation]
A different file inventory does not necessarily mean a different biological
result. For CBIcall WES/WGS reproducibility claims, report whether the normalized
final VCF fingerprint is the same and then describe any differences in runtime
metadata, logs, or auxiliary files.
:::

## Minimal Evidence Table

For a manuscript or reviewer response, summarize the comparison in a compact
matrix rather than listing every file path.

| Dataset | Environments compared | Workflow/resource evidence | Final VCF evidence |
| --- | --- | --- | --- |
| Bundled WES test dataset | local, ARM, Google Cloud, SLURM/CentOS | same workflow and resource fingerprints | same normalized final VCF hash |
| Public 1000G WES sample | local, ARM, Google Cloud, SLURM/CentOS | same workflow and resource fingerprints | same normalized final VCF hash |

Replace the cells with the actual `compare-runs` result labels from your runs.

<details>
<summary>Files to keep for audit</summary>

| File | Why it matters |
| --- | --- |
| `cross-env-compare.txt` | Canonical text comparison report. |
| `cross-env-compare.html` | Browser view of baseline and all-to-all evidence. |
| `run-report.json` | Compact per-run audit report used by `compare-runs`. |
| `log.json` | Resolved parameters, workflow, resource, and runtime evidence. |
| `cbicall-execution-contract.json` | Backend-ready launch plan and normalized execution fingerprint. |
| Workflow log | Backend-specific execution log. |

</details>
