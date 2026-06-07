# Cross-Environment Reproducibility

Use this page to document CBIcall reproducibility across machines, operating
systems, cloud instances, and HPC environments. The goal is to audit the
framework layer: the same parameters, workflow registry entry, resource bundle,
and final-output fingerprint should be recovered across environments.

## Planned Matrix

A compact validation matrix can include two WES datasets across four execution
environments.

| Dataset | Local workstation | ARM workstation | Google Cloud | SLURM/CentOS HPC |
| --- | --- | --- | --- | --- |
| Bundled WES test dataset | pending | pending | pending | pending |
| Public 1000G WES sample | pending | pending | pending | pending |

For each completed run, archive:

- `run-report.json`
- `run-report.html`
- `log.json`
- `cbicall-execution-contract.json`
- the workflow log
- the `compare-runs` text and HTML reports

## What To Compare

Run comparison should report at least:

| Layer | Expected interpretation |
| --- | --- |
| CBIcall version | Same for strict release reproducibility. |
| Python, Java, backend version | Useful runtime evidence; may differ when testing portability. |
| Workflow fingerprint | Same when the same workflow files were executed. |
| Resource fingerprint | Same when the same declared resource bundle was used. |
| File inventory | May differ because logs and backend metadata can differ across environments. |
| Normalized final VCF hash | Primary final-output reproducibility check. |

Example command:

```bash
bin/cbicall compare-runs run_local/ run_cloud/ run_hpc/ \
  --alias local cloud hpc \
  --output cross-env-compare.txt
```

For three or more runs, the report includes both the baseline comparison and an
all-to-all reproducibility matrix.

:::note[Interpretation]
A different file inventory does not necessarily mean a different biological
result. For CBIcall WES/WGS reproducibility claims, report whether the normalized
final VCF fingerprint is the same and then describe any differences in runtime
metadata, logs, or auxiliary files.
:::
