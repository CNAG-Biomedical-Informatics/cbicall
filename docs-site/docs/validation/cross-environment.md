# Cross-Environment Reproducibility

Cross-environment reproducibility checks whether the same CBIcall analysis can
recover the same final variant calls across machines and runtime environments.
It is separate from truth-set benchmarking, which is covered in
[GIAB Benchmarking](giab).

This example uses 1000 Genomes WES sample `HG00103` (`SRR1596639`) with the
native Bash GATK 4.6 WES single-sample workflow on b37 and the
`cbicall-germline-resources-v1` resource bundle.

## Compared Runs

| Run label | Environment role | Notes |
| --- | --- | --- |
| `ws5` | Local x86_64 workstation | Baseline run. |
| `hpc` | SLURM/CentOS HPC | Production-style HPC environment. |
| `ws1` | ARM workstation | Heterogeneous CPU architecture. |
| `gcloud` | Google Cloud VM | Cloud VM run included as an outlier example. |

All runs were compared with `cbicall compare-runs` using the final QC VCF
fingerprints recorded in each `run-report.json`.

## Result

| Comparison target | VCF `calls` fingerprint | VCF `strict records` fingerprint | Interpretation |
| --- | --- | --- | --- |
| `ws5` baseline | `727f877de6ec...0a26a976` | `a659e0105d8b...4a5539c0` | Baseline. |
| `hpc` | `727f877de6ec...0a26a976` | `a659e0105d8b...4a5539c0` | Matches baseline at call and strict-record level. |
| `ws1` | `727f877de6ec...0a26a976` | `5c03b2a73c25...7407487e` | Call-equivalent to baseline; strict-record-only numeric drift. |
| `gcloud` | `d962213762e2...46e18fc7` | `53ca560f371e...6a6455c3` | Call-level outlier in this comparison. |

The `ws5`, `hpc`, and ARM `ws1` runs produced identical final variant calls
under the CBIcall call-level VCF fingerprint. The ARM run differed only under
the stricter full-record fingerprint. Manual inspection localized this to a
minor numeric difference in non-call fields for one record
(`3:196281208 T>C`, with unchanged filter and genotype), consistent with
environment-level numerical drift rather than a changed variant call.

The `gcloud` run differed at the call level in this comparison. Treat this as a
variant-output discrepancy and inspect execution conditions before grouping it
with the call-equivalent runs.

![Screenshot of the 1000 Genomes cross-environment CBIcall compare-runs HTML Outputs matrix.](/img/cross-environment-compare-outputs.png)

![Screenshot of the 1000 Genomes cross-environment MultiQC pairwise summary, restricted to output and final-VCF columns.](/img/cross-environment-multiqc-final-vcf.png)

## Evidence Artifacts

| Artifact | Use |
| --- | --- |
| `compare-runs.txt` | Canonical text audit report for archiving and diffing. |
| `compare-runs.html` | Static browser report; this validation uses the `Outputs` matrix and VCF `calls` / `strict records` rows. |
| `compare-runs_mqc/` | MultiQC custom-content bundle; this validation uses the output and Final VCF `calls` / `strict records` columns. |

The evidence bundle can be regenerated from collected run directories:

```bash
bin/cbicall compare-runs \
  ws5/HG00103/SRR1596639/cbicall_bash_gatk-4.6_wes_single_b37_178178909641565/ \
  gcloud/HG00103/SRR1596639/cbicall_bash_gatk-4.6_wes_single_b37_178153768026219/ \
  hpc/HG00103/SRR1596639/cbicall_bash_gatk-4.6_wes_single_b37_178179425143788/ \
  ws1/HG00103/SRR1596639/cbicall_bash_gatk-4.6_wes_single_b37_178153803413337/ \
  --alias ws5 gcloud hpc ws1 \
  --output compare-runs.txt \
  --html compare-runs.html \
  --multiqc compare-runs_mqc
```
