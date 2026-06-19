# Cross-Environment Reproducibility

Cross-environment reproducibility checks whether the same CBIcall analysis can
recover the **same final variant calls** across machines and runtime
environments.
It is separate from truth-set benchmarking, which is covered in
[GIAB Benchmarking](giab).

To evaluate CBIcall portability and reproducibility, a 1000 Genomes Project WES
sample, `HG00103` (`SRR1596639`), was processed with identical workflow
definitions across four computational environments. The WES sample used the
native Bash GATK 4.6 single-sample workflow on b37 with the
`cbicall-germline-resources-v1` resource bundle.

:::note[Installation smoke test]
CBIcall also includes a small repository reference sample used by the
installation and backend-equivalence checks. That sample is covered in
[Integration Tests](integration-tests); it is useful for confirming that an
installation works before running larger external data.
:::

## Compared Runs

| Run label | Environment | Notes |
| --- | --- | --- |
| `gcloud` | Google Cloud Platform VM, Ubuntu 22.04 | Cloud VM run. |
| `ws5` | x86_64 workstation, Linux Mint 20.3 | Local workstation baseline. |
| `ws1` | macOS ARM64 workstation running Ubuntu 24.04 in a VM | Heterogeneous CPU architecture. |
| `hpc` | SLURM-managed HPC cluster, x86_64 CentOS Linux 7.9.2009 | Production-style HPC environment. |

All runs were compared with `cbicall compare-runs` using the final QC VCF
fingerprints recorded in each `run-report.json`. The **primary evidence** is the
CBIcall text and HTML comparison report, not the optional MultiQC export.

## Result

The final VCF contained **23,562 variant records**: 19,578 `PASS` records
and 3,984 non-PASS records. All four runs produced identical final variant
calls under the CBIcall **call-level** VCF fingerprint.
The `ws5`, `hpc`, and `gcloud` runs also matched under the stricter full-record
fingerprint. The ARM `ws1` run differed only under the strict-record
fingerprint. Manual inspection localized this to a minor numeric difference in
non-call fields for one record (`3:196281208 T>C`), consistent with
environment-level numerical drift rather than a changed variant call.

![Screenshot of the 1000 Genomes cross-environment CBIcall Final VCF calls NxN heatmap, showing identical call-level fingerprints across all four environments.](/img/cross-environment-compare-outputs.png)

_Final VCF `calls` NxN heatmap from `compare-runs.html`; all 23,562 call-level records, including PASS and non-PASS records, match across every environment pair._

| Comparison target | VCF `calls` fingerprint | VCF `strict records` fingerprint | Interpretation |
| --- | --- | --- | --- |
| `ws5` baseline | `727f877de6ec...0a26a976` | `a659e0105d8b...4a5539c0` | Baseline. |
| `hpc` | `727f877de6ec...0a26a976` | `a659e0105d8b...4a5539c0` | Matches baseline at call and strict-record level. |
| `gcloud` | `727f877de6ec...0a26a976` | `a659e0105d8b...4a5539c0` | Matches baseline at call and strict-record level. |
| `ws1` | `727f877de6ec...0a26a976` | `5c03b2a73c25...7407487e` | Call-equivalent to baseline; strict-record-only numeric drift. |

<details className="cbicallInfoDetails">
<summary>Call-level versus strict VCF hashes</summary>

CBIcall records two VCF fingerprints because different questions need different
levels of sensitivity.

| Fingerprint | Fields hashed | Example of change detected | Typical interpretation |
| --- | --- | --- | --- |
| `calls` | `CHROM`, `POS`, `REF`, `ALT`, `FILTER`, and each sample `GT` in VCF sample order, across all final VCF records. | A variant becomes filtered instead of `PASS`, a genotype changes from `0/1` to `0/0`, or a site appears/disappears. | Final reported calls changed. This is the primary reproducibility check for WES/WGS outputs. Because `FILTER` is hashed, PASS and non-PASS records are both audited. |
| `strict records` | Complete sorted non-header VCF records, including `QUAL`, `INFO`, all `FORMAT` fields, `PL`, annotations, and other numeric fields. | Same `PASS` and `GT`, but `QUAL` shifts from `2267.64` to `2266.64`, or `PL` shifts from `2275,0,3196` to `2274,0,3196`. | The full VCF record changed. This can detect small numerical drift even when the final call is unchanged. |

In this validation, `ws1` and `ws5` had the same `FILTER=PASS` and `GT=0/1` for
`3:196281208 T>C`, so the **calls hash matched**. The strict hash changed because
only non-call numeric fields shifted.

</details>


## Evidence Artifacts

| Artifact | Use |
| --- | --- |
| `compare-runs.txt` | Canonical text audit report for archiving and diffing. |
| `compare-runs.html` | Static browser report; this validation uses the Final VCF `calls` heatmap plus the VCF `calls` / `strict records` rows. |
| `compare-runs_mqc/` | Optional MultiQC custom-content export for projects that already aggregate QC with MultiQC. It is not the primary evidence artifact for this validation. |

The evidence bundle can be regenerated from collected run directories:

```bash
bin/cbicall compare-runs \
  ws5/HG00103/SRR1596639/cbicall_bash_gatk-4.6_wes_single_b37_178178909641565/ \
  hpc/HG00103/SRR1596639/cbicall_bash_gatk-4.6_wes_single_b37_178179425143788/ \
  gcloud/HG00103/SRR1596639/cbicall_bash_gatk-4.6_wes_single_b37_178185765566264/ \
  ws1/HG00103/SRR1596639/cbicall_bash_gatk-4.6_wes_single_b37_178153803413337/ \
  --alias ws5 hpc gcloud ws1 \
  --output compare-runs.txt \
  --html compare-runs.html
```

:::tip[Optional MultiQC export]
Add `--multiqc compare-runs_mqc` only if the comparison needs to be included in
a larger MultiQC project report. For interpreting this reproducibility check,
use `compare-runs.txt` and `compare-runs.html`.
:::
