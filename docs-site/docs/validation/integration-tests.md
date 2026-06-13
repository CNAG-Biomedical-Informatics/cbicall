# Integration Tests

CBIcall ships small integration tests that run example workflows and validate
the resulting run directory against small contract fixtures.

They do not replace biological or clinical validation of the underlying variant-calling methods.

<div className="cbicallNotePanel">
  <p><strong>Use this page when the question is:</strong> can this CBIcall installation run the shipped example workflows and reproduce the expected example outputs?</p>
</div>

This page is about execution tests. If the question is whether a parameters YAML
resolves before launch, use `validate-parameters` as described in
[Configuration Reference](../help/configuration-reference). If the question is
whether a custom resource catalog or installed resource directory is valid, use
[Resource Validation](../usage/resource-validation).

## Minimal Test

From the repository root:

```bash
bin/cbicall test --wes-bash -t 1
```

On HPC, pass the same runtime profile you use for normal Bash runs:

```bash
bin/cbicall test --wes-bash -t 1 --runtime-profile cnag-hpc
```

:::note[Bash-specific profile]
`--runtime-profile` selects a registry-resolved Bash `env.sh` mapping. Other
backends use their own configuration/profile mechanisms.
:::

| Check | What it confirms |
| --- | --- |
| `test --wes-bash` | The required bundled Bash WES workflow runs and reproduces the expected normalized VCF hash declared by the contract fixture. |
| `test --backend-equivalence` | Runs Bash plus available native WES backends and checks that their normalized final VCF content matches Bash. |

## Backend Equivalence Check

Use backend-equivalence mode before tagging or publishing an image:

```bash
bin/cbicall test --backend-equivalence -t 1 --runtime-profile local
```

On HPC, use the same Bash runtime profile as normal runs:

```bash
bin/cbicall test --backend-equivalence -t 1 --runtime-profile cnag-hpc
```

`--backend-equivalence` uses Bash as the required native baseline and compares
available Snakemake, Nextflow, and Cromwell WES runs against it. It compares the
**normalized final VCF hash**, not logs, workflow files, execution contracts, or
full output inventories. Missing optional backends are skipped, but at least one
non-Bash backend must be available for the backend-equivalence check to pass.

The final summary is the key audit line:

```text
Backend equivalence summary
========================================
Baseline
  WES Bash      => passed | 621373dda8e7...72559baf | 6 records (...)

Backend equivalence
  WES Snakemake => same final VCF | 621373dda8e7...72559baf | 6 records (...)
  WES Nextflow  => same final VCF | 621373dda8e7...72559baf | 6 records (...)
  WES Cromwell  => same final VCF | 621373dda8e7...72559baf | 6 records (...)

Compared non-Bash backends: 3
Status: PASSED
Exit code: 0
========================================
```

:::tip[Where to go next]
Use this page to run the shipped examples. For parameter-file checks, see
[Configuration Reference](../help/configuration-reference). For comparing
repeated runs, see [Run Comparison](run-comparison).
:::

## Test Matrix

The minimal WES Bash test is the default smoke test. Use the matrix below only
when checking optional backends or external nf-core examples.

<details>
<summary>Show all integration-test commands</summary>

| Command | Workflow path | Needs CBIcall bundle? | Extra requirement | Contract check |
| --- | --- | --- | --- | --- |
| `bin/cbicall test --wes-bash -t 1` | **Native WES**, Bash | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | <span className="cbicallTestBadge cbicallTestBadgeNeutral">none</span> | Run report fields, expected files, **normalized VCF hash** |
| `bin/cbicall test --wes-snakemake -t 1` | **Native WES**, Snakemake | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | `snakemake` on `PATH` | Run report fields, expected files, **normalized VCF hash** |
| `bin/cbicall test --wes-nextflow -t 1` | **Native WES**, Nextflow | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | `nextflow` on `PATH` | Run report fields, expected files, **normalized VCF hash** |
| `bin/cbicall test --wes-cromwell -t 1` | **Native WES**, Cromwell | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | `CROMWELL_JAR` or `cromwell` on `PATH` | Generated inputs/options/metadata, run report fields, expected files, **normalized VCF hash** |
| `bin/cbicall test --mit-bash -t 1` | **Native mtDNA**, Bash | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | x86_64 host | Run report fields, expected files, **prioritized variants hash**, **raw JSON hash** |
| `bin/cbicall test --nf-core-demo -t 4` | **nf-core/demo** | <span className="cbicallTestBadge cbicallTestBadgeNo">X bundle</span> | Nextflow plus selected nf-core runtime profile | Generated params/config, run reports, **pipeline info**, **MultiQC anchors** |
| `bin/cbicall test --nf-core-sarek -t 4` | **nf-core/Sarek** | <span className="cbicallTestBadge cbicallTestBadgeNo">X bundle</span> | Nextflow plus selected nf-core runtime profile and Sarek inputs/resources | Generated params/config, run reports, **pipeline info**, **MultiQC anchors**, declared canonical outputs when produced |
| `bin/cbicall test --backend-equivalence -t 1` | **Native WES backend equivalence** | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | At least one non-Bash native backend available | Bash baseline plus available native WES backends; **same normalized final VCF** required |
| `bin/cbicall test --all -t 1` | **Native tests only** | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | Optional backends skipped if missing | Runs WES Bash, WES Snakemake, WES Nextflow, WES Cromwell, and mtDNA contracts |

</details>

<details>
<summary>Advanced test flags and backend requirements</summary>

| Item | Use |
| --- | --- |
| `--runtime-profile cnag-hpc` | Test the same native environment profile used for production runs. |
| `--keep-external-work` | Keep nf-core/Nextflow `work/` and `.nextflow/` directories for debugging. |
| Snakemake tests | Require `snakemake` on `PATH`. |
| Nextflow and nf-core tests | Require `nextflow` on `PATH`; on CNAG HPC load it with `module load Nextflow/25.10.2`. |
| Cromwell tests | Require `CROMWELL_JAR=/path/to/cromwell.jar` or a `cromwell` launcher on `PATH`. |

</details>

## Outputs

The test command prints the run directory, workflow log, `run-report.json`,
`run-report.html`, launcher log, contract fixture, validation status, and output
hashes when the contract declares them.

:::tip[VCF hashes during tests]
WES tests compute normalized [SHA-256](https://en.wikipedia.org/wiki/SHA-2)
values directly from the newly generated VCF after removing headers and sorting
variant records. The expected value is stored in a small contract fixture, not
in a copied `ref_*` run directory. The test does not rely on stored
`03_stats/*.vcf.sha256.txt` files, because those sidecar files can be absent or
stale; they are kept for audit, run reports, and `compare-runs`.
:::

:::info[nf-core tests]
The nf-core integration tests validate CBIcall's nf-core execution envelope:
parameter resolution, generated Nextflow params/config, workflow logs, run
reports, and stable nf-core output anchors such as `pipeline_info` and MultiQC.
They do not replace the upstream nf-core `nf-test` suite.
:::

For WES, the run directory looks like:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_gatk-4.6_wes_single_b37_*/
  01_bam/
  02_varcall/
  03_stats/
  logs/
```

The optional Snakemake WES test uses the same input and expected VCF records, but
the run directory starts with `cbicall_snakemake_gatk-4.6_wes_single_b37_*`.

The optional Nextflow WES test uses the same input and expected VCF records, but
the run directory starts with `cbicall_nextflow_gatk-4.6_wes_single_b37_*`.

The optional Cromwell WES test uses the same input and expected VCF records, but
the run directory starts with `cbicall_cromwell_gatk-4.6_wes_single_b37_*` and
contains generated `cbicall_cromwell.*.json` launch files.

For mtDNA, the run directory looks like:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_gatk-3.5_mit_single_rsrs_*/
  01_mtoolbox/
  02_browser/
```

Use [Run Comparison](run-comparison) when you want to compare two or more completed runs.
