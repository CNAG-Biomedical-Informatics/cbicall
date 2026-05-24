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
[Resource Validation](resource-validation).

## Minimal Test

From the repository root:

```bash
bin/cbicall test --wes-bash -t 1
```

On HPC, pass the same runtime profile you use for normal runs:

```bash
bin/cbicall test --wes-bash -t 1 --runtime-profile cnag-hpc
```

| Check | What it confirms |
| --- | --- |
| `test --wes-bash` | The required bundled Bash WES workflow runs and reproduces the expected normalized VCF hash declared by the contract fixture. |
| `test --release` | Runs Bash plus available native WES backends and checks that their normalized final VCF content matches Bash. |

## Release Check

Use release mode before tagging or publishing an image:

```bash
bin/cbicall test --release -t 1 --runtime-profile local
```

On HPC, use the same runtime profile as normal runs:

```bash
bin/cbicall test --release -t 1 --runtime-profile cnag-hpc
```

`--release` uses Bash as the required native baseline and compares available
Snakemake, Nextflow, and Cromwell WES runs against it. It compares the
**normalized final VCF hash**, not logs, workflow files, execution contracts, or
full output inventories. Missing optional backends are skipped, but at least one
non-Bash backend must be available for the release check to pass.

The final summary is the key audit line:

```text
Release equivalence summary
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

| Command | Workflow path | Needs CBIcall bundle? | Extra requirement | Contract check |
| --- | --- | --- | --- | --- |
| `bin/cbicall test --wes-bash -t 1` | **Native WES**, Bash | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | <span className="cbicallTestBadge cbicallTestBadgeNeutral">none</span> | Run report fields, expected files, **normalized VCF hash** |
| `bin/cbicall test --wes-snakemake -t 1` | **Native WES**, Snakemake | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | `snakemake` on `PATH` | Run report fields, expected files, **normalized VCF hash** |
| `bin/cbicall test --wes-nextflow -t 1` | **Native WES**, Nextflow | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | `nextflow` on `PATH` | Run report fields, expected files, **normalized VCF hash** |
| `bin/cbicall test --wes-cromwell -t 1` | **Native WES**, Cromwell | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | `CROMWELL_JAR` or `cromwell` on `PATH` | Generated inputs/options/metadata, run report fields, expected files, **normalized VCF hash** |
| `bin/cbicall test --mit-bash -t 1` | **Native mtDNA**, Bash | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | x86_64 host | Run report fields, expected files, **prioritized variants hash**, **raw JSON hash** |
| `bin/cbicall test --nf-core-demo -t 4` | **nf-core/demo** | <span className="cbicallTestBadge cbicallTestBadgeNo">X bundle</span> | Nextflow plus selected nf-core runtime profile | Generated params/config, run reports, **pipeline info**, **MultiQC anchors** |
| `bin/cbicall test --nf-core-sarek -t 4` | **nf-core/Sarek** | <span className="cbicallTestBadge cbicallTestBadgeNo">X bundle</span> | Nextflow plus selected nf-core runtime profile and Sarek inputs/resources | Generated params/config, run reports, **pipeline info**, **MultiQC anchors**, declared canonical outputs when produced |
| `bin/cbicall test --release -t 1` | **Native WES backend equivalence** | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | At least one non-Bash native backend available | Bash baseline plus available native WES backends; **same normalized final VCF** required |
| `bin/cbicall test --all -t 1` | **Native tests only** | <span className="cbicallTestBadge cbicallTestBadgeYes">V bundle</span> | Optional backends skipped if missing | Runs WES Bash, WES Snakemake, WES Nextflow, WES Cromwell, and mtDNA contracts |

Useful variants:

```bash
bin/cbicall test --release -t 1 --runtime-profile cnag-hpc
bin/cbicall test --wes-bash -t 1 --runtime-profile cnag-hpc
bin/cbicall test --nf-core-demo -t 4 --keep-external-work
```

**`--runtime-profile`** is for native CBIcall workflow profiles.
**`--keep-external-work`** keeps heavy Nextflow state such as `work/` and
`.nextflow/` after an nf-core contract test so you can inspect task-level
details.

:::note[Workflow backend dependencies]
Snakemake, Nextflow, and Cromwell are not part of the CBIcall resource bundle.
Install them in the runtime environment before running their backend-specific
tests. For Cromwell, set `CROMWELL_JAR=/path/to/cromwell.jar` or put a
`cromwell` launcher on `PATH`.
For Nextflow or nf-core tests on CNAG HPC, load Nextflow first:

```bash
module load Nextflow/25.10.2
```
:::

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
