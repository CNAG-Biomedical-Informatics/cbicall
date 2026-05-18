# Integration Tests

CBIcall ships small integration tests that run the bundled example workflows and
compare their outputs against shipped references.

They do not replace biological or clinical validation of the underlying variant-calling methods.

<div className="cbicallNotePanel">
  <p><strong>Use this page when the question is:</strong> can this CBIcall installation run the shipped example workflows and reproduce the expected example outputs?</p>
</div>

This page is about execution tests. If the question is whether a parameters YAML
resolves before launch, use `validate-param` as described in
[Configuration Reference](../help/configuration-reference). If the question is
whether a custom resource catalog or installed resource directory is valid, use
[Resource Validation](resource-validation).

## Minimal Test

From the repository root:

```bash
bin/cbicall test --wes-bash -t 1
```

| Check | What it confirms |
| --- | --- |
| `test --wes-bash` | The required bundled Bash WES workflow runs and reproduces the shipped reference VCF under deterministic comparison rules. |

:::tip[Where to go next]
Use this page to run the shipped examples. For parameter-file checks, see
[Configuration Reference](../help/configuration-reference). For comparing
repeated runs, see [Run Comparison](run-comparison).
:::

## Test Commands

```bash
bin/cbicall test --wes-bash -t 1
bin/cbicall test --wes-snakemake -t 1
bin/cbicall test --wes-nextflow -t 1
bin/cbicall test --mit-bash -t 1
bin/cbicall test --all -t 1
```

| Command | Use |
| --- | --- |
| `bin/cbicall test --wes-bash -t 1` | Required Bash WES integration test. Run this first. |
| `bin/cbicall test --wes-snakemake -t 1` | Optional Snakemake WES test. Requires `snakemake` on `PATH` and compares the resulting VCF to the same Bash reference VCF. |
| `bin/cbicall test --wes-nextflow -t 1` | Optional Nextflow WES test. Requires `nextflow` on `PATH` and compares the resulting VCF to the same Bash reference VCF. |
| `bin/cbicall test --mit-bash -t 1` | Optional mtDNA Bash integration test after the WES path is working. |
| `bin/cbicall test --all -t 1` | Run all bundled integration examples. Optional engine tests are skipped when their engine is not installed. |

:::note[Workflow engine dependencies]
Snakemake and Nextflow are not part of the CBIcall resource bundle. Install them
in the runtime environment before running their backend-specific tests.
:::

## Outputs

The test command prints the run directory, workflow log, `run-report.json`, launcher log, and output file used for comparison.

:::tip[VCF hashes during tests]
WES tests compute normalized [SHA-256](https://en.wikipedia.org/wiki/SHA-2)
values directly from the shipped reference
VCF and the newly generated VCF after removing headers and sorting variant
records. The test does not rely on stored `03_stats/*.vcf.sha256.txt` files,
because those sidecar files can be absent or stale; they are kept for audit, run
reports, and `compare-runs`.
:::

For WES, the run directory looks like:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_b37_gatk-4.6_*/
  01_bam/
  02_varcall/
  03_stats/
  logs/
```

The optional Snakemake WES test uses the same input and expected VCF records, but
the run directory starts with `cbicall_snakemake_wes_single_b37_gatk-4.6_*`.

The optional Nextflow WES test uses the same input and expected VCF records, but
the run directory starts with `cbicall_nextflow_wes_single_b37_gatk-4.6_*`.

For mtDNA, the run directory looks like:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_mit_single_rsrs_gatk-3.5_*/
  01_mtoolbox/
  02_browser/
```

Use [Run Comparison](run-comparison) when you want to compare two or more completed runs.
