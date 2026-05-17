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
bin/cbicall test --wes -t 1
```

| Check | What it confirms |
| --- | --- |
| `test --wes` | The bundled WES workflow runs and reproduces the shipped reference VCF under deterministic comparison rules. |

:::tip[Where to go next]
Use this page to run the shipped examples. For parameter-file checks, see
[Configuration Reference](../help/configuration-reference). For comparing
repeated runs, see [Run Comparison](run-comparison).
:::

## Test Commands

```bash
bin/cbicall test --wes -t 1
bin/cbicall test --mit -t 1
bin/cbicall test --all -t 1
```

| Command | Use |
| --- | --- |
| `bin/cbicall test --wes -t 1` | Fast WES integration test. Run this first. |
| `bin/cbicall test --mit -t 1` | mtDNA integration test after the WES path is working. |
| `bin/cbicall test --all -t 1` | Run all bundled integration examples. |

## Outputs

The test command prints the run directory, workflow log, `run-report.json`, launcher log, and output file used for comparison.

For WES, the run directory looks like:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_b37_gatk-4.6_*/
  01_bam/
  02_varcall/
  03_stats/
  logs/
```

For mtDNA, the run directory looks like:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_mit_single_rsrs_gatk-3.5_*/
  01_mtoolbox/
  02_browser/
```

Use [Run Comparison](run-comparison) when you want to compare two or more completed runs.
