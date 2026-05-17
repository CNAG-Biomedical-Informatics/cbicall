# Integration Tests

CBIcall ships small integration tests for framework-level validation. They are intended to confirm that the CLI, YAML resolution, workflow registry, resource catalog, bundled workflows, and deterministic output comparison work together.

They do not replace biological or clinical validation of the underlying variant-calling methods.

## Minimal Check

From the repository root:

```bash
bin/cbicall validate-registry
bin/cbicall validate-resources
bin/cbicall doctor -p examples/input/param.yaml
bin/cbicall test --wes -t 1
```

| Check | What it confirms |
| --- | --- |
| `validate-registry` | The workflow registry conforms to its JSON Schema. |
| `validate-resources` | The resource catalog is well formed and compatible workflow keys exist in the registry. |
| `doctor` | The YAML resolves to a declared workflow, profile, pipeline implementation version, and selected resource. |
| `test --wes` | The bundled WES workflow runs and reproduces the shipped reference VCF under deterministic comparison rules. |

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
