# GIAB Benchmarking

Genome in a Bottle (GIAB) benchmarking is analytical validation of the selected
variant-calling workflow, not validation of CBIcall itself. Use this page to
record truth-set comparisons for native CBIcall WES/WGS workflows.

## Recommended Reporting

Document the benchmark as a reproducible contract:

| Field | Record |
| --- | --- |
| GIAB sample | e.g. HG001, HG002, HG003, or HG004 |
| Sequencing data | FASTQ accessions or local source |
| Genome build | b37, GRCh38, or other |
| Truth VCF and BED | exact release and checksum when available |
| CBIcall workflow | pipeline, mode, workflow backend, software stack, registry version |
| Resource bundle | resource key, version, fingerprint |
| Comparison tool | hap.py, rtg vcfeval, bcftools, or another tool |
| Metrics | precision, recall, F1, Ti/Tv, PASS counts, or other reported values |

## Minimal Evidence

Keep the CBIcall audit files together with the benchmarking output:

```text
run-report.json
run-report.html
log.json
cbicall-execution-contract.json
benchmark command log
benchmark metrics table
```

:::important[Scope]
GIAB benchmarking evaluates the analytical behavior of the selected workflow and
resource bundle against a truth set. It complements, but does not replace,
CBIcall run comparison, which audits whether executions used the same workflow,
resources, runtime context, and final-output fingerprints.
:::
