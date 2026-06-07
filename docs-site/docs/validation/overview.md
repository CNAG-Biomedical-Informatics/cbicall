# Validation & Reproducibility

CBIcall separates framework validation from biological validation of the
underlying variant-calling methods. This section collects the evidence and tools
used to audit CBIcall executions.

| Question | Page |
| --- | --- |
| Can this installation run the shipped example workflows? | [Integration Tests](integration-tests) |
| Do repeated runs match across machines or environments? | [Cross-Environment Reproducibility](cross-environment) |
| How do I compare completed run reports? | [Run Comparison](run-comparison) |
| Are installed resource bundles compatible with selected workflows? | [Resource Validation](../usage/resource-validation) |
| How do native WES calls compare with a truth set? | [GIAB Benchmarking](giab) |

:::important[Scope]
Integration tests and run comparison validate the CBIcall execution contract:
parameters, workflow resolution, resource identity, runtime evidence, and output
fingerprints. They do not replace biological or clinical validation of GATK,
MToolBox, nf-core/Sarek, or other upstream methods.
:::

## Evidence Layers

CBIcall records reproducibility evidence at several layers:

- `log.json` records resolved parameters, runtime profile, workflow selection,
  and resource identity.
- `cbicall-execution-contract.json` records the backend-ready launch plan.
- `run-report.json` records workflow fingerprints, resource fingerprints,
  output inventories, final-output hashes, software/runtime evidence, and status.
- `cbicall compare-runs` compares completed reports across independent runs.

For output reproducibility, prioritize normalized final-output fingerprints. File
inventories and logs are useful audit evidence, but they can differ between
environments even when final VCF records match.
