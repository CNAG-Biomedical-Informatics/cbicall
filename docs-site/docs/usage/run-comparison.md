# Run Comparison

`cbicall compare-runs` compares completed CBIcall run directories or
`run-report.json` files. It is intended for reproducibility checks across
repeated local runs, HPC runs, container runs, or cloud runs.

The command **does not try to prove that two biological analyses are equivalent**.
It gives an audit trail for the framework layer: which CBIcall version ran,
which registered pipeline was resolved, which workflow files were executed,
which resource identity was selected, and whether normalized VCF output
fingerprints and output file inventories match when available.

## Minimal Audit

```bash
bin/cbicall compare-runs run_a/ run_b/ \
  --output compare-report.txt \
  --html compare-report.html
```

Keep these files from each completed run:

| File | Why it matters |
| --- | --- |
| `run-report.json` | Compact provenance report used by `compare-runs`. |
| `log.json` | Full resolved configuration, runtime parameters, and resource details. |
| Workflow log | Execution log for the selected Bash, Snakemake, or Nextflow backend. |
| `03_stats/*.vcf.sha256.txt` | Normalized VCF fingerprint report when produced by the workflow. |

For a concise methods audit, archive `compare-report.txt` together with the two
`run-report.json` files. The optional HTML file is useful for manual browsing
but contains the same comparison content as the text report.

With two runs, CBIcall prints a direct pairwise comparison. With three or more
runs, the first run is used as the baseline and the remaining runs are compared
against it:

```bash
bin/cbicall compare-runs baseline_run/ repeat_1/ repeat_2/ repeat_3/ \
  --output compare-report.txt \
  --html compare-report.html
```

## What Is Compared

| Layer | Fields |
| --- | --- |
| Framework | CBIcall version, Python version, and workflow engine version recorded in `run-report.json`. |
| Pipeline | Workflow key, pipeline implementation version, entrypoint, and workflow fingerprint. |
| Workflow files | Entrypoint and helper/config file paths plus their SHA-256 values. |
| Resources | Resource key and resource fingerprint from the selected resource catalog entry. |
| Outputs | Run-directory file inventory fingerprint and normalized VCF fingerprints when `03_stats/*.vcf.sha256.txt` is present. |

:::note[Runtime fingerprints]
Workflow and resource fingerprints are **computed at runtime** from the files and
catalog entries actually resolved for the run. CBIcall deliberately does not
store expected workflow hashes in the workflow registry or resource catalog:
otherwise every harmless comment or formatting edit in a workflow script would
force metadata churn before the next run.
:::

The workflow fingerprint is computed from the resolved workflow files. Any byte
change in the entrypoint, helpers, Snakefile, or config files changes this
fingerprint, including comment-only edits. This is deliberate: it tells the
auditor that the implementation used for the second run was not exactly the same
implementation used for the first run.

The output fingerprint is different. It is computed from **normalized VCF records**,
not from the raw VCF file bytes. This avoids reporting false differences caused
only by VCF header timestamps, command lines, or compression metadata.

The file inventory fingerprint is also path-based, not content-based. It hashes
the sorted list of relative file paths in the run directory, excluding
`run-report.json` itself. It is useful for spotting layout differences, including
expected differences caused by settings such as `cleanup_bam`.

The status vocabulary is intentionally small:

| Status | Meaning |
| --- | --- |
| `same` | Values or fingerprints match. |
| `different` | Values or fingerprints exist in all compared runs but differ. |
| `missing` | Evidence is present in only some runs. |
| `not available` | Evidence is not recorded in any compared run. |

## How To Read The Result

Use this order when auditing two runs:

1. Check **Framework**. A different CBIcall version means the execution driver
   changed between runs. A different Python or workflow engine version means
   the runtime stack changed.
2. Check **Pipeline** and **Workflow files**. A different workflow fingerprint
   means the resolved workflow implementation changed. Inspect the listed file
   fingerprints to locate the changed file.
3. Check **Resources**. A different resource key or hash means the selected
   external dependency set was not the same.
4. Check **Outputs**. A different file inventory means the run directories do
   not contain the same relative file layout. Matching normalized VCF
   fingerprints indicate that the
   compared variant records match under CBIcall's deterministic VCF comparison
   rules.

If the workflow fingerprint changed but the normalized VCF fingerprint is the
same, the two runs produced the same compared VCF records despite an
implementation text change. That is useful audit evidence, but the change should
still be inspected before claiming full execution identity.

## Reports

The text report is the canonical audit artifact because it is easy to diff,
archive, and attach to review material. The HTML report is a static rendering of
the same information for browsing.

![Screenshot of the CBIcall HTML run comparison report showing summary counters, baseline run metadata, and status pills for framework, pipeline, and workflow file checks.](/img/compare-runs-html-report.png)

## Interpretation

A changed workflow fingerprint means the resolved workflow files are not
byte-identical. That is audit evidence, not automatically a failed analysis. For
example, editing a comment in a Bash workflow changes the workflow hash but may
leave the normalized VCF fingerprint unchanged.

For output reproducibility, prioritize the **normalized VCF fingerprint**. Use
the **file inventory fingerprint** to audit whether the run directory layout
matches. For
implementation provenance, inspect the **workflow and helper file fingerprints**.
