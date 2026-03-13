# Quickstart

This page is the shortest path to a validated local test run using the example data shipped in the repository.

If you first need to decide which installation method or workflow applies to you, start with [Choose your path](choose-your-path.md).

---

## 1. Confirm the CLI is available

You can check the command-line interface using:

```bash
bin/cbicall --help
```

To display the version:

```bash
bin/cbicall --version
```

---

## 2. Run the shipped WES smoke test

This is the most reliable first run because it uses the same example layout used by the integration test script.

```bash
cd examples/input
./run_tests.sh --wes
```

What this does:

- runs `cbicall` on the bundled WES example
- compares the resulting VCF against the shipped reference output
- exits non-zero if the result differs

When the test completes successfully, you should have a run directory like:

```
CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_b37_gatk-4.6_*/
  01_bam/
  02_varcall/
  03_stats/
  logs/
```
The final `VCF` will be located at:

```
02_varcall/CNAG99901P.hc.QC.vcf.gz
```

---

## 3. Run the shipped mtDNA smoke test

Run this after the WES example is working. The mtDNA example depends on the WES/WGS-style project structure and uses the bundled mtDNA example configuration.

```bash
cd examples/input
./run_tests.sh --mit
```

When the test completes successfully, you should have a run directory like:

```text
CNAG999_exome/CNAG99901P_ex/cbicall_bash_mit_single_rsrs_gatk-3.5_*/
  01_mtoolbox/
  02_browser/
```

Start by checking:

```
02_browser/README.txt
```

---

## 4. Run CBIcall directly with a YAML file

Once the smoke tests work, the minimal direct invocation is:

```bash
bin/cbicall -p param.yaml -t 8
```

- `-p` selects the parameter file
- `-t` sets the number of threads

---

## 5. Next steps

- Use [Choose your path](choose-your-path.md) if you are still deciding between installation methods or workflows.
- For a realistic WES walkthrough, see [End-to-end example WES](end-to-end-example-wes.md).
- For mtDNA, see [End-to-end example mtDNA](end-to-end-example-mit.md).

[➡️ End-to-end example WES](end-to-end-example-wes.md){ .md-button .md-button--primary }
