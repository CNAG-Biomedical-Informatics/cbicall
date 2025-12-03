# Quickstart

This page provides a minimal introduction to using CBIcall.  
It focuses on **basic commands**, **help options**, and **a simple built-in test run**.

---

## 1. Basic usage

You can check the command-line interface using:

```bash
bin/cbicall --help
```

To display the version:

```bash
bin/cbicall --version
```

---

## 2. Minimal command example

Once you have a parameters file (YAML), you can run:

```bash
bin/cbicall -p param.yaml -t 8
```

- `-p` selects the param file  
- `-t` sets the number of threads

This is the minimal invocation required for any workflow.

---

## 3. Run a small WES job (workstation)

CBIcall includes small test configurations to validate that your installation works.

Example:

```bash
cd examples/input
../../bin/cbicall -p wes_single.yaml -t 4
```

When the job completes successfully (<2 min), you should see:

```
CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_gatk-4.6_*/
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

## 4. Run a small mtDNA job (workstation)

Example:

```bash
cd examples/input # ommit if you are already there
../../bin/cbicall -p mit_single.yaml -t 4
```

When the job completes successfully (~ 5 min), you should see:

```
CNAG999_exome/CNAG99901P_ex/cbicall_bash_mit_single_gatk-3.5_*/
  01_mtoolbox/
  02_browser/
```

We recommend you to check:

```
02_browser/README.txt
```

---


## 5. Next steps

For a complete, real-case workflow demonstration, see:

[➡️ End-to-end example WES](end-to-end-example-wes-single.md){ .md-button .md-button--primary }
