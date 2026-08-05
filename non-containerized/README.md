# Non-containerized installation

Install the CBIcall Python package first. The large reference and software
bundle is separate and is required only for the bundled `cbicall-core`
workflows.

## 1. Install CBIcall

PyPI is recommended for normal use:

```bash
python3 -m pip install --upgrade cbicall
cbicall --version
```

Install optional Snakemake and MultiQC dependencies when needed:

```bash
python3 -m pip install --upgrade "cbicall[all]"
```

Nextflow and Cromwell are external executables and are not installed as Python
dependencies. Bash workflows do not require either engine.

## 2. Choose the workflow source

| Workflow source | CBIcall bundle required? | Additional runtime |
| --- | --- | --- |
| `workflow_provider: cbicall-core` (default) | Yes | The selected backend: Bash, Snakemake, Nextflow, or Cromwell. mtDNA also requires x86_64 Linux. |
| `workflow_provider: nf-core` | No | Nextflow and a working nf-core runtime profile, such as Docker or Apptainer. |

To test external nf-core support without installing the CBIcall bundle:

```bash
cbicall test --nf-core-demo -t 4
```

## 3. Install resources for `cbicall-core`

Choose a persistent directory and install the registered bundle:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
mkdir -p "$CBICALL_DATA"
cbicall install-resources --outdir "$CBICALL_DATA"
```

The installer downloads and assembles the bundle, verifies its catalog-declared
checksums, extracts `Databases/` and `NGSutils/`, and writes
`cbicall-resource-installation.json` with installation provenance.

Set `CBICALL_DATA` in every shell or scheduler environment used to run bundled
workflows:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
cbicall doctor
```

The same environment variable is propagated to all four bundled WES/WGS
backends and recorded in the execution contract.

## 4. Run an integration test

```bash
cbicall test --wes-bash -t 1
```

On x86_64, the mtDNA workflow can also be tested with:

```bash
cbicall test --mit-bash -t 1
```

Installed-package tests use packaged fixtures and create a temporary workspace.
Use `--workspace /path/to/new-or-empty-directory` to retain the test outputs in
a chosen location.

## Source checkout for development

```bash
git clone https://github.com/CNAG-Biomedical-Informatics/cbicall.git
cd cbicall
python3 -m pip install -e ".[all,test]"
pytest
```

Use `cbicall` after installation. The repository launcher `bin/cbicall` remains
available inside a source checkout.

<details>
<summary>Manual resource-download recovery</summary>

Google Drive can restrict large automated downloads. Print the registered file
list when automatic download fails:

```bash
cbicall install-resources \
  --outdir "$CBICALL_DATA" \
  --print-manual-download
```

Place every listed file in `$CBICALL_DATA`, then resume verification, assembly,
and extraction:

```bash
cbicall install-resources \
  --outdir "$CBICALL_DATA" \
  --skip-download
```

Use `--verify-resource-id-only` to verify the small catalog-pinned resource
identifier before downloading the archive. Add `--remove-parts` after successful
verification when disk space is limited.

</details>

## System requirements

- Linux on amd64 or arm64; macOS can be used through a Linux VM.
- Python 3.8 or newer.
- Java 17 for GATK 4.6 workflows.
- `libncurses.so.5` and `libtinfo.so.5` for bundled legacy tools.
- Sufficient storage for the resource bundle and analysis inputs; allow at least
  100 GB for installation and the shipped WES test.

The bundled MToolBox mtDNA workflow supports x86_64 only.
