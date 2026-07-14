# Non-containerized installation

CBIcall and its external bioinformatics resource bundle are installed
separately. The Python package contains the driver, workflow definitions,
schemas, catalogs, report assets, and integration-test fixtures. Reference
genomes and third-party bioinformatics tools remain in the external bundle.

## Install CBIcall

### Method 1: PyPI

Install the core command:

```bash
python3 -m pip install cbicall
cbicall --version
```

Install the optional Snakemake and MultiQC Python dependencies when needed:

```bash
python3 -m pip install "cbicall[all]"
```

Nextflow and Cromwell are external executables and are not installed as Python
dependencies. Bash workflows need neither engine.

### Method 2: Source checkout

Use an editable install for development:

```bash
git clone https://github.com/CNAG-Biomedical-Informatics/cbicall.git
cd cbicall
python3 -m pip install -e ".[all,test]"
pytest
```

After either installation method, use `cbicall` directly. The repository
launcher `bin/cbicall` remains available for source-checkout compatibility.

## Choose a workflow path

| Workflow path | CBIcall bundle required? | Additional runtime |
| --- | --- | --- |
| External nf-core provider | No | Nextflow and the selected nf-core profile, such as Docker or Apptainer. |
| Native Bash, Snakemake, Nextflow, or Cromwell WES/WGS | Yes | The selected workflow backend and the CBIcall resource bundle. |
| Native Bash mtDNA | Yes | x86_64 host and the CBIcall resource bundle. |

The nf-core demo can be validated and launched without the native bundle:

```bash
cbicall validate-parameters -p examples/input/nf-core-demo.yaml --no-color
cbicall run -p examples/input/nf-core-demo.yaml -t 4 --no-color
```

## Install resources for native workflows

Choose a persistent directory for the external bundle:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
mkdir -p "$CBICALL_DATA"
cbicall install-resources --outdir "$CBICALL_DATA"
```

To verify only the small catalog-pinned resource identifier before downloading
the large archive:

```bash
cbicall install-resources \
  --outdir "$CBICALL_DATA" \
  --verify-resource-id-only
```

Google Drive may restrict large automated downloads. Print the manual download
list when necessary:

```bash
cbicall install-resources \
  --outdir "$CBICALL_DATA" \
  --print-manual-download
```

Place every listed file in `$CBICALL_DATA`, then resume assembly, checksum
verification, extraction, and manifest creation:

```bash
cbicall install-resources \
  --outdir "$CBICALL_DATA" \
  --skip-download
```

The installer:

- downloads missing split files when possible;
- reassembles the archive;
- verifies the catalog-declared checksums;
- validates the optional SHA-256-pinned resource identifier;
- extracts `Databases/` and `NGSutils/`;
- writes `cbicall-resource-installation.json` with installation provenance.

Add `--remove-parts` after successful verification when disk space is limited.

## Configure the resource location

Set `CBICALL_DATA` in every shell or scheduler environment used to launch
native workflows:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
cbicall validate-resources
cbicall validate-parameters -p examples/input/param.yaml
```

The same variable is passed consistently to native Bash, Snakemake, Nextflow,
and Cromwell workflows and is recorded in the execution contract. Checked-in
backend paths remain fallbacks for existing source and institutional profiles.

## Run integration tests

```bash
cbicall test --wes-bash -t 1
cbicall test --mit-bash -t 1
```

Source-checkout tests write under `examples/input`. Installed-package tests
stage the packaged fixtures in a temporary directory and print that directory
so generated reports can be inspected.

## System requirements

- Linux on amd64 or arm64; macOS can be used with a Linux VM.
- Python 3.8 or newer.
- Java 17 for current GATK 4.6 workflows.
- `libncurses.so.5` and `libtinfo.so.5` compatibility libraries for bundled
  legacy tools.
- At least 16 GB RAM and 100 GB disk for the native resource bundle and test
  workflows.

The mtDNA workflow uses MToolBox and is supported on x86_64 only.
