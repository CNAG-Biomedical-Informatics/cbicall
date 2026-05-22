# Non-containerized installation

Feel free to work with your preferred virtual environment. For this document, we'll move directly to the setup steps.

### Method 1: Download from GitHub

Use `git clone` to get the latest (stable) version:

```bash
git clone https://github.com/CNAG-Biomedical-Informatics/cbicall.git
cd cbicall
```

If you only need to update to the latest version do:

```bash
git pull
```

Install dependencies for Python 3:

```
pip3 install -r requirements.txt
```

> **Note:** If you are installing `cbicall` in an HPC environment for shared use, we recommend installing the required Python 3 modules in a central location. This allows users to simply do:

```bash
# Load Python + modules
module load Python/3.10.8-GCCcore-12.2.0
export PYTHONPATH="/software/biomed/cbi_py3/lib/python3.10/site-packages:${PYTHONPATH}"
``` 

Testing the deployment:

```bash
pytest
```

### Choose a workflow path

CBIcall can now be installed and used before downloading the large CBIcall
germline resource bundle.

| Workflow path | Resource bundle required? | Extra runtime |
| --- | --- | --- |
| `workflow_provider: nf-core` | No | Nextflow plus the selected nf-core runtime profile, such as Docker or Singularity/Apptainer. |
| Native CBIcall Bash/Snakemake/Nextflow WES/WGS/mtDNA | Yes | The CBIcall resource bundle installed as `DATADIR`. |

For a quick nf-core test from a source checkout:

```bash
cd examples/input
../../bin/cbicall validate-parameters -p nf-core-demo.yaml --no-color
../../bin/cbicall run -p nf-core-demo.yaml -t 4 --no-color
```

This uses nf-core's own test data and does not require the CBIcall-provided
bundle. Install or load Nextflow and the selected container/runtime profile
before running it.

### Download the Resource Bundle for Native Workflows

> Note: this process can be lengthy.

This section is required for the native CBIcall WES/WGS/mtDNA workflows. It is
not required for nf-core provider workflows where resources are managed by
Nextflow/nf-core.

Choose a directory where the databases and bundled external tools should be installed. This directory will become your `DATADIR`.

```bash
mkdir -p /absolute/path/to/cbicall-data
python3 $path_to_cbicall/scripts/download_cbicall_bundle.py --outdir /absolute/path/to/cbicall-data
```

Replace `$path_to_cbicall` with your CBIcall installation path.

To verify only the catalog-to-Google-Drive bundle identity before starting the large archive download:

```bash
python3 $path_to_cbicall/scripts/download_cbicall_bundle.py \
  --outdir /absolute/path/to/cbicall-data \
  --verify-resource-id-only
```

Google Drive can be restrictive with large files. If the Python download stalls or fails, print the manual download list:

```bash
python3 $path_to_cbicall/scripts/download_cbicall_bundle.py \
  --outdir /absolute/path/to/cbicall-data \
  --print-manual-download
```

Download every listed file into `/absolute/path/to/cbicall-data`, then let the script continue from those files:

```bash
python3 $path_to_cbicall/scripts/download_cbicall_bundle.py \
  --outdir /absolute/path/to/cbicall-data \
  --skip-download
```

The script will:

- download missing split files when possible
- reassemble `data.tar.gz`
- verify the split parts or assembled archive with `data.tar.gz.md5`
- load the CBIcall resource catalog, locally or from the catalog URL
- optionally verify a small GDrive resource identifier file such as `cbicall-resource-id.json`
- rename the verified archive using the bundle identity, for example `cbicall-germline-resources-v1.tar.gz`
- extract the archive into `DATADIR`
- write `cbicall-resource-installation.json` with the installed bundle provenance

If disk space is tight and the checksum has passed, add `--remove-parts` to remove `data.tar.gz.part-*` after assembly.

CBIcall keeps the rich resource registry in `resources/cbicall-resource-catalog.json`. The GDrive bundle only needs a small identifier file, for example `cbicall-resource-id.json` containing `{"resource_key": "cbicall-germline-resources-v1"}`. When that identifier file is available, the registry can store its Google Drive file ID and SHA-256 so the downloader can verify that the remote bundle matches the local CBIcall catalog entry.

### Point Native Workflows to your resource directory

Native CBIcall workflows read resource paths from Bash `env.sh` files and from
Snakemake/Nextflow `config.yaml` files. In a non-containerized installation,
point those files
to the host directory where you installed the CBIcall-provided resource bundle:

```bash
export CBICALL_DATA="/absolute/path/to/cbicall-data"

sed -i "s|^DATADIR=.*|DATADIR=${CBICALL_DATA}|" workflows/bash/gatk-4.6/env.sh
sed -i "s|^DATADIR=.*|DATADIR=${CBICALL_DATA}|" workflows/bash/gatk-3.5/env.sh
sed -i "s|^datadir:.*|datadir: \"${CBICALL_DATA}\"|" workflows/snakemake/gatk-4.6/config.yaml
```

The native Nextflow config is a symlink to this shared GATK 4.6 backend config, so one edit updates both Snakemake and Nextflow native workflows.

Confirm that CBIcall sees the configured resources:

```bash
bin/cbicall validate-resources
bin/cbicall validate-parameters -p examples/input/param.yaml
```

Ok, finally we are going to install `Java 8` in case you don't have it already:

```
sudo apt install openjdk-8-jdk # In some systems you might need Java 17 -> openjdk-17-jre
```

## Performing integration tests

Once you are in the root directory of the repo:

**WES**:

```bash
bin/cbicall test --wes-bash -t 1
```

**mtDNA**:

```bash
bin/cbicall test --mit-bash -t 1
```

## System requirements

- OS/ARCH supported: **linux/amd64** and **linux/arm64**.
- Ideally a Debian-based distribution (Ubuntu or Mint), but any other (e.g., CentOS, OpenSUSE) should do as well (untested).
- Python >= 3.8
- Java 8
- 16GB of RAM
- \>= 1 core (ideally i7 or Xeon).
- At least 100GB HDD.

## Platform Compatibility
This distribution is written in Python 3 and is intended to run on any platform supported by Python 3. It has been tested on Debian Linux and macOS. Please report any issues.
