# Non-containerized installation

Feel free to work with your preferred virtual environment. For this document, we'll move directly to the setup steps.

### Method 1: Download from GitHub

Use `git clone` to get the latest (stable) version:

```bash
git clone https://github.com/CNAG-Biomedical-Informatics/cbicall.git
cd cbicall
```

If you only new to update to the lastest version do:

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

### Downloading Required Databases and Software

> Note: this process can be lenghty.

Choose a directory where the databases and bundled external tools should be installed. This directory will become your `DATADIR`.

```bash
mkdir -p /absolute/path/to/cbicall-data
python3 $path_to_cbicall/scripts/01_download_external_data.py --outdir /absolute/path/to/cbicall-data
```

Replace `$path_to_cbicall` with your CBIcall installation path.

Google Drive can be restrictive with large files. If the Python download stalls or fails, print the manual download list:

```bash
python3 $path_to_cbicall/scripts/01_download_external_data.py \
  --outdir /absolute/path/to/cbicall-data \
  --print-manual-download
```

Download every listed file into `/absolute/path/to/cbicall-data`, then let the script continue from those files:

```bash
python3 $path_to_cbicall/scripts/01_download_external_data.py \
  --outdir /absolute/path/to/cbicall-data \
  --skip-download
```

The script will:

- download missing split files when possible
- reassemble `data.tar.gz`
- verify `data.tar.gz` with `data.tar.gz.md5`
- read the local CBIcall resource registry
- optionally verify a small GDrive bundle identifier file such as `cbicall-bundle-id.json`
- rename the verified archive using the bundle identity, for example `cbicall-germline-resources-v1.tar.gz`
- extract the archive into `DATADIR`
- write `cbicall-resource-installation.json` with the installed bundle provenance

If disk space is tight and the checksum has passed, add `--remove-parts` to remove `data.tar.gz.part-*` after assembly.

CBIcall keeps the rich resource registry in `resources/cbicall-resource-catalog.json`. The GDrive bundle only needs a small identifier file, for example `cbicall-bundle-id.json` containing `{"bundle": "cbicall-germline-resources-v1"}`. When that identifier file is available, the registry can store its Google Drive file ID and SHA-256 so the downloader can verify that the remote bundle matches the local CBIcall catalog entry.

Finally, in the `cbicall` repo:

Change the `DATADIR` variable in `workflows/bash/*/env.sh` and `workflows/snakemake/*/config.yaml` so that it matches the location of your downloaded data.

Ok, finally we are going to install `Java 8` in case you don't have it already:

```
sudo apt install openjdk-8-jdk # In some systems you might need Java 17 -> openjdk-17-jre
```

## Performing integration tests

Once you are in the root directory of the repo:

**WES**:

```bash
cd examples/input
./run_tests.sh --wes
```

**mtDNA**:

```bash
cd examples/input
./run_tests.sh --mit
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
