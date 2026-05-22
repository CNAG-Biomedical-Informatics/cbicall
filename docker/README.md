# Containerized Installation

## Install CBIcall with Docker

CBIcall installation and CBIcall resource installation are separate steps. Pull
or build the CBIcall image first; download the native CBIcall resource bundle
only if the workflow path you choose needs it.

### Method 1: Installing from Docker Hub (fast)

Pull the latest Docker image from [Docker Hub](https://hub.docker.com/r/manuelrueda/cbicall):

```bash
docker pull manuelrueda/cbicall:latest
docker image tag manuelrueda/cbicall:latest cnag/cbicall:latest
```

### Method 2: Installing from Dockerfile (slow)

Download the `Dockerfile` from [GitHub](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/Dockerfile):

```bash
wget https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/main/docker/Dockerfile
```

Then build the container:

The Dockerfile installs the CBIcall Python dependencies, including
Snakemake, and a pinned Nextflow launcher. Bash workflows do not need either
engine, but the image can run the packaged Snakemake and Nextflow WES/WGS
workflows without extra engine installation.

- **For Docker version 19.03 and above (supports buildx):**

  ```bash
  docker buildx build --no-cache -t cnag/cbicall:latest .
  ```

- **For Docker versions older than 19.03 (no buildx support):**

  ```bash
  docker build --no-cache -t cnag/cbicall:latest .
  ```

## Choose a Workflow Path

| Workflow path | Needs `/cbicall-data`? | Notes |
| --- | --- | --- |
| `workflow_provider: nf-core` | No CBIcall bundle required | Nextflow/nf-core manage their own test data, references, and containers. The selected runtime profile must still work in your environment. |
| Native CBIcall Bash/Snakemake/Nextflow WES/WGS/mtDNA | Yes | Mount the CBIcall resource bundle as `/cbicall-data` and configure `DATADIR`. |

For Docker-based nf-core tests, make sure the selected nf-core profile can run
from your container setup. Many users run nf-core directly on the host with
Docker or on HPC with Singularity/Apptainer.

## Download the Resource Bundle for Native Workflows

> Note: this process can be lengthy.

Begin by downloading the required databases and software for native CBIcall
workflows. Save the data outside the container; this preserves it across
container restarts and lets you update the software without downloading the data
again. Skip this section when you only want to validate or run registered
nf-core provider workflows.

Install dependencies for Python 3:

```
pip3 install gdown
```

Finally, navigate to a directory where you want the databases stored and execute:

```bash
mkdir -p /absolute/path/to/cbicall-data
wget https://raw.githubusercontent.com/mrueda/cbicall/refs/heads/main/scripts/download_cbicall_bundle.py
python3 ./download_cbicall_bundle.py --outdir /absolute/path/to/cbicall-data
```

To verify only the catalog-to-Google-Drive bundle identity before starting the large archive download:

```bash
python3 ./download_cbicall_bundle.py \
  --outdir /absolute/path/to/cbicall-data \
  --verify-resource-id-only
```

Google Drive can be restrictive with large files. If the Python download stalls or fails, print the manual download list:

```bash
python3 ./download_cbicall_bundle.py \
  --outdir /absolute/path/to/cbicall-data \
  --print-manual-download
```

Download every listed file into `/absolute/path/to/cbicall-data`, then let the script continue from those files:

```bash
python3 ./download_cbicall_bundle.py \
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
Snakemake/Nextflow `config.yaml` files. In Docker, mount your host resource
directory as
`/cbicall-data` and point CBIcall to that container path:

```bash
sed -i 's|^DATADIR=.*|DATADIR=/cbicall-data|' workflows/bash/gatk-4.6/env.sh
sed -i 's|^DATADIR=.*|DATADIR=/cbicall-data|' workflows/bash/gatk-3.5/env.sh
sed -i 's|^datadir:.*|datadir: "/cbicall-data"|' workflows/snakemake/gatk-4.6/config.yaml
```

The native Nextflow config is a symlink to this shared GATK 4.6 backend config, so one edit updates both Snakemake and Nextflow native workflows.

## Running and Interacting with the Container

For native CBIcall workflows:

```bash
# Please update '/absolute/path/to/cbicall-data' with your actual local data path
docker run -tid --volume /absolute/path/to/cbicall-data:/cbicall-data -e USERNAME=root --name cbicall cnag/cbicall:latest
```

For nf-core-only validation, the CBIcall bundle mount is not required:

```bash
docker run -tid -e USERNAME=root --name cbicall cnag/cbicall:latest
```

The container/runtime setup must still support the nf-core profile you select.

To connect to the container:

```bash
docker exec -ti cbicall bash
```

Inside the container, confirm that CBIcall sees the mounted resources when
running native workflows:

```bash
bin/cbicall validate-resources
bin/cbicall validate-parameters -p examples/input/param.yaml
```

## Performing integration tests

Inside the container, from the CBIcall repository root:

**WES**

```bash
bin/cbicall test --wes-bash -t 1
```

**mtDNA**

```bash
bin/cbicall test --mit-bash -t 1
```

## System requirements

- OS/ARCH supported: **linux/amd64** and **linux/arm64**.
- Ideally a Debian-based distribution (Ubuntu or Mint), but any other (e.g., CentOS, OpenSUSE) should do as well.
- 16GB of RAM
- \>= 1 core (ideally i7 or Xeon).
- At least 100GB HDD.

## Common errors: Symptoms and treatment

  * Dockerfile:

          * DNS errors

            - Error: Temporary failure resolving 'foo'

              Solution: https://askubuntu.com/questions/91543/apt-get-update-fails-to-fetch-files-temporary-failure-resolving-error
