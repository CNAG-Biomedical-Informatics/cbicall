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

### Method 2: Building from the source checkout (slow)

Clone CBIcall so Docker can build the exact checked-out revision:

```bash
git clone https://github.com/CNAG-Biomedical-Informatics/cbicall.git
cd cbicall
```

For a released version, check out its tag before building. For example:

```bash
git checkout 1.1.0
```

Then build the container from the repository root. The source revision and
CBIcall version are recorded in the image metadata:

The Dockerfile installs the CBIcall Python dependencies, including Snakemake,
plus pinned Nextflow, Cromwell, and WOMtool launchers. Bash workflows do not
need those engines, but the image can run the packaged Snakemake, Nextflow, and
Cromwell WES/WGS workflows without extra engine installation.

- **For Docker version 19.03 and above (supports buildx):**

  ```bash
  docker buildx build --load \
    -f docker/Dockerfile \
    --build-arg CBICALL_VERSION="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' src/cbicall/__about__.py)" \
    --build-arg VCS_REF="$(git rev-parse HEAD)" \
    -t cnag/cbicall:latest .
  ```

- **For Docker versions older than 19.03 (no buildx support):**

  ```bash
  docker build \
    -f docker/Dockerfile \
    --build-arg CBICALL_VERSION="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' src/cbicall/__about__.py)" \
    --build-arg VCS_REF="$(git rev-parse HEAD)" \
    -t cnag/cbicall:latest .
  ```

Inspect the revision packaged in an image with:

```bash
docker inspect cnag/cbicall:latest \
  --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

The Docker Hub multi-architecture image is built in the same way: GitHub
Actions checks out one repository commit and passes that checkout as the Docker
build context.

## Choose a Workflow Path

| Workflow path | Needs `/cbicall-data`? | Notes |
| --- | --- | --- |
| `workflow_provider: nf-core` | No CBIcall bundle required | Nextflow/nf-core manage their own test data, references, and containers. The selected runtime profile must still work in your environment. |
| Native CBIcall Bash/Snakemake/Nextflow/Cromwell WES/WGS/mtDNA | Yes | Mount the CBIcall resource bundle as `/cbicall-data` and set `CBICALL_DATA=/cbicall-data`. |

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

Run the installer from the CBIcall image while mounting a persistent host
directory:

```bash
mkdir -p /absolute/path/to/cbicall-data
docker run --rm \
  --volume /absolute/path/to/cbicall-data:/cbicall-data \
  cnag/cbicall:latest \
  cbicall install-resources --outdir /cbicall-data
```

To verify only the catalog-to-Google-Drive bundle identity before starting the large archive download:

```bash
docker run --rm --volume /absolute/path/to/cbicall-data:/cbicall-data \
  cnag/cbicall:latest cbicall install-resources --outdir /cbicall-data \
  --verify-resource-id-only
```

Google Drive can be restrictive with large files. If the Python download stalls or fails, print the manual download list:

```bash
docker run --rm --volume /absolute/path/to/cbicall-data:/cbicall-data \
  cnag/cbicall:latest cbicall install-resources --outdir /cbicall-data \
  --print-manual-download
```

Download every listed file into `/absolute/path/to/cbicall-data`, then let the script continue from those files. **This step can take time because it assembles, verifies, and extracts the full resource bundle.** On a typical VM or workstation disk, expect roughly 20-50 minutes after all parts are present; faster disks may be shorter.

```bash
docker run --rm --volume /absolute/path/to/cbicall-data:/cbicall-data \
  cnag/cbicall:latest cbicall install-resources --outdir /cbicall-data \
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

## Running and Interacting with the Container

For native CBIcall workflows:

```bash
# Please update '/absolute/path/to/cbicall-data' with your actual local data path
docker run -tid \
  --volume /absolute/path/to/cbicall-data:/cbicall-data \
  --env CBICALL_DATA=/cbicall-data \
  --env USERNAME=root \
  --name cbicall cnag/cbicall:latest
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

### Validate the mounted resources

The container was started with `CBICALL_DATA=/cbicall-data`, so all native
backends receive the mounted resource location without modifying packaged
workflow files.

Then confirm that CBIcall sees the mounted resources:

```bash
cbicall doctor
cbicall validate-parameters -p examples/input/param.yaml
```

## Performing integration tests

Inside the container, from the CBIcall repository root:

**WES**

```bash
cbicall test --wes-bash -t 1
```

**mtDNA**

```bash
cbicall test --mit-bash -t 1
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
