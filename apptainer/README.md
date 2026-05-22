# Containerized Installation (HPC with Apptainer / Singularity)

This section describes how to install and run **CBIcall** on High-Performance Computing (HPC) systems where Docker is not available, using **Apptainer (formerly Singularity)**.

Apptainer can execute Docker images directly, but container images are **read-only by design**. As a result, databases, configuration files, and workflows must be stored outside the container and bind-mounted at runtime.

## Storage and Environment Model

An Apptainer setup has three separate paths. They can live in `$HOME`, a project
directory, or a shared software area depending on your site policy and quota.

| Variable | Host path | Bound inside container | Purpose |
| --- | --- | --- | --- |
| `SIF_IMAGE` | Path to the CBIcall `.sif` file you pulled | image itself | Immutable CBIcall application image |
| `CBICALL_WRITABLE` | Writable copy of the CBIcall installation | `/usr/share/cbicall` | Place where profiles, examples, `env.sh`, configs, and run files can be edited |
| `CBICALL_DATA` | Installed CBIcall resource bundle | `/cbicall-data` | Native CBIcall WES/WGS/mtDNA databases and external tools |

For nf-core workflows, Nextflow may also use a Singularity/Apptainer cache for
task containers. That cache is separate from `SIF_IMAGE`; it stores containers
used by nf-core pipeline tasks, not the CBIcall application image.

```bash
# Choose where you want to keep the reusable CBIcall Apptainer image.
export SIF_IMAGE=/path/you/choose/cbicall_latest.sif

# Per-user writable CBIcall copy. This path is bind-mounted over /usr/share/cbicall.
export CBICALL_WRITABLE=$HOME/cbicall

# Native CBIcall resource bundle. Not required for nf-core-only runs.
export CBICALL_DATA=/path/to/cbicall-data

# Optional: Nextflow task-container cache for nf-core workflows.
export NXF_SINGULARITY_CACHEDIR=/path/to/nextflow-singularity-cache
export NXF_SINGULARITY_LIBRARYDIR=/path/to/nextflow-singularity-cache
```

`SIF_IMAGE` can be stored wherever you pull it. `CBICALL_WRITABLE` is commonly
per-user because it contains editable configuration. `CBICALL_DATA` is often
shared by an administrator for native workflows.

## Installing the Container Image

### Loading Apptainer

If your HPC system uses environment modules, load Apptainer first:

```bash
module load apptainer
```

> **Note:** On some HPC systems, `apptainer` is already available by default and does not need to be loaded explicitly.

---

### Obtaining the CBIcall Container Image

There are two ways to run the CBIcall container with Apptainer.  
For most users and all production workflows, **explicitly pulling the image is recommended**.

---

#### Option A: Pull the image explicitly (`apptainer pull`) — **recommended**

Pull the latest CBIcall image from Docker Hub and save it as a Singularity Image File (`.sif`):

```bash
apptainer pull "$SIF_IMAGE" docker://manuelrueda/cbicall:latest
```

This command:

- downloads the image once
- creates a reusable `.sif` file
- avoids reliance on the user cache
- is suitable for batch scheduling and offline execution

This step only needs to be performed once.

---

#### Option B: Run directly from Docker Hub (`apptainer run`)

Apptainer can execute the CBIcall container directly from Docker Hub without explicitly pulling the image first:

```bash
apptainer run docker://manuelrueda/cbicall:latest
```

On first use, the image is automatically downloaded, converted to a Singularity Image File (`.sif`), and cached locally. Subsequent runs reuse the cached image.

This method does **not** create a persistent `.sif` file in the working directory. The cached image is stored in Apptainer’s user cache directory.

This option is suitable for:
- quick testing
- exploratory or interactive use

For reproducible workflows, offline execution, or batch jobs, explicitly pulling the image (Option A) is strongly recommended.

---

### Preparing a Writable CBIcall Directory

The CBIcall installation inside the container (`/usr/share/cbicall`) is **read-only**.  
To allow configuration changes and execution of tests, create a writable copy in your home directory:

```bash
apptainer exec "$SIF_IMAGE" \
  bash -lc 'mkdir -p "$CBICALL_WRITABLE" && cp -a /usr/share/cbicall/. "$CBICALL_WRITABLE"/'
```

This step only needs to be done once.

---

## Choose a Workflow Path

CBIcall can be used with registered nf-core workflows before downloading the
large CBIcall germline resource bundle.

| Workflow path | Resource bundle required? | Runtime requirement |
| --- | --- | --- |
| `workflow_provider: nf-core` | No | Nextflow plus Singularity/Apptainer, usually provided by the HPC environment. |
| Native CBIcall Bash/Snakemake/Nextflow WES/WGS/mtDNA | Yes | Bind the CBIcall resource bundle and configure `DATADIR`. |

Use the resource download steps below for native CBIcall workflows. For nf-core
provider runs, start with the checked-in `examples/input/nf-core-demo.yaml` and
the nf-core Provider page in the online documentation.

---

## Download the Resource Bundle for Native Workflows

> **Note:** this process can be lengthy.

Begin by downloading the required databases and software for native CBIcall
workflows. Save the data **outside the container**; this preserves it across
container runs and avoids repeated downloads. Skip this section for nf-core
provider runs that use Nextflow/nf-core-managed resources.

Install the Python dependency (on the host):

```bash
pip3 install --user gdown
```

Choose a location on the host filesystem where the databases will be stored:

```bash
mkdir -p "$CBICALL_DATA"
cd "$CBICALL_DATA"
```

Download the data preparation script and execute it:

```bash
wget https://raw.githubusercontent.com/mrueda/cbicall/refs/heads/main/scripts/download_cbicall_bundle.py
python3 ./download_cbicall_bundle.py --outdir "$CBICALL_DATA"
```

To verify only the catalog-to-Google-Drive bundle identity before starting the large archive download:

```bash
python3 ./download_cbicall_bundle.py \
  --outdir "$CBICALL_DATA" \
  --verify-resource-id-only
```

Google Drive may occasionally block or stall automated downloads. If this happens, print the manual download list:

```bash
python3 ./download_cbicall_bundle.py \
  --outdir "$CBICALL_DATA" \
  --print-manual-download
```

Download every listed file into `$CBICALL_DATA`, then let the script continue from those files:

```bash
python3 ./download_cbicall_bundle.py \
  --outdir "$CBICALL_DATA" \
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

---

## Running and Interacting with the Container

Start an interactive shell inside the container, overlaying the writable CBIcall
copy. For native CBIcall workflows, also bind the external data directory:

```bash
apptainer shell \
  --pwd /usr/share/cbicall \
  --bind "$CBICALL_WRITABLE":/usr/share/cbicall \
  --bind "$CBICALL_DATA":/cbicall-data \
  "$SIF_IMAGE"
```

You will start directly in the CBIcall working directory.

For nf-core-only validation, the CBIcall data bind is not required:

```bash
apptainer shell \
  --pwd /usr/share/cbicall \
  --bind "$CBICALL_WRITABLE":/usr/share/cbicall \
  "$SIF_IMAGE"
```

### Point Native Workflows to `/cbicall-data`

Native CBIcall workflows read resource paths from Bash `env.sh` files and from
Snakemake/Nextflow `config.yaml` files. In Apptainer, bind your external
resource directory as
`/cbicall-data` and point CBIcall to that container path:

```bash
sed -i 's|^DATADIR=.*|DATADIR=/cbicall-data|' workflows/bash/gatk-4.6/env.sh
sed -i 's|^DATADIR=.*|DATADIR=/cbicall-data|' workflows/bash/gatk-3.5/env.sh
sed -i 's|^datadir:.*|datadir: "/cbicall-data"|' workflows/snakemake/gatk-4.6/config.yaml
```

The native Nextflow config is a symlink to this shared GATK 4.6 backend config, so one edit updates both Snakemake and Nextflow native workflows.

Confirm that CBIcall sees the mounted resources:

```bash
bin/cbicall validate-resources
bin/cbicall validate-parameters -p examples/input/param.yaml
```

---

## Performing Integration Tests

Inside the container, from the CBIcall repository root:

### WES

```bash
bin/cbicall test --wes-bash -t 1
```

### mtDNA

```bash
bin/cbicall test --mit-bash -t 1
```

---

## Notes

- `$HOME` inside the container corresponds to your **host home directory**.
- All configuration changes are performed on the writable copy in `CBICALL_WRITABLE`.
- `SIF_IMAGE` is only the CBIcall application image. It is independent from
  `NXF_SINGULARITY_CACHEDIR` and `NXF_SINGULARITY_LIBRARYDIR`, which are used by
  Nextflow for nf-core task containers.
- The container image itself remains immutable, ensuring reproducibility and HPC safety.

## How to run a job in **Slurm**

You can find an example here:

[Slurm Job Example](https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/refs/heads/main/examples/scripts/run_cbicall_apptainer_slurm.sh)

## Cleaning Up (Apptainer / Singularity)

Unlike Docker, Apptainer does not manage running containers or image caches.
A container image is simply a regular file on disk.

### Removing the Container Image

To delete a CBIcall container image (`.sif` file), simply remove it:

```bash
cd  # home directory
rm "$SIF_IMAGE"
```

If the image is currently in use by a running job, the removal will fail until the job finishes.

---

### Removing the Writable CBIcall Copy (Optional)

If you created a writable copy of the CBIcall installation and want to reset it:

```bash
rm -rf "$CBICALL_WRITABLE"
```

This does **not** affect the container image itself.

---

### Removing Downloaded Databases (Optional)

To remove the externally downloaded CBIcall databases:

```bash
rm -rf $CBICALL_DATA
```

⚠️ This operation is irreversible and will require re-downloading the data.

---

### Summary

- Apptainer images (`.sif`) are removed with `rm`
- There is no equivalent to `docker stop`, `docker rm`, or `docker rmi`
- User data and configuration live outside the container and must be cleaned separately
