# Containerized Installation (HPC with Apptainer / Singularity)

This section describes how to install and run **CBIcall** on High-Performance Computing (HPC) systems where Docker is not available, using **Apptainer (formerly Singularity)**.

Apptainer can execute Docker images directly, but container images are **read-only by design**. As a result, databases, configuration files, and workflows must be stored outside the container and bind-mounted at runtime.

---

## Downloading Required Databases and Software

> **Note:** this process can be lengthy.

Begin by downloading the required databases and software. Save the data **outside the container**; this preserves it across container runs and avoids repeated downloads.

Install the Python dependency (on the host):

```bash
pip3 install --user gdown
```

Choose a location on the host filesystem where the databases will be stored:

```bash
CBICALL_DATA=/software/biomed/cbicall-data
mkdir -p "$CBICALL_DATA"
cd "$CBICALL_DATA"
```

Download the data preparation script and execute it:

```bash
wget https://raw.githubusercontent.com/mrueda/cbicall/refs/heads/main/scripts/01_download_external_data.py
python3 ./01_download_external_data.py --outdir "$CBICALL_DATA"
```

To verify only the catalog-to-Google-Drive bundle identity before starting the large archive download:

```bash
python3 ./01_download_external_data.py \
  --outdir "$CBICALL_DATA" \
  --verify-bundle-id-only
```

Google Drive may occasionally block or stall automated downloads. If this happens, print the manual download list:

```bash
python3 ./01_download_external_data.py \
  --outdir "$CBICALL_DATA" \
  --print-manual-download
```

Download every listed file into `$CBICALL_DATA`, then let the script continue from those files:

```bash
python3 ./01_download_external_data.py \
  --outdir "$CBICALL_DATA" \
  --skip-download
```

The script will:

- download missing split files when possible
- reassemble `data.tar.gz`
- verify the split parts or assembled archive with `data.tar.gz.md5`
- load the CBIcall resource catalog, locally or from the catalog URL
- optionally verify a small GDrive bundle identifier file such as `cbicall-bundle-id.json`
- rename the verified archive using the bundle identity, for example `cbicall-germline-resources-v1.tar.gz`
- extract the archive into `DATADIR`
- write `cbicall-resource-installation.json` with the installed bundle provenance

If disk space is tight and the checksum has passed, add `--remove-parts` to remove `data.tar.gz.part-*` after assembly.

CBIcall keeps the rich resource registry in `resources/cbicall-resource-catalog.json`. The GDrive bundle only needs a small identifier file, for example `cbicall-bundle-id.json` containing `{"bundle": "cbicall-germline-resources-v1"}`. When that identifier file is available, the registry can store its Google Drive file ID and SHA-256 so the downloader can verify that the remote bundle matches the local CBIcall catalog entry.

---

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
apptainer pull cbicall_latest.sif docker://manuelrueda/cbicall:latest
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
apptainer exec cbicall_latest.sif \
  bash -lc 'mkdir -p $HOME/cbicall && cp -a /usr/share/cbicall/. $HOME/cbicall/'
```

This step only needs to be done once.

---

### Running and Interacting with the Container

Start an interactive shell inside the container, overlaying the writable CBIcall copy and binding the external data directory:

```bash
apptainer shell \
  --pwd /usr/share/cbicall \
  --bind "$HOME/cbicall":/usr/share/cbicall \
  --bind "$CBICALL_DATA":/cbicall-data \
  cbicall_latest.sif
```

You will start directly in the CBIcall working directory.

### Point CBIcall to `/cbicall-data`

CBIcall workflows read resource paths from Bash `env.sh` files and the
Snakemake `config.yaml`. In Apptainer, bind your external resource directory as
`/cbicall-data` and point CBIcall to that container path:

```bash
sed -i 's|^DATADIR=.*|DATADIR=/cbicall-data|' workflows/bash/gatk-4.6/env.sh
sed -i 's|^DATADIR=.*|DATADIR=/cbicall-data|' workflows/bash/gatk-3.5/env.sh
sed -i 's|^datadir:.*|datadir: "/cbicall-data"|' workflows/snakemake/gatk-4.6/config.yaml
```

Confirm that CBIcall sees the mounted resources:

```bash
bin/cbicall validate-resources
bin/cbicall doctor -p examples/input/param.yaml
```

---

## Performing Integration Tests

Inside the container, from the CBIcall repository root:

### WES

```bash
bin/cbicall test --wes -t 1
```

### mtDNA

```bash
bin/cbicall test --mit -t 1
```

---

## Notes

- `$HOME` inside the container corresponds to your **host home directory**.
- All configuration changes are performed on the writable copy in `$HOME/cbicall`.
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
rm cbicall_latest.sif
```

If the image is currently in use by a running job, the removal will fail until the job finishes.

---

### Removing the Writable CBIcall Copy (Optional)

If you created a writable copy of the CBIcall installation in your home directory and want to reset it:

```bash
rm -rf $HOME/cbicall
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
