# HPC installation with Apptainer

Use Apptainer when Docker is unavailable on an HPC system. The CBIcall image is
immutable; project files and the optional `cbicall-core` resource bundle remain
on the host and are bind-mounted at runtime.

## 1. Pull a released image

Load Apptainer if your site provides it through environment modules, then pull
a versioned image:

```bash
module load apptainer 2>/dev/null || true

export SIF_IMAGE=/absolute/path/to/cbicall_1.2.0.sif
apptainer pull "$SIF_IMAGE" docker://manuelrueda/cbicall:1.2.0
```

Keep the `.sif` file for later interactive and scheduled runs. Pinning the image
tag avoids an unnoticed change from `latest`.

## 2. Choose the workflow source

| Workflow source | CBIcall bundle required? | Recommended HPC setup |
| --- | --- | --- |
| `workflow_provider: cbicall-core` (default) | Yes | Run the CBIcall SIF and bind the resource and project directories as shown below. |
| `workflow_provider: nf-core` | No | Install CBIcall on the host, then use the site's Nextflow and Apptainer modules so nf-core can launch its task containers normally. |

The remaining commands on this page cover bundled `cbicall-core` workflows.

## 3. Install the resource bundle

Choose a persistent host directory:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
mkdir -p "$CBICALL_DATA"

apptainer exec \
  --bind "$CBICALL_DATA":/cbicall-data \
  "$SIF_IMAGE" \
  cbicall install-resources --outdir /cbicall-data
```

Verify the mounted installation:

```bash
apptainer exec \
  --bind "$CBICALL_DATA":/cbicall-data \
  --env CBICALL_DATA=/cbicall-data \
  "$SIF_IMAGE" \
  cbicall doctor
```

<details>
<summary>Manual resource-download recovery</summary>

If the automatic Google Drive download fails, print the registered file list:

```bash
apptainer exec \
  --bind "$CBICALL_DATA":/cbicall-data \
  "$SIF_IMAGE" \
  cbicall install-resources --outdir /cbicall-data --print-manual-download
```

Place every listed file in `$CBICALL_DATA`, then resume:

```bash
apptainer exec \
  --bind "$CBICALL_DATA":/cbicall-data \
  "$SIF_IMAGE" \
  cbicall install-resources --outdir /cbicall-data --skip-download
```

Use `--verify-resource-id-only` to verify the small catalog-pinned identifier
before downloading the archive. Add `--remove-parts` after successful
verification when disk space is limited.

</details>

## 4. Run an analysis

Keep the parameters YAML and input data under one project directory. Bind that
directory at the same absolute path so paths in the YAML remain valid:

```bash
export PROJECT_DIR=/absolute/path/to/project

apptainer exec \
  --bind "$CBICALL_DATA":/cbicall-data \
  --bind "$PROJECT_DIR":"$PROJECT_DIR" \
  --env CBICALL_DATA=/cbicall-data \
  --pwd "$PROJECT_DIR" \
  "$SIF_IMAGE" \
  cbicall run -p "$PROJECT_DIR/parameters.yaml" -t 4
```

The run directory is written directly to the host project directory. Normal
execution does not require a writable copy of `/usr/share/cbicall`.

## 5. Run an integration test

Use a new or empty host directory to retain the test outputs:

```bash
export TEST_DIR=/absolute/path/to/cbicall-wes-test
mkdir -p "$TEST_DIR"

apptainer exec \
  --bind "$CBICALL_DATA":/cbicall-data \
  --bind "$TEST_DIR":"$TEST_DIR" \
  --env CBICALL_DATA=/cbicall-data \
  "$SIF_IMAGE" \
  cbicall test --wes-bash -t 1 --workspace "$TEST_DIR"
```

## Slurm

Use the same binds inside the scheduled job. A complete template is available
in [run_cbicall_apptainer_slurm.sh](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/examples/scripts/run_cbicall_apptainer_slurm.sh).

<details>
<summary>Site-specific workflow development</summary>

The packaged workflow definitions are intentionally immutable. Developers who
need to modify a workflow, registry, or Bash runtime profile should use a
version-matched source checkout and bind that checkout explicitly. Normal users
need only the SIF, project directory, and resource bundle.

</details>

## Notes

- The host home directory is normally visible inside Apptainer unless restricted by site policy.
- The CBIcall SIF is separate from any Nextflow task-container cache used by external nf-core workflows.
- The bundled MToolBox mtDNA workflow runs on x86_64 only.
