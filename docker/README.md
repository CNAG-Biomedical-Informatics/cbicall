# Docker installation

The Docker image contains CBIcall and all four supported workflow backends. The
large reference and software bundle remains outside the image and is mounted
only for bundled `cbicall-core` workflows.

## 1. Get the image

For reproducible analyses, use a released version rather than `latest`:

```bash
export CBICALL_IMAGE=manuelrueda/cbicall:1.2.1
docker pull "$CBICALL_IMAGE"
```

<details>
<summary>Build the checked-out source instead</summary>

```bash
git clone https://github.com/CNAG-Biomedical-Informatics/cbicall.git
cd cbicall
git checkout v1.2.1

docker buildx build --load \
  -f docker/Dockerfile \
  --build-arg CBICALL_VERSION="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' src/cbicall/__about__.py)" \
  --build-arg VCS_REF="$(git rev-parse HEAD)" \
  -t cbicall:local .
```

The Docker build uses the current checkout as its build context; it does not
clone a different revision. Inspect the recorded commit with:

```bash
docker inspect cbicall:local \
  --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

</details>

## 2. Choose the workflow source

| Workflow source | CBIcall bundle required? | Notes |
| --- | --- | --- |
| `workflow_provider: cbicall-core` (default) | Yes | Bundled native WES/WGS workflows support Bash, Snakemake, Nextflow, and Cromwell. mtDNA uses Bash and requires x86_64. |
| `workflow_provider: nf-core` | No | The selected nf-core task runtime must be available from the container setup. Running CBIcall directly on the host is often simpler for Docker-based nf-core runs. |

## 3. Install resources for `cbicall-core`

Keep the bundle in a persistent host directory:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
mkdir -p "$CBICALL_DATA"

docker run --rm \
  --volume "$CBICALL_DATA":/cbicall-data \
  "$CBICALL_IMAGE" \
  cbicall install-resources --outdir /cbicall-data
```

Confirm that the image can see the installation:

```bash
docker run --rm \
  --volume "$CBICALL_DATA":/cbicall-data \
  --env CBICALL_DATA=/cbicall-data \
  "$CBICALL_IMAGE" \
  cbicall doctor
```

Skip this step when using only an external nf-core workflow.

<details>
<summary>Manual resource-download recovery</summary>

If the automatic Google Drive download fails, print the registered file list:

```bash
docker run --rm \
  --volume "$CBICALL_DATA":/cbicall-data \
  "$CBICALL_IMAGE" \
  cbicall install-resources --outdir /cbicall-data --print-manual-download
```

Place every listed file in `$CBICALL_DATA`, then resume:

```bash
docker run --rm \
  --volume "$CBICALL_DATA":/cbicall-data \
  "$CBICALL_IMAGE" \
  cbicall install-resources --outdir /cbicall-data --skip-download
```

Use `--verify-resource-id-only` to verify the small catalog-pinned identifier
before downloading the archive. Add `--remove-parts` after successful
verification when disk space is limited.

</details>

## 4. Run an analysis

Keep the parameters YAML and input data under one host project directory. Bind
that directory at the same absolute path so paths in the YAML remain valid:

```bash
export PROJECT_DIR=/absolute/path/to/project

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --env CBICALL_DATA=/cbicall-data \
  --volume "$CBICALL_DATA":/cbicall-data \
  --volume "$PROJECT_DIR":"$PROJECT_DIR" \
  --workdir "$PROJECT_DIR" \
  "$CBICALL_IMAGE" \
  cbicall run -p "$PROJECT_DIR/parameters.yaml" -t 4
```

The run directory is written to the mounted project directory and remains on
the host after the container exits.

## 5. Run an integration test

Use a new or empty host directory to retain the generated reports:

```bash
export TEST_DIR=/absolute/path/to/cbicall-wes-test
mkdir -p "$TEST_DIR"

docker run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp \
  --env CBICALL_DATA=/cbicall-data \
  --volume "$CBICALL_DATA":/cbicall-data \
  --volume "$TEST_DIR":/work \
  "$CBICALL_IMAGE" \
  cbicall test --wes-bash -t 1 --workspace /work
```

## System requirements

- Docker with Buildx support for local multi-platform builds.
- A Linux amd64 or arm64 host; mtDNA workflows require amd64/x86_64.
- Sufficient CPU, memory, and storage for the selected workflow. Allow at least
  16 GB RAM and 100 GB disk for the bundle and shipped WES test.
