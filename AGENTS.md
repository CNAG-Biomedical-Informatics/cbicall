# Repository maintenance notes

## Backend-local configuration copies

Do not replace the stack-local environment files with symbolic links. Python
wheels and other installation formats do not preserve symlink behavior
reliably, so CBIcall intentionally ships regular files at both locations.

The following pairs must remain byte-for-byte identical:

- `workflows/bash/gatk-3.5/env.sh` and `workflows/bash/gatk-4.6/env.sh`
- `workflows/bash/gatk-3.5/cnag-hpc-env.sh` and `workflows/bash/gatk-4.6/cnag-hpc-env.sh`

Treat the GATK 4.6 file as the editing source, then copy it to GATK 3.5. For
the CNAG profile, use:

```bash
cp workflows/bash/gatk-4.6/cnag-hpc-env.sh \
   workflows/bash/gatk-3.5/cnag-hpc-env.sh
```

The Snakemake, Nextflow, and Cromwell `gatk-4.6/config.yaml` files follow the
same packaging rule: they remain registry-local regular files and must stay
byte-for-byte synchronized. The test suite enforces these invariants.

The `cnag-hpc` runtime profile is deliberately supported only by the Bash
backend.
