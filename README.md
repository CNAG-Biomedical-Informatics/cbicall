<div align="center">
  <a href="https://github.com/mrueda/cbicall">
    <img src="https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/main/docs-site/static/img/cbicall-logo.png"
         width="300" alt="CBIcall">
  </a>
  <p><em>CNAG Biomedical Informatics framework for variant calling</em></p>
</div>


[![Build](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/build-and-test.yml)
![version](https://img.shields.io/badge/version-1.0.1--beta.1-28a745)
[![Coverage Status](https://coveralls.io/repos/github/CNAG-Biomedical-Informatics/cbicall/badge.svg?branch=main)](https://coveralls.io/github/CNAG-Biomedical-Informatics/cbicall?branch=main)
[![Docker Build](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/docker-build-multi-arch.yml/badge.svg?branch=main)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/docker-build-multi-arch.yml)
[![Docker Image Size](https://img.shields.io/docker/image-size/manuelrueda/cbicall/latest?logo=docker&label=image%20size)](https://hub.docker.com/r/manuelrueda/cbicall/)
[![Docker Pulls](https://img.shields.io/docker/pulls/manuelrueda/cbicall.svg?logo=docker&label=pulls)](https://hub.docker.com/r/manuelrueda/cbicall/)
[![Documentation Status](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/documentation.yml/badge.svg)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/documentation.yml)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

**CBIcall** (**C**NAG **B**iomedical **I**nformatics framework for variant **call**ing) is a
**configuration-driven framework** for reproducible variant calling in large sequencing cohorts.
It validates analysis intent, resolves approved workflow backends and resource definitions,
and records **structured run reports** for auditable production runs and run-to-run comparison.

**📘 Documentation:** <a href="https://cnag-biomedical-informatics.github.io/cbicall" target="_blank">https://cnag-biomedical-informatics.github.io/cbicall</a>

**🐳 Docker Hub Image:** <a href="https://hub.docker.com/r/manuelrueda/cbicall/tags" target="_blank">https://hub.docker.com/r/manuelrueda/cbicall/tags</a>


# Table of contents
- [Installation](#installation)
  - [Non-Containerized](non-containerized/README.md)
  - Containerized
    - [Docker](docker/README.md)
    - [Apptainer](apptainer/README.md)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Citation](#citation)
  - [Author](#author)
- [License](#copyright-and-license)


CBIcall orchestrates germline variant calling workflows for Illumina sequencing data.
It does **not** implement variant calling algorithms itself. Instead, it validates
parameters, resolves workflows from a versioned registry, checks resource compatibility,
launches native CBIcall workflows through **Bash, Snakemake, Nextflow, or Cromwell**
backends, and captures audit artifacts for traceability. Selected external nf-core
workflows can also run through the same validation and provenance layer.

CBIcall uses a three-part execution contract:

| Contract layer | Purpose |
| --- | --- |
| Parameters YAML | User analysis intent: inputs, pipeline, mode, genome, backend, and runtime options. |
| Workflow registry | Approved workflow implementations and backend-specific entrypoints. |
| Resource catalog | External references, tool/resource bundles, compatibility rules, and resource identity. |

Key points:

- Configuration-driven execution from a YAML parameter file
- Native CBIcall workflow support through Bash, Snakemake, Nextflow, and Cromwell backends
- Support for WES, WGS, and mtDNA analysis modes
- Registry-backed support for selected external nf-core/Nextflow workflows
- Structured audit artifacts: `log.json`, `run-report.json`, optional `run-report.html`, workflow fingerprints, resource identity, output inventories, and normalized VCF hashes
- Programmatic run comparison with `cbicall compare-runs`
- Optional partial workflow starts for supported backends

Workflow sources:

| Source | Role |
| --- | --- |
| Native CBIcall workflows | Packaged WES/WGS/mtDNA pipelines with CBIcall validation, logging, and output structure. |
| External nf-core workflows | Selected registry-backed Nextflow workflows executed with CBIcall validation and provenance. |

# Quick Start

    bin/cbicall run -p params.yaml -t 8

Runnable examples and sample inputs are available under `examples/`.

# Documentation

The full technical reference lives in the documentation site and repository docs:

- General usage and parameter reference: [docs-site/docs/usage/usage.md](docs-site/docs/usage/usage.md)
- Quick start: [docs-site/docs/usage/quickstart.md](docs-site/docs/usage/quickstart.md)
- End-to-end examples: [WES/WGS](docs-site/docs/usage/end-to-end-example-wes.md), [mtDNA](docs-site/docs/usage/end-to-end-example-mit.md)
- Run comparison and audit reports: [docs-site/docs/usage/run-comparison.md](docs-site/docs/usage/run-comparison.md)
- Technical details: [docs-site/docs/technical-details/architecture.md](docs-site/docs/technical-details/architecture.md)

# Citation

CBIcall: a configuration-driven framework for variant calling in large sequencing cohorts. [Preprint DOI](https://doi.org/10.64898/2026.03.23.713646).

# Author

Written by Manuel Rueda (mrueda). GitHub repository: [https://github.com/CNAG-Biomedical-Informatics/cbicall](https://github.com/CNAG-Biomedical-Informatics/cbicall).

# Copyright and license

Please see the included LICENSE file for distribution and usage terms.

