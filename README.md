<div align="center">
  <a href="https://github.com/mrueda/cbicall">
    <img src="https://raw.githubusercontent.com/mrueda/cbicall/main/docs/img/cbicall-logo.png"
         width="300" alt="CBIcall">
  </a>
  <p><em>CNAG Biomedical Informatics framework for variant calling</em></p>
</div>


[![Build](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/build-and-test.yml)
![version](https://img.shields.io/badge/version-1.0.0-28a745)
[![Coverage Status](https://coveralls.io/repos/github/CNAG-Biomedical-Informatics/cbicall/badge.svg?branch=main)](https://coveralls.io/github/CNAG-Biomedical-Informatics/cbicall?branch=main)
[![Docker Build](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/docker-build-multi-arch.yml/badge.svg?branch=main)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/docker-build-multi-arch.yml)
[![Docker Image Size](https://img.shields.io/docker/image-size/manuelrueda/cbicall/latest?logo=docker&label=image%20size)](https://hub.docker.com/r/manuelrueda/cbicall/)
[![Docker Pulls](https://img.shields.io/docker/pulls/manuelrueda/cbicall.svg?logo=docker&label=pulls)](https://hub.docker.com/r/manuelrueda/cbicall/)
[![Documentation Status](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/documentation.yml/badge.svg)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/documentation.yml)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

**CBIcall** (**C**NAG **B**iomedical **I**nformatics framework for variant **call**ing) is a
configuration-driven workflow dispatcher for Illumina sequencing variant calling,
designed to support reproducible, auditable, and consistent execution of analysis pipelines
in production environments.

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
parameters, resolves workflows from a versioned registry, dispatches curated Bash
or Snakemake pipelines, and captures logs and structured metadata for traceability.

Key points:

- Configuration-driven execution from a YAML parameter file
- Curated Bash and Snakemake workflow support
- Support for WES, WGS, and mtDNA analysis modes
- Structured run logging and traceable runtime metadata
- Optional partial workflow starts for supported engines

# Quick Start

    bin/cbicall -p params.yaml -t 8

Runnable examples and sample inputs are available under `examples/`.

# Documentation

The full technical reference lives in the documentation site and repository docs:

- General usage and parameter reference: [docs/usage/usage.md](docs/usage/usage.md)
- Quick start: [docs/usage/quickstart.md](docs/usage/quickstart.md)
- End-to-end examples: [docs/usage/end-to-end-example-wes.md](docs/usage/end-to-end-example-wes.md), [docs/usage/end-to-end-example-mit.md](docs/usage/end-to-end-example-mit.md)
- Technical details: [docs/technical-details/architecture.md](docs/technical-details/architecture.md)

# Citation

CBIcall: a configuration-driven framework for variant calling in large sequencing cohorts. _Manuscript In preparation._

# Author

Written by Manuel Rueda (mrueda). GitHub repository: [https://github.com/CNAG-Biomedical-Informatics/cbicall](https://github.com/CNAG-Biomedical-Informatics/cbicall).

# Copyright and license

Please see the included LICENSE file for distribution and usage terms.

