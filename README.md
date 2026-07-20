<div align="center">
  <a href="https://github.com/CNAG-Biomedical-Informatics/cbicall">
    <img src="https://raw.githubusercontent.com/CNAG-Biomedical-Informatics/cbicall/main/docs-site/static/img/cbicall-logo.png"
         width="300" alt="CBIcall">
  </a>
  <p><em>CNAG Biomedical Informatics framework for variant calling</em></p>
</div>


[![Build](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/build-and-test.yml)
[![PyPI](https://img.shields.io/pypi/v/cbicall.svg)](https://pypi.org/project/cbicall/)
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

**Documentation:** <https://cnag-biomedical-informatics.github.io/cbicall/>


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

# Installation

Install the released command from PyPI (recommended):

```bash
python3 -m pip install --upgrade cbicall
```

Install all optional Python integrations when needed:

```bash
python3 -m pip install --upgrade "cbicall[all]"
```

Use an editable source checkout only for development:

```bash
git clone https://github.com/CNAG-Biomedical-Informatics/cbicall.git
cd cbicall
python3 -m pip install -e ".[all,test]"
```

Native workflows also require the separately distributed resource bundle. See
the documentation for Python, Docker, Apptainer, cloud, and HPC installation
instructions.

Point all native backends to an installed bundle with:

```bash
export CBICALL_DATA=/absolute/path/to/cbicall-data
cbicall doctor
```

# Quick Start

Generate a self-contained tour of the WES audit report and interactive mtDNA
browser without installing the external resource bundle:

```bash
cbicall demo
```

This uses precomputed outputs from the packaged integration fixture. It does
not execute BWA, GATK, or MToolBox. For a real analysis, use a parameters YAML:

```bash
cbicall run -p params.yaml -t 8
```

Runnable examples and sample inputs are available under `examples/`.

# Citation

CBIcall: a configuration-driven framework for variant calling in large sequencing cohorts. [Preprint DOI](https://doi.org/10.64898/2026.03.23.713646).

# Author

Written by Manuel Rueda (mrueda). GitHub repository: [https://github.com/CNAG-Biomedical-Informatics/cbicall](https://github.com/CNAG-Biomedical-Informatics/cbicall).

# Copyright and license

Please see the [GPLv3 license](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/LICENSE) for distribution and usage terms.
