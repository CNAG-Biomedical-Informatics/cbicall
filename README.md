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
[![Documentation Status](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/documentation.yml/badge.svg)](https://github.com/cnag-biomedical-informatics/cbicall/actions/workflows/documentation.yml)
[![License](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

**CBIcall** (**C**NAG **B**iomedical **I**nformatics framework for variant **call**ing) is a
configuration-driven framework for running and auditing variant-calling
workflows across heterogeneous computing environments.

Users describe an analysis in one YAML file. CBIcall validates the request,
resolves compatible workflows and resources, launches the selected backend, and
records structured evidence for reproducibility and run comparison. The bundled
workflow collection is identified as `cbicall-core`; selected external nf-core
workflows are also supported.

**Documentation:** <https://cnag-biomedical-informatics.github.io/cbicall/>

## Installation

Install CBIcall from PyPI:

```bash
python3 -m pip install --upgrade cbicall
```

Ready-to-run workflows require additional tools and reference resources. The
documentation covers the resource bundle, optional Python integrations, source
installation, Docker, and Apptainer.

## Quick Start

Generate example WES and mtDNA reports without installing workflow dependencies
or the external resource bundle:

```bash
cbicall demo
```

For a configured analysis:

```bash
cbicall run -p parameters.yaml -t 4
```

## Citation

Rueda M, Fernandez-Orth D, Gut IG. CBIcall: a configuration-driven framework for variant calling in large sequencing cohorts. *Bioinformatics Advances*. 2026. [doi:10.1093/bioadv/vbag232](https://doi.org/10.1093/bioadv/vbag232).

## Author

Manuel Rueda, PhD. CNAG: <https://www.cnag.eu>

## License

CBIcall is distributed under the
[GPLv3 license](https://github.com/CNAG-Biomedical-Informatics/cbicall/blob/main/LICENSE).
