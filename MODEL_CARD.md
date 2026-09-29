# CBIcall

> **Research use only. Not a medical device.**
> This software was developed in the HEREDITARY project as a research prototype. It has not been assessed or certified under the EU Medical Device Regulation (EU) 2017/745. It must not be used to diagnose, treat, monitor or make any other decision about individual patients.

## 1. Overview

| Field | Value |
|---|---|
| Name | CBIcall |
| Version | 1.3.0 |
| Type | Software or library |
| Owner partner and contact | CNAG / Manuel Rueda / [manuel.rueda@cnag.eu](mailto:manuel.rueda@cnag.eu) |
| License | GPL-3.0-or-later |
| Repository and DOI | [GitHub repository](https://github.com/CNAG-Biomedical-Informatics/cbicall); [doi:10.1093/bioadv/vbag232](https://doi.org/10.1093/bioadv/vbag232) |
| Release date | 2026-08-22 |
| Related task or deliverable | WP 3.2 |
| How to cite | Rueda M, Fernandez-Orth D, Gut IG. CBIcall: a configuration-driven framework for variant calling in large sequencing cohorts. *Bioinformatics Advances*. 2026;6(1):vbag232. [doi:10.1093/bioadv/vbag232](https://doi.org/10.1093/bioadv/vbag232). |

## 2. Intended use

**Primary intended use.** CBIcall is a configuration-driven framework for reproducible germline variant calling. It validates configurations, runs WES, WGS and mitochondrial DNA workflows from paired-end FASTQ data, and records provenance for analysis-ready VCFs. Within HEREDITARY, it supports standardized processing of human sequencing cohorts.

**Intended users.** Biomedical researchers, bioinformaticians and genomics facilities authorized to process the input data.

**Out-of-scope uses.** This release must not be used for:

- clinical diagnosis, prognosis, treatment, monitoring or other decisions about individuals;
- re-identification of individuals or unauthorized linkage with other data sources;
- processing human genomic data without an appropriate consent, authorization or other valid legal basis;
- unsupported data types, reference assemblies or configurations; or
- interpreting variant calls as clinically significant without independent validation and qualified expert review.

## 3. Data provenance

CBIcall contains no trained model. The following external datasets were used for benchmarking, reproducibility testing or cohort-scale demonstration; they are not distributed with CBIcall.

| Dataset | Source partner | Type (public, synthetic, pseudonymized, anonymized) | Population (size, age, sex, country, condition) | Ethics approval or access conditions |
|---|---|---|---|---|
| GIAB Ashkenazim trio, v4.2.1 (GRCh38) | NIST Genome in a Bottle Consortium | Public benchmark data | HG002, HG003 and HG004; Ashkenazi Jewish family trio; both sexes | Public benchmark data; evaluation restricted to GIAB high-confidence regions. |
| 1000 Genomes Project Phase 3 | 1000 Genomes Project / IGSR | Public research data | HG00103 for reproducibility testing; 503 European-ancestry participants, both sexes, as a reference cohort | Publicly available under 1000 Genomes/IGSR data-use conditions. |
| NINDS Exome Sequencing in Parkinson's Disease, dbGaP phs001172.v1.p2 | NINDS / dbGaP | Controlled-access, pseudonymized research data | 618 North American Parkinson's disease WES samples processed; 599 self-reported North American Caucasian participants retained after quality control | dbGaP controlled access; Data Access Request 142786-2 and applicable consent restrictions. |

The repository and released files contain no personal data, pseudonymized records, identifiers or pseudonymization keys. Users must obtain source datasets separately and follow their access conditions.

## 4. Method

- **Algorithm or architecture.** CBIcall is a Python validation, execution and audit layer, not a new variant caller. It validates YAML configurations against workflow and resource definitions, then runs approved WES, WGS or mitochondrial DNA workflows through Bash, Cromwell, Nextflow or Snakemake. It records configuration, software, resources and outputs for reproducibility.
- **Training setting.** Not applicable. CBIcall contains no trained statistical or machine-learning model and performs no model training.
- **Privacy techniques.** None. CBIcall runs in the user's computing environment and does not require transfer to an external service. Users remain responsible for access control and data governance.

## 5. Performance

Performance here refers to benchmarking of workflows executed through CBIcall, not to predictive-model performance.

- **Evaluation data.** Accuracy was benchmarked against GIAB v4.2.1 truth sets for the HG002/HG003/HG004 WGS trio at approximately 30× coverage. Reproducibility used 1000 Genomes WES sample HG00103. Cohort-scale operation was demonstrated on 1,102 WES samples.
- **Metrics and results.** GIAB SNV F1 scores were 99.12%–99.15% for all calls and 94.91%–94.93% for PASS calls; indel F1 scores were 98.52%–98.60% and 95.39%–95.71%, respectively. Four computing environments produced identical call-level fingerprints for the reproducibility sample. Mitochondrial calling succeeded for 1,037 of 1,102 WES samples (94.1%).
- **Results by subgroup.** Not available. The cohorts demonstrated workflow operation, not association, diagnosis or clinical prediction.

## 6. Known limitations and risks

- **Generalizability.** Validation is limited to the short-read data, reference assemblies, workflows, resources and datasets described in the publication and documentation. Performance is not established for long-read or somatic data or unsupported configurations.
- **Known biases and confounders.** Results may be affected by sequencing quality, coverage, capture design, batch effects, ancestry, genomic context, reference resources and limitations of the underlying tools. The two demonstration cohorts may retain technical or population differences.
- **Failure modes.** Poor or incomplete FASTQs, low coverage, incorrect metadata, incompatible resources, missing dependencies or unsupported parameters can cause failure or inaccurate calls. Outputs require quality control and, when appropriate, independent validation.
- **Privacy risks.** Model-specific attacks are not applicable. Inputs and outputs may contain sensitive genomic data, while logs may expose user-supplied identifiers or paths; all artifacts require appropriate data protection and access control.

## 7. Responsible use and reporting

- **Terms of use.** CBIcall is released under GPL-3.0-or-later for research use only. It is not validated for clinical decisions. Users are responsible for lawful data access, secure processing, compliance with source-data terms and appropriate validation.
- **Contact for misuse or vulnerability reports.** [manuel.rueda@cnag.eu](mailto:manuel.rueda@cnag.eu) or the [repository issue tracker](https://github.com/CNAG-Biomedical-Informatics/cbicall/issues). Do not disclose personal, patient or sensitive security information in public issues.

## 8. Acknowledgement

Funded by the European Union under grant agreement No 101137074 (HEREDITARY). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Health and Digital Executive Agency (HaDEA). Neither the European Union nor the granting authority can be held responsible for them.
