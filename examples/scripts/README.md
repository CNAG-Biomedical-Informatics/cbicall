# CBIcall script examples

These files are readable launch recipes, not a second configuration layer.
Each script shows the complete YAML submitted to CBIcall followed by the
command that runs it.

## Local examples

| Script | Analysis |
| --- | --- |
| `run_wes_single.sh` | Single-sample WES |
| `run_wgs_single.sh` | Single-sample WGS |
| `run_wes_cohort.sh` | WES joint genotyping from a sample map |
| `run_wgs_cohort.sh` | WGS joint genotyping from a sample map |
| `run_mit_single.sh` | Single-sample mtDNA from a completed WES/WGS run |
| `run_mit_cohort.sh` | Cohort mtDNA from completed WES/WGS runs |

Open the required script and change the path marked near the top. Then run it
without arguments. For example:

```bash
vim examples/scripts/run_wes_single.sh
./examples/scripts/run_wes_single.sh
```

The WES example writes `wes_single.yaml` and launches the visible command:

```bash
cbicall run -p wes_single.yaml -t 4
```

When using a source checkout, either replace `cbicall` with
`/path/to/cbicall/bin/cbicall` or add that directory to `PATH`.

The cohort scripts expect a two-column, tab-separated sample map:

```text
sample01  /data/sample01/02_varcall/sample01.hc.g.vcf.gz
sample02  /data/sample02/02_varcall/sample02.hc.g.vcf.gz
```

## Slurm examples

The Slurm launchers accept a sample ID, pipeline, and FASTQ directory because
those values normally change for every submitted job. They write both the YAML
and generated Slurm job file before calling `sbatch`.

The source-checkout example is configured for the CNAG GenE cluster:

```bash
./examples/scripts/run_cbicall_slurm.sh sample01 wes /data/sample01
```

It shows the selected partition, requested CPUs and memory, Python module,
`PYTHONPATH`, runtime profile, and final `cbicall run` command in the generated
job file. WES uses `research`; WGS uses `research_long`.

For Apptainer, first edit `SIF_IMAGE` and `CBICALL_DATA` near the top of the
script, then submit in the same way:

```bash
./examples/scripts/run_cbicall_apptainer_slurm.sh sample01 wes /data/sample01
```

Both launchers print the YAML and job-script paths before submission so they
can be inspected directly.
