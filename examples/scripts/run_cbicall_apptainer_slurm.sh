#!/bin/bash
#
# run_cbicall_slurm_apptainer.sh
# usage: ./run_cbicall_slurm_apptainer.sh <sample_id> <pipeline: wes|wgs>

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <sample_id> <pipeline: wes|wgs>"
  exit 1
fi

SAMPLE_ID=$1
PIPELINE=$2

if [[ "$PIPELINE" != "wes" && "$PIPELINE" != "wgs" ]]; then
  echo "Error: pipeline must be 'wes' or 'wgs'"
  exit 1
fi

# choose SLURM settings based on pipeline
if [ "$PIPELINE" = "wes" ]; then
  QUEUE="normal"
  TIME="10:00:00"
elif [ "$PIPELINE" = "wgs" ]; then
  QUEUE="vlong"
  TIME="2-00:00:00"
fi

# Uppercase version of pipeline
PIPELINE_UC=${PIPELINE^^}

# where your data and logs live
WORKDIR="/scratch_isilon/projects/0012-hereditary/dbgap/fastq/phs001585/${PIPELINE_UC}/${SAMPLE_ID}"

# name the generated job script
JOB_SCRIPT="job_${SAMPLE_ID}_${PIPELINE}.slurm"

# Number of threads
THREADS=4

# RAM (x1.5 to help prevent oom-kills)
MEM="24G"

# Apptainer settings (edit as needed)
SIF_IMAGE="/software/biomed/containers/cbicall_1.2.0.sif"
CBICALL_DATA="/software/biomed/cbicall-data"

cat > "${JOB_SCRIPT}" <<EOF
#!/bin/bash
#SBATCH --job-name=cbicall
#SBATCH -q ${QUEUE}
#SBATCH -D ${WORKDIR}
#SBATCH -e ${WORKDIR}/slurm-%N.%j.err
#SBATCH -o ${WORKDIR}/slurm-%N.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=${MEM}
#SBATCH -t ${TIME}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=manuel.rueda@cnag.eu

set -euo pipefail

module load apptainer 2>/dev/null || true

# Sanity checks
if [ ! -f "${SIF_IMAGE}" ]; then
  echo "ERROR: SIF image not found: ${SIF_IMAGE}"
  exit 2
fi

if [ ! -d "${CBICALL_DATA}" ]; then
  echo "ERROR: CBICALL_DATA directory not found: ${CBICALL_DATA}"
  exit 2
fi

cd "${WORKDIR}"

# write a pipeline-specific yaml
YAML_FILE="${SAMPLE_ID}_${PIPELINE}_param.yaml"
cat <<YAML > "\${YAML_FILE}"
mode: single
pipeline: ${PIPELINE}
workflow_backend: bash
software_stack: gatk-4.6
resource: "cbicall-germline-resources-v1"
input_dir: ${WORKDIR}
project_dir: ${SAMPLE_ID}_cbicall
cleanup_bam: false
YAML

# Run cbicall inside the container
# - Bind databases to /cbicall-data
# - Bind WORKDIR so paths referenced in the YAML exist inside the container

srun apptainer exec \\
  --bind "${CBICALL_DATA}":/cbicall-data \\
  --bind "${WORKDIR}":"${WORKDIR}" \\
  --env CBICALL_DATA=/cbicall-data \\
  --pwd "${WORKDIR}" \\
  "${SIF_IMAGE}" \\
  cbicall \\
    run \\
    -p "\${YAML_FILE}" \\
    -t ${THREADS}
EOF

# submit it
sbatch "${JOB_SCRIPT}"
