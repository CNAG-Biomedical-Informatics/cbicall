#!/bin/bash
#
# run_cbicall_apptainer_slurm.sh
# usage: ./run_cbicall_apptainer_slurm.sh <sample_id> <pipeline: wes|wgs>

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

# Choose GenE SLURM settings based on the expected runtime. The standard
# research partition is limited to 12 hours; QoS is selected automatically.
if [ "$PIPELINE" = "wes" ]; then
  PARTITION="research"
  TIME="10:00:00"
elif [ "$PIPELINE" = "wgs" ]; then
  PARTITION="research_long"
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
SIF_IMAGE="/scratch_isilon/projects/0012-hereditary/software/containers/cbicall_1.2.0.sif"
CBICALL_DATA="/scratch_isilon/projects/0012-hereditary/software/cbicall-data"

cat > "${JOB_SCRIPT}" <<EOF
#!/bin/bash
#SBATCH --job-name=cbicall
#SBATCH --partition=${PARTITION}
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

apptainer exec \\
  --bind "${CBICALL_DATA}":/cbicall-data \\
  --bind "${WORKDIR}":"${WORKDIR}" \\
  --env CBICALL_DATA=/cbicall-data \\
  --pwd "${WORKDIR}" \\
  "${SIF_IMAGE}" \\
  cbicall \\
    run \\
    -p "\${YAML_FILE}" \\
    --runtime-profile cnag-hpc \\
    -t "\${SLURM_CPUS_PER_TASK}"
EOF

# submit it
sbatch "${JOB_SCRIPT}"
