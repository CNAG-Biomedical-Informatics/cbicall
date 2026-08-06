#!/usr/bin/env bash
set -euo pipefail

# Submit one containerized WES/WGS run through Slurm.
# Usage: ./run_cbicall_apptainer_slurm.sh <sample_id> <wes|wgs> <FASTQ directory>

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <sample_id> <wes|wgs> <FASTQ directory>" >&2
  exit 2
fi

SAMPLE_ID=$1
PIPELINE=$2
INPUT_DIR=$(cd "$3" && pwd -P)

if [ "$PIPELINE" != "wes" ] && [ "$PIPELINE" != "wgs" ]; then
  echo "ERROR: pipeline must be 'wes' or 'wgs'" >&2
  exit 2
fi

# Edit these two paths.
SIF_IMAGE="/path/to/cbicall.sif"
CBICALL_DATA="/path/to/cbicall-data"

if [ ! -f "$SIF_IMAGE" ]; then
  echo "ERROR: CBIcall image not found: $SIF_IMAGE" >&2
  exit 2
fi
if [ ! -d "$CBICALL_DATA" ]; then
  echo "ERROR: resource directory not found: $CBICALL_DATA" >&2
  exit 2
fi

THREADS=4
MEM=24G
if [ "$PIPELINE" = "wes" ]; then
  PARTITION=research
  TIME=10:00:00
  GENOME=b37
else
  PARTITION=research_long
  TIME=2-00:00:00
  GENOME=hg38
fi

PARAM_FILE="$INPUT_DIR/${SAMPLE_ID}_${PIPELINE}_single_apptainer.yaml"
JOB_SCRIPT="$INPUT_DIR/job_${SAMPLE_ID}_${PIPELINE}_apptainer.slurm"

cat > "$PARAM_FILE" <<EOF
mode: single
pipeline: $PIPELINE
workflow_provider: cbicall-core
workflow_backend: bash
software_stack: gatk-4.6
genome: $GENOME
resource: cbicall-germline-resources-v1
input_dir: "$INPUT_DIR"
cleanup_bam: false
EOF

cat > "$JOB_SCRIPT" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=cbicall-${SAMPLE_ID}
#SBATCH --partition=$PARTITION
#SBATCH --chdir=$INPUT_DIR
#SBATCH --error=$INPUT_DIR/slurm-%N.%j.err
#SBATCH --output=$INPUT_DIR/slurm-%N.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$MEM
#SBATCH --time=$TIME

set -euo pipefail

module load apptainer 2>/dev/null || true

apptainer exec \
  --bind "$CBICALL_DATA":/cbicall-data \
  --bind "$INPUT_DIR":"$INPUT_DIR" \
  --env CBICALL_DATA=/cbicall-data \
  --pwd "$INPUT_DIR" \
  "$SIF_IMAGE" \
  cbicall run \
    -p "$PARAM_FILE" \
    --runtime-profile cnag-hpc \
    -t "\${SLURM_CPUS_PER_TASK}"
EOF

echo "Parameters: $PARAM_FILE"
echo "Job script: $JOB_SCRIPT"
sbatch "$JOB_SCRIPT"
