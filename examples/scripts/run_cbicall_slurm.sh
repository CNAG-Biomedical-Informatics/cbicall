#!/usr/bin/env bash
set -euo pipefail

# Submit one native Bash WES/WGS run to the CNAG GenE cluster.
# Usage: ./run_cbicall_slurm.sh <sample_id> <wes|wgs> <FASTQ directory>

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

# CNAG GenE installation. Edit these paths for another installation.
CBICALL_DIR="/scratch_isilon/projects/0012-hereditary/software/cbicall"
CBICALL_PYTHON_PREFIX="/scratch_isilon/projects/0012-hereditary/software/cbi_py3"

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

PARAM_FILE="$INPUT_DIR/${SAMPLE_ID}_${PIPELINE}_single.yaml"
JOB_SCRIPT="$INPUT_DIR/job_${SAMPLE_ID}_${PIPELINE}.slurm"

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

module load Python/3.13.5-GCCcore-14.3.0
export PYTHONPATH="$CBICALL_PYTHON_PREFIX/lib/python3.13/site-packages\${PYTHONPATH:+:\$PYTHONPATH}"

"$CBICALL_DIR/bin/cbicall" run \
  -p "$PARAM_FILE" \
  --runtime-profile cnag-hpc \
  -t "\${SLURM_CPUS_PER_TASK}"
EOF

echo "Parameters: $PARAM_FILE"
echo "Job script: $JOB_SCRIPT"
sbatch "$JOB_SCRIPT"
