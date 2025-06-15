#!/bin/bash
#
# run_cbicall_slurm.sh
# usage: ./run_cbicall_slurm.sh <sample_id> <pipeline: wes|wgs>

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

# use a simple ASCII locale
export LANG=C
export LC_ALL=C

module load Perl/5.36.0-GCCcore-12.2.0
eval "\$(perl -Mlocal::lib=/software/biomed/cbi_perl5)"

CBICALL_DIR="/software/biomed/cbicall"
CBICALL="\$CBICALL_DIR/bin/cbicall"

cd \$SLURM_SUBMIT_DIR

# write a pipeline‐specific yaml
YAML_FILE="${SAMPLE_ID}_${PIPELINE}_param.yaml"
cat <<YAML > "\${YAML_FILE}"
mode: single
pipeline: ${PIPELINE}
workflow_engine: bash
gatk_version: gatk-4.6
sample: ${WORKDIR}
projectdir: ${SAMPLE_ID}_cbicall
YAML

srun "\$CBICALL" \\
     -p "\$YAML_FILE" \\
     -t $THREADS \\
     --no-color \\
     --no-emoji
EOF

# submit it
sbatch "${JOB_SCRIPT}"
