#!/bin/bash
#
# submit_cbicall.sh
# usage: submit_cbicall.sh <sample_id> <pipeline: wes|wgs>

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

# where your data and logs live
WORKDIR="/software/biomed/test/${PIPELINE}/${SAMPLE_ID}"

# name the generated job script
JOB_SCRIPT="job_${SAMPLE_ID}_${PIPELINE}.slurm"

cat > "${JOB_SCRIPT}" <<EOF
#!/bin/bash
#SBATCH --job-name=cbicall
#SBATCH -q normal
#SBATCH -D ${WORKDIR}
#SBATCH -e ${WORKDIR}/slurm-%N.%j.err
#SBATCH -o ${WORKDIR}/slurm-%N.%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -t 10:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=manuel.rueda@cnag.eu

module load Perl/5.36.0-GCCcore-12.2.0
eval "\$(perl -Mlocal::lib=/software/biomed/cbi_perl5)"

CBICALL_DIR="/software/biomed/cbicall"
CBICALL="\$CBICALL_DIR/bin/cbicall"
THREADS=4

cd \$SLURM_SUBMIT_DIR

# write a pipeline‐specific yaml
YAML_FILE="test.bash_${PIPELINE}_single.yaml"
cat <<YAML > "\${YAML_FILE}"
mode: single
pipeline: ${PIPELINE}
workflow_engine: bash
gatk_version: gatk-4.6
sample: ${WORKDIR}
YAML

srun "\$CBICALL" \\
     -p "\${YAML_FILE}" \\
     -t \$THREADS \\
     --no-color \\
     --no-emoji
EOF

# submit it
sbatch "${JOB_SCRIPT}"
