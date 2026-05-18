#!/usr/bin/env bash
set -euo pipefail
LC_ALL=C

############################################
# Config (can be overridden via env vars)  #
############################################

CBICALL=${CBICALL:-'../../bin/cbicall'}
THREADS=${THREADS:-1}
PATTERN=${PATTERN:-'^#'}
LAUNCHER_LOG_DIR=${LAUNCHER_LOG_DIR:-$(pwd -P)}

RUN_WES_BASH=0
RUN_WES_SNAKEMAKE=0
RUN_WES_NEXTFLOW=0
RUN_MIT_BASH=0
SKIP_MISSING_OPTIONAL=${CBICALL_TEST_SKIP_MISSING_OPTIONAL:-0}

############################################
# Usage / argument parsing                 #
############################################

usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run cbicall integration tests.
Each requested test prints the new run directory, workflow log, run report,
launcher log, and outputs used for comparison.

You MUST specify at least one of:
  --wes-bash       Run the Bash WES test.
  --wes-snakemake  Run the Snakemake WES test.
  --wes-nextflow   Run the Nextflow WES test.
  --mit-bash       Run the Bash MIT test.

Options:
  --wes-bash       Run the Bash WES test.
  --wes-snakemake  Run the Snakemake WES test. Requires snakemake on PATH.
  --wes-nextflow   Run the Nextflow WES test. Requires nextflow on PATH.
  --mit-bash       Run the Bash MIT test.
  -h, --help       Show this help message and exit.

Environment variables:
  CBICALL                         Path to cbicall executable (default: ../../bin/cbicall)
  THREADS                         Number of threads to use (default: 1)
  PATTERN                         Regex pattern to filter out lines before comparison (default: '^#')
  CBICALL_TEST_SKIP_MISSING_OPTIONAL
                                  Skip optional engine tests when dependencies are missing.
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --wes-bash)
      RUN_WES_BASH=1
      shift
      ;;
    --wes-snakemake)
      RUN_WES_SNAKEMAKE=1
      shift
      ;;
    --wes-nextflow)
      RUN_WES_NEXTFLOW=1
      shift
      ;;
    --mit-bash)
      RUN_MIT_BASH=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

# Require at least one test to be selected
if [ "$RUN_WES_BASH" -eq 0 ] && [ "$RUN_WES_SNAKEMAKE" -eq 0 ] && [ "$RUN_WES_NEXTFLOW" -eq 0 ] && [ "$RUN_MIT_BASH" -eq 0 ]; then
  echo "Error: you must specify at least one of --wes-bash, --wes-snakemake, --wes-nextflow, or --mit-bash." >&2
  usage >&2
  exit 1
fi

############################################
# Helper: check reference directory        #
############################################
# Args:
#   1: path to reference file
#
# This checks that the directory containing the reference
# file exists before we try to run or compare anything.

check_ref_dir() {
  local ref_path="$1"
  local ref_dir

  ref_dir=$(dirname "$ref_path")

  if [ ! -d "$ref_dir" ]; then
    echo "ERROR: Reference directory does not exist: $ref_dir" >&2
    echo "       (from reference path: $ref_path)" >&2
    exit 1
  fi
}

############################################
# Helper: compare files                    #
############################################
# Args:
#   1: reference file
#   2: test file
#   3: filter tool (zgrep or grep)

compare_files() {
  local ref="$1"
  local tst="$2"
  local filter_tool="$3"
  local ref_filtered
  local tst_filtered

  echo "Comparing:"
  echo "  REF : $ref"
  echo "  TEST: $tst"

  if [ ! -f "$tst" ]; then
    echo "ERROR: Test output does not exist: $tst"
    return 1
  fi

  ref_filtered=$(mktemp)
  tst_filtered=$(mktemp)
  "$filter_tool" -v -E "$PATTERN" "$ref" | sort > "$ref_filtered"
  "$filter_tool" -v -E "$PATTERN" "$tst" | sort > "$tst_filtered"

  echo "Normalized SHA-256 (pattern: $PATTERN):"
  printf "  REF : "
  sha256sum "$ref_filtered" | awk '{print $1}'
  printf "  TEST: "
  sha256sum "$tst_filtered" | awk '{print $1}'

  if diff -u "$ref_filtered" "$tst_filtered"; then
    rm -f "$ref_filtered" "$tst_filtered"
    echo "SUCCESS: Files match."
    return 0
  else
    rm -f "$ref_filtered" "$tst_filtered"
    echo "ERROR: Differences found (see unified diff above)."
    return 1
  fi
}

############################################
# Helper: run directory tracking           #
############################################

list_run_dirs() {
  local base_dir="$1"
  local run_glob="$2"

  find "$base_dir" -maxdepth 1 -type d -name "$run_glob" -print 2>/dev/null | sort
}

new_run_dir() {
  local before_file="$1"
  local after_file="$2"

  comm -13 "$before_file" "$after_file" | tail -n 1
}

run_dir_from_launcher_log() {
  local launcher_log="$1"
  local path

  path=$(sed -n 's/^Working directory: //p' "$launcher_log" | tail -n 1)
  if [ -n "$path" ]; then
    echo "$path"
    return 0
  fi

  path=$(sed -n 's/^  Report[[:space:]]*=> //p' "$launcher_log" | tail -n 1)
  if [ -n "$path" ]; then
    dirname "$path"
    return 0
  fi

  path=$(sed -n 's/^  Log[[:space:]]*=> //p' "$launcher_log" | tail -n 1)
  if [ -n "$path" ]; then
    dirname "$path"
    return 0
  fi
}

print_run_artifacts() {
  local label="$1"
  local run_dir="$2"
  local workflow_log="$3"
  local report="$4"
  local launcher_log="$5"
  shift 5

  echo
  echo "$label run artifacts:"
  print_artifact "Run dir" "$run_dir"
  print_artifact "Workflow log" "$workflow_log"
  print_artifact "Run report" "$report"
  print_artifact "Launcher log" "$launcher_log"

  for output in "$@"; do
    print_artifact "Output" "$output"
  done
  echo
}

print_artifact() {
  local label="$1"
  local path="$2"

  if [ -e "$path" ]; then
    printf "  %-12s : %s\n" "$label" "$path"
  else
    printf "  %-12s : %s (not created)\n" "$label" "$path"
  fi
}

append_summary() {
  local label="$1"
  local status="$2"
  local detail="${3:-}"

  if [ -n "$detail" ]; then
    summary="${summary}${label}: ${status} (${detail})"$'\n'
  else
    summary="${summary}${label}: ${status}"$'\n'
  fi
}

############################################
# Overall status                           #
############################################

overall_status=0
summary=""

############################################
# TEST 1: WES                              #
############################################

if [ "$RUN_WES_BASH" -eq 1 ]; then
  echo "========================================"
  echo "TEST: WES Bash"
  echo "========================================"

  REF_VCF='CNAG999_exome/CNAG99901P_ex/ref_cbicall_bash_wes_single_b37_gatk-4.6_765963065360466/02_varcall/CNAG99901P.hc.QC.vcf.gz'
  PARAM_WES='param.yaml'
  BASE_DIR='CNAG999_exome/CNAG99901P_ex'
  RUN_GLOB='cbicall_bash_wes_single_b37_gatk-4.6_*'
  WORKFLOW_LOG_NAME='bash_wes_single_b37_gatk-4.6.log'

  # Ensure reference directory exists
  check_ref_dir "$REF_VCF"

  echo "Running WES integration test..."
  WES_STATUS=0
  WES_BEFORE=$(mktemp)
  WES_AFTER=$(mktemp)
  WES_LAUNCHER_LOG=$(mktemp -p "$LAUNCHER_LOG_DIR" cbicall-test-wes.XXXXXX.log)
  list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_BEFORE"

  if ! "$CBICALL" run -p "$PARAM_WES" -t "$THREADS" > "$WES_LAUNCHER_LOG" 2>&1; then
    list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_AFTER"
    WES_RUN_DIR=$(run_dir_from_launcher_log "$WES_LAUNCHER_LOG")
    if [ -z "${WES_RUN_DIR:-}" ]; then
      WES_RUN_DIR=$(new_run_dir "$WES_BEFORE" "$WES_AFTER" || true)
    fi
    echo "ERROR: WES Bash cbicall command failed. Last launcher log lines:"
    tail -n 40 "$WES_LAUNCHER_LOG" || true
    if [ -n "${WES_RUN_DIR:-}" ]; then
      print_run_artifacts "WES Bash" "$WES_RUN_DIR" "$WES_RUN_DIR/$WORKFLOW_LOG_NAME" "$WES_RUN_DIR/run-report.json" "$WES_LAUNCHER_LOG"
    else
      echo "Launcher log: $WES_LAUNCHER_LOG"
    fi
    WES_STATUS=1
    append_summary "WES Bash" "failed" "$WES_LAUNCHER_LOG"
  else
    list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_AFTER"
    WES_RUN_DIR=$(run_dir_from_launcher_log "$WES_LAUNCHER_LOG")
    if [ -z "${WES_RUN_DIR:-}" ]; then
      WES_RUN_DIR=$(new_run_dir "$WES_BEFORE" "$WES_AFTER" || true)
    fi

    if [ -z "${WES_RUN_DIR:-}" ]; then
      echo "ERROR: WES Bash command finished but no new run directory was found."
      echo "Launcher log: $WES_LAUNCHER_LOG"
      WES_STATUS=1
      append_summary "WES Bash" "failed" "no new run directory"
    else
      TEST_RESULT_WES="$WES_RUN_DIR/02_varcall/CNAG99901P.hc.QC.vcf.gz"
      print_run_artifacts "WES Bash" "$WES_RUN_DIR" "$WES_RUN_DIR/$WORKFLOW_LOG_NAME" "$WES_RUN_DIR/run-report.json" "$WES_LAUNCHER_LOG" "$TEST_RESULT_WES"

      if compare_files "$REF_VCF" "$TEST_RESULT_WES" zgrep; then
        append_summary "WES Bash" "passed" "$WES_RUN_DIR"
      else
        WES_STATUS=1
        append_summary "WES Bash" "failed" "$WES_RUN_DIR"
      fi
    fi
  fi

  if [ "$WES_STATUS" -ne 0 ]; then
    overall_status=1
  fi
  rm -f "$WES_BEFORE" "$WES_AFTER"
else
  echo "Skipping WES Bash test (not requested)."
  append_summary "WES Bash" "skipped"
fi

############################################
# TEST 1B: WES Snakemake                   #
############################################

if [ "$RUN_WES_SNAKEMAKE" -eq 1 ]; then
  echo "========================================"
  echo "TEST: WES Snakemake"
  echo "========================================"

  REF_VCF='CNAG999_exome/CNAG99901P_ex/ref_cbicall_bash_wes_single_b37_gatk-4.6_765963065360466/02_varcall/CNAG99901P.hc.QC.vcf.gz'
  BASE_DIR='CNAG999_exome/CNAG99901P_ex'
  RUN_GLOB='cbicall_snakemake_wes_single_b37_gatk-4.6_*'
  WORKFLOW_LOG_NAME='snakemake_wes_single_b37_gatk-4.6.log'

  check_ref_dir "$REF_VCF"

  if ! command -v snakemake >/dev/null 2>&1; then
    if [ "$SKIP_MISSING_OPTIONAL" -eq 1 ]; then
      echo "SKIP: Snakemake is not installed; install it to run Snakemake backend integration tests."
      append_summary "WES Snakemake" "skipped" "snakemake not found"
    else
      echo "ERROR: Snakemake is not installed; install it or omit --wes-snakemake." >&2
      append_summary "WES Snakemake" "failed" "snakemake not found"
      overall_status=1
    fi
  else
    echo "Running WES Snakemake integration test..."
    WES_SMK_STATUS=0
    WES_SMK_BEFORE=$(mktemp)
    WES_SMK_AFTER=$(mktemp)
    WES_SMK_PARAM=$(mktemp -p "$LAUNCHER_LOG_DIR" cbicall-wes-snakemake.XXXXXX.yaml)
    WES_SMK_LAUNCHER_LOG=$(mktemp -p "$LAUNCHER_LOG_DIR" cbicall-test-wes-snakemake.XXXXXX.log)
    list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_SMK_BEFORE"

    cat > "$WES_SMK_PARAM" <<EOF
mode: single
pipeline: wes
workflow_engine: snakemake
profile: local
gatk_version: gatk-4.6
resource: "cbicall-germline-resources-v1"
input_dir: CNAG999_exome/CNAG99901P_ex
EOF

    if ! "$CBICALL" run -p "$WES_SMK_PARAM" -t "$THREADS" > "$WES_SMK_LAUNCHER_LOG" 2>&1; then
      list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_SMK_AFTER"
      WES_SMK_RUN_DIR=$(run_dir_from_launcher_log "$WES_SMK_LAUNCHER_LOG")
      if [ -z "${WES_SMK_RUN_DIR:-}" ]; then
        WES_SMK_RUN_DIR=$(new_run_dir "$WES_SMK_BEFORE" "$WES_SMK_AFTER" || true)
      fi
      echo "ERROR: WES Snakemake cbicall command failed. Last launcher log lines:"
      tail -n 40 "$WES_SMK_LAUNCHER_LOG" || true
      if [ -n "${WES_SMK_RUN_DIR:-}" ]; then
        print_run_artifacts "WES Snakemake" "$WES_SMK_RUN_DIR" "$WES_SMK_RUN_DIR/$WORKFLOW_LOG_NAME" "$WES_SMK_RUN_DIR/run-report.json" "$WES_SMK_LAUNCHER_LOG"
      else
        echo "Launcher log: $WES_SMK_LAUNCHER_LOG"
      fi
      WES_SMK_STATUS=1
      append_summary "WES Snakemake" "failed" "$WES_SMK_LAUNCHER_LOG"
    else
      list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_SMK_AFTER"
      WES_SMK_RUN_DIR=$(run_dir_from_launcher_log "$WES_SMK_LAUNCHER_LOG")
      if [ -z "${WES_SMK_RUN_DIR:-}" ]; then
        WES_SMK_RUN_DIR=$(new_run_dir "$WES_SMK_BEFORE" "$WES_SMK_AFTER" || true)
      fi

      if [ -z "${WES_SMK_RUN_DIR:-}" ]; then
        echo "ERROR: WES Snakemake command finished but no new run directory was found."
        echo "Launcher log: $WES_SMK_LAUNCHER_LOG"
        WES_SMK_STATUS=1
        append_summary "WES Snakemake" "failed" "no new run directory"
      else
        TEST_RESULT_WES_SMK="$WES_SMK_RUN_DIR/02_varcall/CNAG99901P.hc.QC.vcf.gz"
        print_run_artifacts "WES Snakemake" "$WES_SMK_RUN_DIR" "$WES_SMK_RUN_DIR/$WORKFLOW_LOG_NAME" "$WES_SMK_RUN_DIR/run-report.json" "$WES_SMK_LAUNCHER_LOG" "$TEST_RESULT_WES_SMK"

        if compare_files "$REF_VCF" "$TEST_RESULT_WES_SMK" zgrep; then
          append_summary "WES Snakemake" "passed" "$WES_SMK_RUN_DIR"
        else
          WES_SMK_STATUS=1
          append_summary "WES Snakemake" "failed" "$WES_SMK_RUN_DIR"
        fi
      fi
    fi

    if [ "$WES_SMK_STATUS" -ne 0 ]; then
      overall_status=1
    fi
    rm -f "$WES_SMK_BEFORE" "$WES_SMK_AFTER" "$WES_SMK_PARAM"
  fi
else
  echo "Skipping WES Snakemake test (not requested)."
  append_summary "WES Snakemake" "skipped"
fi

############################################
# TEST 1C: WES Nextflow                    #
############################################

if [ "$RUN_WES_NEXTFLOW" -eq 1 ]; then
  echo "========================================"
  echo "TEST: WES Nextflow"
  echo "========================================"

  REF_VCF='CNAG999_exome/CNAG99901P_ex/ref_cbicall_bash_wes_single_b37_gatk-4.6_765963065360466/02_varcall/CNAG99901P.hc.QC.vcf.gz'
  BASE_DIR='CNAG999_exome/CNAG99901P_ex'
  RUN_GLOB='cbicall_nextflow_wes_single_b37_gatk-4.6_*'
  WORKFLOW_LOG_NAME='nextflow_wes_single_b37_gatk-4.6.log'

  check_ref_dir "$REF_VCF"

  if ! command -v nextflow >/dev/null 2>&1; then
    if [ "$SKIP_MISSING_OPTIONAL" -eq 1 ]; then
      echo "SKIP: Nextflow is not installed; install it to run Nextflow backend integration tests."
      append_summary "WES Nextflow" "skipped" "nextflow not found"
    else
      echo "ERROR: Nextflow is not installed; install it or omit --wes-nextflow." >&2
      append_summary "WES Nextflow" "failed" "nextflow not found"
      overall_status=1
    fi
  else
    echo "Running WES Nextflow integration test..."
    WES_NF_STATUS=0
    WES_NF_BEFORE=$(mktemp)
    WES_NF_AFTER=$(mktemp)
    WES_NF_PARAM=$(mktemp -p "$LAUNCHER_LOG_DIR" cbicall-wes-nextflow.XXXXXX.yaml)
    WES_NF_LAUNCHER_LOG=$(mktemp -p "$LAUNCHER_LOG_DIR" cbicall-test-wes-nextflow.XXXXXX.log)
    list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_NF_BEFORE"

    cat > "$WES_NF_PARAM" <<EOF
mode: single
pipeline: wes
workflow_engine: nextflow
profile: local
gatk_version: gatk-4.6
resource: "cbicall-germline-resources-v1"
input_dir: CNAG999_exome/CNAG99901P_ex
EOF

    if ! "$CBICALL" run -p "$WES_NF_PARAM" -t "$THREADS" > "$WES_NF_LAUNCHER_LOG" 2>&1; then
      list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_NF_AFTER"
      WES_NF_RUN_DIR=$(run_dir_from_launcher_log "$WES_NF_LAUNCHER_LOG")
      if [ -z "${WES_NF_RUN_DIR:-}" ]; then
        WES_NF_RUN_DIR=$(new_run_dir "$WES_NF_BEFORE" "$WES_NF_AFTER" || true)
      fi
      echo "ERROR: WES Nextflow cbicall command failed. Last launcher log lines:"
      tail -n 40 "$WES_NF_LAUNCHER_LOG" || true
      if [ -n "${WES_NF_RUN_DIR:-}" ]; then
        print_run_artifacts "WES Nextflow" "$WES_NF_RUN_DIR" "$WES_NF_RUN_DIR/$WORKFLOW_LOG_NAME" "$WES_NF_RUN_DIR/run-report.json" "$WES_NF_LAUNCHER_LOG"
      else
        echo "Launcher log: $WES_NF_LAUNCHER_LOG"
      fi
      WES_NF_STATUS=1
      append_summary "WES Nextflow" "failed" "$WES_NF_LAUNCHER_LOG"
    else
      list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$WES_NF_AFTER"
      WES_NF_RUN_DIR=$(run_dir_from_launcher_log "$WES_NF_LAUNCHER_LOG")
      if [ -z "${WES_NF_RUN_DIR:-}" ]; then
        WES_NF_RUN_DIR=$(new_run_dir "$WES_NF_BEFORE" "$WES_NF_AFTER" || true)
      fi

      if [ -z "${WES_NF_RUN_DIR:-}" ]; then
        echo "ERROR: WES Nextflow command finished but no new run directory was found."
        echo "Launcher log: $WES_NF_LAUNCHER_LOG"
        WES_NF_STATUS=1
        append_summary "WES Nextflow" "failed" "no new run directory"
      else
        TEST_RESULT_WES_NF="$WES_NF_RUN_DIR/02_varcall/CNAG99901P.hc.QC.vcf.gz"
        print_run_artifacts "WES Nextflow" "$WES_NF_RUN_DIR" "$WES_NF_RUN_DIR/$WORKFLOW_LOG_NAME" "$WES_NF_RUN_DIR/run-report.json" "$WES_NF_LAUNCHER_LOG" "$TEST_RESULT_WES_NF"

        if compare_files "$REF_VCF" "$TEST_RESULT_WES_NF" zgrep; then
          append_summary "WES Nextflow" "passed" "$WES_NF_RUN_DIR"
        else
          WES_NF_STATUS=1
          append_summary "WES Nextflow" "failed" "$WES_NF_RUN_DIR"
        fi
      fi
    fi

    if [ "$WES_NF_STATUS" -ne 0 ]; then
      overall_status=1
    fi
    rm -f "$WES_NF_BEFORE" "$WES_NF_AFTER" "$WES_NF_PARAM"
  fi
else
  echo "Skipping WES Nextflow test (not requested)."
  append_summary "WES Nextflow" "skipped"
fi

############################################
# TEST 2: MIT                              #
############################################

if [ "$RUN_MIT_BASH" -eq 1 ]; then
  echo "========================================"
  echo "TEST: MIT Bash"
  echo "========================================"

  REF_MIT='CNAG999_exome/CNAG99901P_ex/ref_cbicall_bash_mit_single_rsrs_gatk-3.5_649547582283533/01_mtoolbox/mit_prioritized_variants.txt'
  REF_MIT_JSON='CNAG999_exome/CNAG99901P_ex/ref_cbicall_bash_mit_single_rsrs_gatk-3.5_649547582283533/01_mtoolbox/mit.raw.json'
  PARAM_MIT='mit_single.yaml'
  BASE_DIR='CNAG999_exome/CNAG99901P_ex'
  RUN_GLOB='cbicall_bash_mit_single_rsrs_gatk-3.5_*'
  WORKFLOW_LOG_NAME='bash_mit_single_rsrs_gatk-3.5.log'

  # Ensure reference files exist
  check_ref_dir "$REF_MIT"
  check_ref_dir "$REF_MIT_JSON"

  echo "Running MIT integration test..."
  MIT_STATUS=0
  MIT_BEFORE=$(mktemp)
  MIT_AFTER=$(mktemp)
  MIT_LAUNCHER_LOG=$(mktemp -p "$LAUNCHER_LOG_DIR" cbicall-test-mit.XXXXXX.log)
  list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$MIT_BEFORE"

  if ! "$CBICALL" run -p "$PARAM_MIT" -t "$THREADS" > "$MIT_LAUNCHER_LOG" 2>&1; then
    list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$MIT_AFTER"
    MIT_RUN_DIR=$(run_dir_from_launcher_log "$MIT_LAUNCHER_LOG")
    if [ -z "${MIT_RUN_DIR:-}" ]; then
      MIT_RUN_DIR=$(new_run_dir "$MIT_BEFORE" "$MIT_AFTER" || true)
    fi
    echo "ERROR: MIT Bash cbicall command failed. Last launcher log lines:"
    tail -n 40 "$MIT_LAUNCHER_LOG" || true
    if [ -n "${MIT_RUN_DIR:-}" ]; then
      print_run_artifacts "MIT Bash" "$MIT_RUN_DIR" "$MIT_RUN_DIR/$WORKFLOW_LOG_NAME" "$MIT_RUN_DIR/run-report.json" "$MIT_LAUNCHER_LOG"
    else
      echo "Launcher log: $MIT_LAUNCHER_LOG"
    fi
    MIT_STATUS=1
    append_summary "MIT Bash" "failed" "$MIT_LAUNCHER_LOG"
  else
    list_run_dirs "$BASE_DIR" "$RUN_GLOB" > "$MIT_AFTER"
    MIT_RUN_DIR=$(run_dir_from_launcher_log "$MIT_LAUNCHER_LOG")
    if [ -z "${MIT_RUN_DIR:-}" ]; then
      MIT_RUN_DIR=$(new_run_dir "$MIT_BEFORE" "$MIT_AFTER" || true)
    fi

    if [ -z "${MIT_RUN_DIR:-}" ]; then
      echo "ERROR: MIT Bash command finished but no new run directory was found."
      echo "Launcher log: $MIT_LAUNCHER_LOG"
      MIT_STATUS=1
      append_summary "MIT Bash" "failed" "no new run directory"
    else
      TEST_RESULT_MIT="$MIT_RUN_DIR/01_mtoolbox/mit_prioritized_variants.txt"
      TEST_RESULT_MIT_JSON="$MIT_RUN_DIR/01_mtoolbox/mit.raw.json"
      print_run_artifacts "MIT Bash" "$MIT_RUN_DIR" "$MIT_RUN_DIR/$WORKFLOW_LOG_NAME" "$MIT_RUN_DIR/run-report.json" "$MIT_LAUNCHER_LOG" "$TEST_RESULT_MIT" "$TEST_RESULT_MIT_JSON"

      if ! compare_files "$REF_MIT" "$TEST_RESULT_MIT" grep; then
        MIT_STATUS=1
      fi

      if [ ! -f "$TEST_RESULT_MIT_JSON" ]; then
        echo "ERROR: No MIT raw JSON result file found: $TEST_RESULT_MIT_JSON"
        MIT_STATUS=1
      elif ! diff -u "$REF_MIT_JSON" "$TEST_RESULT_MIT_JSON" >/dev/null; then
        echo "DIFF: mit.raw.json differs from reference:"
        diff -u "$REF_MIT_JSON" "$TEST_RESULT_MIT_JSON" || true
        MIT_STATUS=1
      else
        echo "SUCCESS: mit.raw.json matches."
      fi

      if [ "$MIT_STATUS" -eq 0 ]; then
        append_summary "MIT Bash" "passed" "$MIT_RUN_DIR"
      else
        append_summary "MIT Bash" "failed" "$MIT_RUN_DIR"
      fi
    fi
  fi

  if [ "$MIT_STATUS" -ne 0 ]; then
    overall_status=1
  fi
  rm -f "$MIT_BEFORE" "$MIT_AFTER"

else
  echo "Skipping MIT Bash test (not requested)."
  append_summary "MIT Bash" "skipped"
fi

echo "========================================"
echo "Integration test summary"
echo "========================================"
printf "%s" "$summary"
echo "========================================"
echo "All requested tests finished."
echo "Exit code: $overall_status"
echo "========================================"

exit "$overall_status"
