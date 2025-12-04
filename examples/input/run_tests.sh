#!/usr/bin/env bash
set -euo pipefail

############################################
# Config (can be overridden via env vars)  #
############################################

CBICALL=${CBICALL:-'../../bin/cbicall'}
THREADS=${THREADS:-1}
PATTERN=${PATTERN:-'#'}

RUN_WES=0
RUN_MIT=0

############################################
# Usage / argument parsing                 #
############################################

usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run cbicall unit tests.

You MUST specify at least one of:
  --wes      Run the WES test.
  --mit      Run the MIT test.

Options:
  --wes        Run the WES test.
  --mit        Run the MIT test.
  -h, --help   Show this help message and exit.

Environment variables:
  CBICALL   Path to cbicall executable (default: ../../bin/cbicall)
  THREADS   Number of threads to use (default: 1)
  PATTERN   Regex pattern to filter out lines before comparison (default: '#')
EOF
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --wes)
      RUN_WES=1
      shift
      ;;
    --mit)
      RUN_MIT=1
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
if [ "$RUN_WES" -eq 0 ] && [ "$RUN_MIT" -eq 0 ]; then
  echo "Error: you must specify at least one of --wes or --mit." >&2
  usage >&2
  exit 1
fi

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

  echo "Comparing:"
  echo "  REF : $ref"
  echo "  TEST: $tst"

  if diff -u \
      <("$filter_tool" -v -E "$PATTERN" "$ref" | LC_ALL=C sort) \
      <("$filter_tool" -v -E "$PATTERN" "$tst" | LC_ALL=C sort); then
    echo "SUCCESS: Files match."
    return 0
  else
    echo "ERROR: Differences found (see unified diff above)."
    return 1
  fi
}

############################################
# Overall status                           #
############################################

overall_status=0

############################################
# TEST 1: WES                              #
############################################

if [ "$RUN_WES" -eq 1 ]; then
  echo "========================================"
  echo "TEST: WES"
  echo "========================================"

  REF_VCF='CNAG999_exome/CNAG99901P_ex/ref_cbicall_bash_wes_single_gatk-4.6_648480891882296/02_varcall/CNAG99901P.hc.QC.vcf.gz'
  PARAM_WES='param.yaml'

  echo "Running WES unit test..."
  "$CBICALL" -p "$PARAM_WES" -t "$THREADS" > /dev/null 2>&1

  # Find latest WES result
  TEST_RESULT_WES=$(ls -t -- CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_gatk-4.6*/02_varcall/CNAG99901P.hc.QC.vcf.gz 2>/dev/null | head -n 1 || true)

  if [ -z "${TEST_RESULT_WES:-}" ]; then
    echo "ERROR: No WES result file found."
    overall_status=1
  else
    if ! compare_files "$REF_VCF" "$TEST_RESULT_WES" zgrep; then
      overall_status=1
    fi
  fi
else
  echo "Skipping WES test (not requested)."
fi

############################################
# TEST 2: MIT                              #
############################################

if [ "$RUN_MIT" -eq 1 ]; then
  echo "========================================"
  echo "TEST: MIT"
  echo "========================================"

  REF_MIT='CNAG999_exome/CNAG99901P_ex/ref_cbicall_bash_mit_single_gatk-3.5_648496901937882/01_mtoolbox/mit_prioritized_variants.txt'
  PARAM_MIT='mit_single.yaml'

  echo "Running MIT unit test..."
  "$CBICALL" -p "$PARAM_MIT" -t "$THREADS" > /dev/null 2>&1

  # Find latest MIT result
  TEST_RESULT_MIT=$(ls -t -- CNAG999_exome/CNAG99901P_ex/cbicall_bash_mit_single_gatk-3.5*/01_mtoolbox/mit_prioritized_variants.txt 2>/dev/null | head -n 1 || true)

  if [ -z "${TEST_RESULT_MIT:-}" ]; then
    echo "ERROR: No MIT result file found."
    overall_status=1
  else
    # MIT outputs are plain text, use grep
    if ! compare_files "$REF_MIT" "$TEST_RESULT_MIT" grep; then
      overall_status=1
    fi
  fi
else
  echo "Skipping MIT test (not requested)."
fi

echo "========================================"
echo "All requested tests finished."
echo "Exit code: $overall_status"
echo "========================================"

exit "$overall_status"

