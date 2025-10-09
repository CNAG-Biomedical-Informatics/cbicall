#!/usr/bin/env bash
set -euo pipefail

# Fixed reference file (left as-is per your request)
REF='CNAG999_exome/CNAG99901P_ex/ref_cbicall_bash_wes_single_gatk-4.6_175801291279888/02_varcall/CNAG99901P.hc.QC.vcf.gz'
PAT='#'
CBICALL='../../bin/cbicall'
PARAM='param.yaml'
THREADS=1

# Run Unit test
echo "Running unit test..."
$CBICALL -p $PARAM -t $THREADS > /dev/null 2>&1

# Identify the latest result
TEST_RESULT=$(ls -t -- CNAG999_exome/CNAG99901P_ex/cbicall_bash_wes_single_gatk-4.6*/02_varcall/CNAG99901P.hc.QC.vcf.gz | head -n 1)

echo "Comparing:"
echo "$REF"
echo "$TEST_RESULT"

# Order-independent compare: exclude PAT, sort, then unified diff
if diff -u <(zgrep -v -E "$PAT" $REF | LC_ALL=C sort) <(zgrep -v -E "$PAT" $TEST_RESULT | LC_ALL=C sort); then
  echo "SUCCESS: Test ran fine."
else
  echo "ERROR: Differences found (see unified diff above)."
fi
