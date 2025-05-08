#!/bin/bash
set -eu 

dir=/media/mrueda/2TBS/CNAG/Project_CBI_Call
cbicall=/media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/bin/cbicall
ncpu=4

for dirname in MA99999_exome
do
 cd $dir/$dirname
 echo $dirname
  echo "...$dirname"
  cat<<EOF> $dirname.wes_cohort.yaml
mode: cohort
pipeline: wes
sample: $dir/$dirname
EOF

$cbicall -t $ncpu -p $dirname.wes_cohort.yaml 
cd ..
done
