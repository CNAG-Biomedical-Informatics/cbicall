#!/bin/bash
set -eu 

dir=/media/mrueda/2TBS/CNAG/Project_CBI_Call
cbicall=/media/mrueda/2TBS/CNAG/Project_CBI_Call/cbicall/bin/cbicall
ncpu=4

for dirname in MA99999_exome
do
 cd $dir/$dirname
 echo $dirname
 for sample in MA*ex
 do
  echo "...$sample"
  cd $sample
  cat<<EOF>$sample.mit_single.yaml
mode: single
pipeline: mit
sample: $dir/$dirname/$sample
EOF
$cbicall -t $ncpu -p $sample.mit_single.yaml > $sample.mit_single.log 2>&1
  cd ..
 done
cd ..
done
