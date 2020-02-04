#!/bin/bash
# this script run prep_test.in script file (that had been prepared befor) to run analysis in loop with cpptraj for any number of md trajectories
# see the tutorial https://ambermd.org/tutorials/TrajectoryAnalysis.php
output=$(pwd)

for sim in "${output}"/* ; do
 if [[ -d $sim ]]; then
  simulation=$(basename "${sim}")
  cd "${sim}"
  cpptraj -i prep_*
  sleep 1
  echo "${simulation} has been post-processed"
  cd -
 fi
done
