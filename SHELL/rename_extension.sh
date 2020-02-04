#!/bin/bash
# rename extension of MD trajectory from mdcrd to netcdf
output=$(pwd)
for sim in "${output}"/* ; do
 if [[ -d $sim ]]; then
  simulation2=$(basename "${sim}")
  simulation="${simulation2/${project}}"
  cd "${sim}"
  rename "s/mdcrd/netcdf/" *.mdcrd
  echo "${simulation} has been post-processed! WELL DONE !"
  cd -
 fi
done
