#!/bin/bash
# this script extract info from logs of MD simulations
# it could be adapted to extract any info from MD output in txt format
output=$(pwd)

for sim in "${output}"/* ; do
 if [[ -d $sim ]]; then
  simulation=$(basename "${sim}")
  cd "${sim}"
  # obtain minimized step in pdb format
  ambpdb -p ${sim}/protein.prmtop -c ${sim}/mini_finish*.rst > ${sim}/${simulation}_after_minimization.pdb 
  rm ${sim}/run_job_cuda.sh
  rm ${sim}/slurm*.out
  echo "${simulation}:" >> ${output}/final_energies.txt
  echo '' >> ${output}/final_energies.txt
  echo '             IN A MOMENT JuST BEFORE 4 ROUNDS OF ENERGY MINIMIZATION' >> ${output}/final_energies.txt
  echo '' >> ${output}/final_energies.txt
  echo '' >> ${output}/final_energies.txt
  sed -n 250,255p ${sim}/temp_amber/mini25p.out >> ${output}/final_energies.txt
  echo '' >> ${output}/final_energies.txt
  echo '' >> ${output}/final_energies.txt
  grep -A9 FINAL ${sim}/mini_finish*.out >> ${output}/final_energies.txt
  #sed -n 3465,3470p ${sim}/mini_finish*.out >> ${output}/final_energies.txt
  echo '' >> ${output}/final_energies.txt
  echo '********' >> ${output}/final_energies.txt
  echo '' >> ${output}/final_energies.txt
  rm -r ${sim}/temp_amber
  echo "${simulation} has been post-processed"
  cd -
 fi
done
