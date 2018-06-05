#!/bin/bash

#Comment:
# this file 
#runs minimization/equilibration on lomonosov
#

echo "Do not forget to copy namd2 executable in the directory!!!"

mkdir -p ~/_scratch/nucleosome_amber/1kx5nt_simul

cp -r ~/nucleosome_amber/1kx5nt_simul/* ~/_scratch/nucleosome_amber/1kx5nt_simul

cd ~/_scratch/nucleosome_amber/1kx5nt_simul

cd output

sbatch -n64 -p regular4 -o min_eq.log impi ./namd2 ../input/min_equil.conf





