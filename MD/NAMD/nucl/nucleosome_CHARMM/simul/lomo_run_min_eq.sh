#!/bin/bash

#Comment:
# this file 
#runs minimization/equilibration on lomonosov
#

echo "Do not forget to copy namd2 executable in the directory!!!"

mkdir -p ~/_scratch/nucleosome/1kx5_simul_link

cp -r ~/nucleosome/1kx5_simul_link/* ~/_scratch/nucleosome/1kx5_simul_link

cd ~/_scratch/nucleosome/1kx5_simul_link
mkdir -p output
cd output
cp ~/software/source/NAMD_2.9_Source_CUDA/Linux-x86_64-icc/namd2_cuda .

sbatch --reservation=shaitan_5 -n1280 -N160 -p gpu -o min_eq.log impi ./namd2_cuda +idlepoll ../input/min_equil.conf





