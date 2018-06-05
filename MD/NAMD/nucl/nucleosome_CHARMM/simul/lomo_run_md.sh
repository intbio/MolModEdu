#!/bin/bash

#Comment:
# this file 
#runs md on lomonosov
#

mkdir -p ~/_scratch/nucleosome/1kx5_simul_link

cp -r ~/nucleosome/1kx5_simul_link/* ~/_scratch/nucleosome/1kx5_simul_link

cd ~/_scratch/nucleosome/1kx5_simul_link
mkdir -p output
cd output
cp ~/software/source/NAMD_2.9_Source_CUDA/Linux-x86_64-icc/namd2_cuda .

sbatch --reservation=shaitan_5 -n1280 -N160 -p gpu -o md.log impi ./namd2_cuda +idlepoll ../input/md.conf



#sbatch -n256 -p regular4 -o md.log impi ./namd2 ../input/md.conf





