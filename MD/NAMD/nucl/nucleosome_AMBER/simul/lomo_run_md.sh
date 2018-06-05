#!/bin/bash

#Comment:
# this file 
#runs md on lomonosov
#



mkdir -p ~/_scratch/nucleosome_amber/1kx5nt_simul

cp -r ~/nucleosome_amber/1kx5nt_simul/* ~/_scratch/nucleosome_amber/1kx5nt_simul

cd ~/_scratch/nucleosome_amber/1kx5nt_simul

cd output

sbatch -n256 -p regular4 -o md.log impi ./namd2 ../input/md.conf





