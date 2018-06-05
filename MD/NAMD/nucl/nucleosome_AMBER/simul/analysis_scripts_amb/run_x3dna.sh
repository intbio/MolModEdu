#!/bin/bash
#first param is file to get pairs, second is to analyze
#here we do the trick
source ~/.profile
cd ../analysis_data/dna_params
cp $1.pdb temp.pdb
find_pair temp.pdb temp.inp
cp $2.pdb temp.pdb
analyze temp.inp
mv bp_step.par $2.dat
find_pair $2.pdb $2.inp
analyze -t=$2.tor $2.pdb







