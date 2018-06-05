#!/bin/bash
#This is for broken line approx DNA analysis
#first param is file to get pairs, second is the file to analyze and ref_frames
#here we do the trick
source ~/.profile
cd ../analysis_data/dna_params_bl
cp $1.pdb temp.pdb
find_pair temp.pdb temp.inp
cp $2.pdb temp.pdb
analyze temp.inp
mv ref_frames.dat $2.dat








