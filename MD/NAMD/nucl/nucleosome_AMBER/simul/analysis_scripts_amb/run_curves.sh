#!/bin/bash
#This is specifically tailored for 147 bp in 1kx5, change last lines of 
#input file if it is different
#the first parameter is name without pdb

source ~/.profile
cd ../analysis_data/dna_params_cur
Cur+ <<!
 &inp file=$1, lis=$1, 
 lib=/Users/alexeyshaytan/soft/curves+/standard
 &end
2 1 -1 0 0
1:147
294:148
!



