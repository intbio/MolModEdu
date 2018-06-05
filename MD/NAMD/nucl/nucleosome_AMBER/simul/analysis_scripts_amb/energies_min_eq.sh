#!/bin/bash
python2.7 namdplot.py -o ../analysis_data/energies_min_eq.dat ../output/min_eq.log 
python2.7 myplot_ener.py -first 11001 ../analysis_data/energies_min_eq.dat
open ../analysis_data/energies_min_eq.png