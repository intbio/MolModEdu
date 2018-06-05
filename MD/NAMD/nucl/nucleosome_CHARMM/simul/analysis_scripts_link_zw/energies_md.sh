#!/bin/bash
python2.7 namdplot.py -o ../analysis_data/energies_md.dat ../output/md.log -s 5
python2.7 myplot_ener.py -w 1 ../analysis_data/energies_md.dat
open ../analysis_data/energies_md.png