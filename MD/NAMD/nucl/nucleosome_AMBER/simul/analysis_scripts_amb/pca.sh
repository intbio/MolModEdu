#!/bin/bash
python2.7 pca_dna.py
python2.7 pca_hfolds.py

#vmd -e visual_pca.tcl
#python2.7 myplot_rmsf.py ../analysis_data/rmsf_chains.dat
#open ../analysis_data/rmsf_chains.png