#!/bin/bash
vmdtpy -e rmsf.tcl
vmdpy -e visual_rmsf.tcl # Adjust scale minmax values here to get best view
python2.7 myplot_rmsf.py ../analysis_data/rmsf_chains.dat --max_y 2 2 2 2 5 --legend "upper center"
python2.7 myplot_rmsf.py ../analysis_data/rmsf_chains_sch.dat --max_y 3.5 3.5 3.5 3.5 3.5 --legend "upper center"

##### We need to compare our RMSF plot with experiment
#
# B-factor = 8*pi^2*RMSF^2
#
python2.7 myplot_rmsf_exp.py ../analysis_data/rmsf_chains.dat --max_y 1.6 1.6 1.6 1.6 4.7 --legend "lower left"
python2.7 myplot_rmsf_exp_sch.py ../analysis_data/rmsf_chains_sch.dat --max_y 2.5 2.5 3.5 2.5 3.7 --legend "lower left"
python2.7 myplot_rmsf_exp_dna.py ../analysis_data/rmsf_chains.dat ../analysis_data/rmsf_chains_sch.dat --max_y 4.5 3.0  --legend "upper center"

open ../analysis_data/rmsf_chains.png


#Now we do all the same but use crystal structure as a reference
#This was found not to be intresting. crystal is far from average.

# vmdt -e rmsf_cryst.tcl
# vmd -e visual_rmsf_cryst.tcl
# python2.7 myplot_rmsf.py ../analysis_data/rmsf_cryst_chains.dat
# python2.7 myplot_rmsf_exp.py ../analysis_data/rmsf_cryst_chains.dat
