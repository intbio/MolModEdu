#!/bin/bash
#vmdt -e extr_nucl_dcd_md.tcl
vmdt -e rmsd_md_cryst.tcl #get the RMSD relative to crystal
vmdt -e rmsd_md_av.tcl #get the RMSD relative to average structure, see inside to specify for what frames to calculate vaerage
python2.7 myplot.py ../analysis_data/rmsd_md_cryst.dat --x_r 0 1000
open ../analysis_data/rmsd_md_cryst.png

python2.7 myplot.py ../analysis_data/rmsd_md_av.dat
open ../analysis_data/rmsd_md_av.png