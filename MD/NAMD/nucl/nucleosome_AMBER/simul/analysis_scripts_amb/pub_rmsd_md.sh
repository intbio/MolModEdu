#!/bin/bash
#vmdt -e extr_nucl_dcd_md.tcl
vmdt -e pub_rmsd_md_cryst.tcl #get the RMSD relative to crystal
python2.7 myplot.py ../analysis_data/pub_rmsd_md_cryst.dat --x_r 0 1000 --max_y 5 --grid True
open ../analysis_data/pub_rmsd_md_cryst.png

