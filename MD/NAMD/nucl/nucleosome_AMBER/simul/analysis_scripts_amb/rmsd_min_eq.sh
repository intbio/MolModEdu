#!/bin/bash
vmdt -e extr_nucl_dcd_min_eq.tcl
vmdt -e rmsd_min_eq.tcl
python2.7 myplot.py ../analysis_data/rmsd_min_eq.dat
open ../analysis_data/rmsd_min_eq.png