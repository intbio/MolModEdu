#!/bin/bash
#vmdt -e extr_nucl_dcd_md.tcl
vmdt -e rmsd_md.tcl
python2.7 myplot.py ../analysis_data/rmsd_md.dat
open ../analysis_data/rmsd_md.png

vmdt -e rmsd_md_core.tcl
python2.7 myplot.py ../analysis_data/rmsd_md_core.dat
open ../analysis_data/rmsd_md_core.png