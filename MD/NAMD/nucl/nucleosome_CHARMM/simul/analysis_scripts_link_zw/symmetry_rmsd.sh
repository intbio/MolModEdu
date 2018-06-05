#!/bin/bash
vmdt -e symmetry_rmsd.tcl
python2.7 myplot.py ../analysis_data/symmetry_rmsd.dat
open ../analysis_data/symmetry_rmsd.png