#!/bin/bash
vmdt -e parts_rmsd.tcl
/usr/bin/R --vanilla --slave < parts_rmsd.r
open ../analysis_data/parts_rmsd.png