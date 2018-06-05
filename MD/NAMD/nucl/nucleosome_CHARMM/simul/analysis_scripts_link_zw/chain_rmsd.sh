#!/bin/bash
vmdt -e chain_rmsd.tcl
/usr/bin/R --vanilla --slave < chain_rmsd.r
open ../analysis_data/chain_rmsd.png