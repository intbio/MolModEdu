#!/bin/bash
vmdt -e dna_hbonds_md.tcl
/usr/bin/R --vanilla --slave < dna_hbonds_md.r
open ../analysis_data/dna_hbonds_hist_md.png