#!/bin/bash
vmdt -e dna_prot_hbonds_md.tcl
/usr/bin/R --vanilla --slave < dna_prot_hbonds_md.r
open ../analysis_data/dna_prot_hbonds_hist_md.png