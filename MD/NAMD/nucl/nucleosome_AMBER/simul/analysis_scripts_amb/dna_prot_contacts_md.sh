#!/bin/bash
vmdt -e dna_prot_contacts_md.tcl
/usr/bin/R --vanilla --slave < dna_prot_contacts_md.r
open ../analysis_data/dna_prot_contacts_hist_md.png