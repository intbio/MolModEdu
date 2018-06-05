#!/bin/bash
#Averaging
/usr/bin/R --vanilla --slave < int_raw_to_avr_dna_prot.r
/usr/bin/R --vanilla --slave < int_raw_to_avr_prot_prot.r
/usr/bin/R --vanilla --slave < int_raw_to_avr_wat.r
/usr/bin/R --vanilla --slave < int_raw_to_avr_wat_dna.r

/usr/bin/R --vanilla --slave < int_raw_to_avr_ions.r


#Analysis and plotting

/usr/bin/R --vanilla --slave < int_total_stat.r
/usr/bin/R --vanilla --slave < int_dna_prot.r
/usr/bin/R --vanilla --slave < int_dna_prot_var.r
/usr/bin/R --vanilla --slave < int_prot_prot.r
/usr/bin/R --vanilla --slave < int_wat.r
/usr/bin/R --vanilla --slave < int_wat_dna.r
/usr/bin/R --vanilla --slave < int_ions.r






