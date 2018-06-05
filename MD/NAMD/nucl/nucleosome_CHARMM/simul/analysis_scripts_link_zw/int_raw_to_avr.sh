#!/bin/bash
/usr/bin/R --vanilla --slave < int_raw_to_avr_prot_prot.r
/usr/bin/R --vanilla --slave < int_raw_to_avr_dna_prot.r

/usr/bin/R --vanilla --slave < int_raw_to_avr_wat.r
/usr/bin/R --vanilla --slave < int_raw_to_avr_wat_dna.r
/usr/bin/R --vanilla --slave < int_raw_to_avr_ions.r

