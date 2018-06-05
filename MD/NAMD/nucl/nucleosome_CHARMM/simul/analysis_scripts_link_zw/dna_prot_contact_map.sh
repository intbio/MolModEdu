#!/bin/bash
mkdir -p dat

vmdt -e dna_prot_contact_map.tcl -args 1 -73 -53 dat/1.dat &
vmdt -e dna_prot_contact_map.tcl -args 0 -52 -32 dat/2.dat &
vmdt -e dna_prot_contact_map.tcl -args 0 -31 -11 dat/3.dat &
vmdt -e dna_prot_contact_map.tcl -args 0 -10 10 dat/4.dat &
vmdt -e dna_prot_contact_map.tcl -args 0 11 31 dat/5.dat &
vmdt -e dna_prot_contact_map.tcl -args 0 32 52 dat/6.dat &
vmdt -e dna_prot_contact_map.tcl -args 0 53 73 dat/7.dat &

wait # this potentially can loose data if the last process exits earlier than other
sleep 5

cat dat/1.dat dat/2.dat dat/3.dat dat/4.dat dat/5.dat dat/6.dat dat/7.dat > ../analysis_data/dna_prot_contact_map.dat

#rm -rf dat