#!/bin/bash
# vmdtpy -python -e int_get_ions_raw_cryst.vmdpy
# vmdtpy -python -e int_get_nucl_prot_raw_cryst.vmdpy
# vmdtpy -python -e int_get_prot_prot_raw_cryst.vmdpy
# vmdtpy -python -e int_get_wat_raw_cryst.vmdpy



vmdtpy -python -e int_get_ions_raw.vmdpy
vmdtpy -python -e int_get_nucl_prot_raw.vmdpy
vmdtpy -python -e int_get_wat_raw.vmdpy
vmdtpy -python -e int_get_prot_prot_raw.vmdpy
vmdtpy -python -e int_get_wat_dna_raw.vmdpy


