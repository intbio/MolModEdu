#!/bin/bash
vmdtpy -python -e int_get_raw_prot_prot_cryst.vmdpy
vmdtpy -python -e int_get_raw_prot_prot.vmdpy

vmdtpy -python -e int_get_raw_nucl_prot_cryst.vmdpy
vmdtpy -python -e int_get_raw_nucl_prot.vmdpy


vmdtpy -python -e int_get_raw_ions_cryst.vmdpy
vmdtpy -python -e int_get_raw_wat_cryst.vmdpy
vmdtpy -python -e int_get_raw_wat_dna_cryst.vmdpy


vmdtpy -python -e int_get_raw_ions.vmdpy
vmdtpy -python -e int_get_raw_wat.vmdpy
vmdtpy -python -e int_get_raw_wat_dna.vmdpy




