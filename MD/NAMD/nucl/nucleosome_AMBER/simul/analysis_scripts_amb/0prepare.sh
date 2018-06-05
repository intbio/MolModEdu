#!/bin/bash
# vmdt -e create_aligned_pdb.tcl
vmdtpy -python -e create_aligned_pdb2.vmdpy
vmdtpy -e extr_nucl_dcd_min_eq.tcl
vmdtpy -e extr_nucl_dcd_md.tcl
vmdtpy -e extr_nucl_solv_dcd_md.tcl
