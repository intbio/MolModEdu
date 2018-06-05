#!/bin/bash

#generate publication graphs

./pub_rmsd_md.sh
vmdtpy -e pub_prot_rmsd_sch.tcl
vmdtpy -e pub_prot_rmsd.tcl
vmdtpy -e pub_dna_rmsd.tcl

/usr/bin/R --vanilla --slave < pub_dna2_path.r
/usr/bin/R --vanilla --slave < pub_dna2_path_ensemble.r
/usr/bin/R --vanilla --slave < pub_dna_rmsd.r
/usr/bin/R --vanilla --slave < pub_dna_rmsf.r
/usr/bin/R --vanilla --slave < pub_prot_rmsf.r
/usr/bin/R --vanilla --slave < pub_prot_rmsd.r
/usr/bin/R --vanilla --slave < pub_prot_rmsd_sch.r

cd ../analysis_data

mkdir -p pub_figs

cp pub_dna_rmsd.png pub_figs/rNCPamb3-dna_rmsd_evol.png
cp pub_dna_rmsf.png pub_figs/rNCPamb8-dna_rmsf.png
cp pub_dna2_geom_avr_ideal.png pub_figs/rNCPamb7-dna_geom_avr.png
cp pub_dna2_geom_fluct_ideal.png pub_figs/rNCPamb6-dna_geom_fluct.png
cp pub_prot_rmsd_sch.png pub_figs/rNCPamb5-prot_rmsd_sch_evol.png
cp pub_prot_rmsd.png pub_figs/rNCPamb4-prot_rmsd_bb_evol.png
cp pub_rmsd_md_cryst.png pub_figs/rNCPamb2-rmsd_md_vs_cryst.png
cp pub_prot_rmsf.png pub_figs/rNCPamb9-prot_rmsf.png



