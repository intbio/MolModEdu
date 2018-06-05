#!/bin/bash
DIRS=( "/Users/alexeyshaytan/work/1kx5_simul/analysis_scripts/" \
	"/Users/alexeyshaytan/work/1kx5nt_simul/analysis_scripts/" \
	"/Volumes/MDBD/Dropbox/work/NIH/histones_work/6md_nucl_expl/6md_1kx5_tails/simul/analysis_scripts" \
	"/Volumes/MDBD/Dropbox/work/NIH/histones_work/6md_nucl_expl/6md_1kx5_notails/simul/analysis_scripts" )

for i in "${DIRS[@]}"
do
echo "Syncing to " $i
rsync -uvrlptgo --exclude render_scripts --exclude view_nucl_md.sh --exclude extr_nucl_solv_dcd_md.tcl --exclude extr_nucl_dcd_min_eq.tcl --exclude extr_nucl_dcd_md.tcl * $i
#rsync -uvrlptgo  --exclude view_nucl_md.sh --exclude extr_nucl_solv_dcd_md.tcl --exclude extr_nucl_dcd_min_eq.tcl --exclude extr_nucl_dcd_md.tcl * $i
done

