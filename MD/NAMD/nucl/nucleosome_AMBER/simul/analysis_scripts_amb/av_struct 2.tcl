#We generate the average structure

mol load psf ../analysis_data/only_nucl_init.psf
mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
#mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd waitfor all


set sel_all [atomselect top "all"]



set av_pos [measure avpos $sel_all first 2500 last 10000 step 1]


set frame_av_all [atomselect top "all" frame 1]	
$frame_av_all set {x y z} $av_pos 

#superimpose on crystal just in case
set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
set frame_ref [atomselect top $alpha_core_ca frame 0]	 
set sel [atomselect top $alpha_core_ca frame 1]	
$frame_av_all move [measure fit $sel $frame_ref]

$frame_av_all writepdb ../analysis_data/only_nucl_average.pdb

exit


