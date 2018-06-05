#this script will put ions in the box back,
#fit nuclesome to initial position and cut ions and waters 1 nm apart from nucl

mol load psf ../input/1kx5_ready.psf
mol addfile ../analysis_data/system_aligned.pdb waitfor all
mol addfile ../output/md.dcd step 10 waitfor all
mol addfile ../output/md2.dcd step 10 waitfor all

#mol addfile md.dcd step 10 waitfor all

package require pbctools
set nframes [expr  [molinfo top get numframes] - 1 ]
#set ions [atomselect top "resname SOD"]

### The following must be simplified in a case where nucleosome is fixed!!!
### we can delete water, and set center to origin!!!
pbc wrap -all -center com -compound res  -centersel "protein" -sel "resname SOD or resname CLA or resname TIP3"

set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"


set frame0_nucl_alpha [atomselect top $alpha_core_ca frame 0]	 
set nucl_alpha [atomselect top $alpha_core_ca]

set all [atomselect top "all"]

#-1 because the frames are numbered from 0
set nframes [expr  [molinfo top get numframes] - 1 ]

for { set i 1 } { $i<=$nframes } { incr i } {	 
$nucl_alpha frame $i
$all frame $i
$all move [measure fit $nucl_alpha $frame0_nucl_alpha]
 
}




animate write dcd ../analysis_data/md_nucl_solv.dcd beg 1 end $nframes waitfor all sel $all


exit


