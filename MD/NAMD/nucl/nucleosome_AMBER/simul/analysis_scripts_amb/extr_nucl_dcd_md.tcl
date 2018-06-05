# mol load psf ../input/1kx5nt_ready.psf
mol load pdb ../analysis_data/system_aligned.pdb 
mol addfile ../output/md.dcd waitfor all
#mol addfile md.dcd step 10 waitfor all
source def_struct_param.tcl

# set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"

set nucl [atomselect top "protein or nucleic"]
set nucl0 [atomselect top "protein or nucleic" frame 0]

set frame0 [atomselect top $alpha_hfolds_ca frame 0]	 
set sel [atomselect top $alpha_hfolds_ca]	 

# $nucl0 writepsf ../analysis_data/only_nucl_init.psf
$nucl0 writepdb ../analysis_data/only_nucl_init.pdb

set nframes [expr  [molinfo top get numframes] - 1 ]

for { set i 1 } { $i<=$nframes } { incr i } {	 
$sel frame $i
$nucl frame $i
$nucl move [measure fit $sel $frame0]
 
}	 


animate write dcd ../analysis_data/md_nucl.dcd beg 1 end $nframes waitfor all sel $nucl


exit


