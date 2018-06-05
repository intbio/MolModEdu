mol load psf ../input/1kx5_ready.psf
mol addfile ../analysis_data/system_aligned.pdb waitfor all
mol addfile ../output/min_eq.dcd waitfor all
#mol addfile md.dcd step 10 waitfor all

# we need to set secondary structure at first and correct it a bit, so that it would be symmtric
# we will choose to stick with stride results for first chains
set alpha_core_bb "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA O C N"
set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"

set dna_bb "segname CHI CHJ and name P O5' C5' C4' C3' O3' and resid '-73' to 73"
set dna_p "segname CHI CHJ and name P and resid '-73' to 73"

#set p [atomselect top "	"]
set alphabb0 [atomselect top $alpha_core_bb frame 0]
set alphaca0 [atomselect top $alpha_core_ca frame 0]

set dnabb0 [atomselect top $dna_bb frame 0]
set dnap0 [atomselect top $dna_p frame 0]


set nucl [atomselect top "protein or nucleic"]
set nucl0 [atomselect top "protein or nucleic" frame 0]

set frame0 [atomselect top "backbone and protein and alpha_helix and name CA" frame 0]	 
set sel [atomselect top "backbone and protein and alpha_helix and name CA"]	 

$nucl0 writepsf ../analysis_data/only_nucl_init.psf
$nucl0 writepdb ../analysis_data/only_nucl_init.pdb

$alphabb0 writepsf ../analysis_data/only_alpha_core_bb.psf
$alphaca0 writepsf ../analysis_data/only_alpha_core_ca.psf

$dnap0 writepsf ../analysis_data/only_dna_p.psf
$dnabb0 writepsf ../analysis_data/only_dna_bb.psf

$dnap0 writepdb ../analysis_data/only_dna_p_init.pdb
$dnabb0 writepdb ../analysis_data/only_dna_bb_init.pdb


set nframes [expr  [molinfo top get numframes] - 1 ]

for { set i 1 } { $i<=$nframes } { incr i } {	 
$sel frame $i
$nucl frame $i
$nucl move [measure fit $sel $frame0]
 
}


animate write dcd ../analysis_data/min_eq_nucl.dcd beg 1 end $nframes waitfor all sel $nucl


exit


