# mol load psf ../analysis_data/only_nucl_init.psf
mol load pdb ../analysis_data/only_nucl_init.pdb 
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd waitfor all

set outfile [open ../analysis_data/rmsd_md_cryst.dat w]	 
set nf [molinfo top get numframes]	 

set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"

set rf 0

set frame0_all [atomselect top "all and noh" frame $rf]	 
set sel_all [atomselect top "all and noh"]	

set frame0_helix [atomselect top "alpha_helix and noh" frame $rf]	 
set sel_helix [atomselect top "alpha_helix and noh"]	

set frame0_helix_CA [atomselect top "alpha_helix and name CA" frame $rf]	 
set sel_helix_CA [atomselect top "alpha_helix and name CA"]	

set frame0_helix_sidech [atomselect top "sidechain and alpha_helix and noh" frame $rf]	 
set sel_helix_sidech [atomselect top "sidechain and alpha_helix and noh"]	

set frame0_dna [atomselect top "nucleic and noh" frame $rf]	 
set sel_dna [atomselect top "nucleic and noh"]	

set frame0_dna_bb [atomselect top "backbone and nucleic and noh" frame $rf]	 
set sel_dna_bb [atomselect top "backbone and nucleic and noh"]

set frame0_alpha_core_ca [atomselect top $alpha_core_ca frame $rf]	 
set sel_alpha_core_ca [atomselect top $alpha_core_ca]	 
# rmsd calculation loop

puts $outfile "RMSD of nucleosome during MD simulations relative to crystal stucture"
# puts $outfile "RMSD of nucleosome during MD simulations relative to 500 ns frame"
puts $outfile "Align using helix CA structures"
puts $outfile "Time, ns"
puts $outfile "RMSD, A"
puts $outfile "Time\tAll\t$\\alpha$-helices\t$\\alpha$-helix C$\\alpha$\t$\\alpha$-helix sch\tDNA\tDNA bb\t$\\alpha$-hHfolds C$\\alpha$"
for { set i 0 } { $i<$nf } { incr i } {	 
$sel_all frame $i
$sel_helix frame $i
$sel_helix_CA frame $i
$sel_helix_sidech frame $i
$sel_dna frame $i
$sel_dna_bb frame $i
$sel_alpha_core_ca frame $i

#if not crystal as reference let's realign
# $sel_all move [measure fit $sel_alpha_core_ca $frame0_alpha_core_ca]
# $sel_all move [measure fit $sel_dna_bb $frame0_dna_bb]


set rmsd_all [measure rmsd $sel_all $frame0_all]
set rmsd_helix [measure rmsd $sel_helix $frame0_helix]
set rmsd_helix_CA [measure rmsd $sel_helix_CA $frame0_helix_CA]
set rmsd_helix_sidech [measure rmsd $sel_helix_sidech $frame0_helix_sidech]
set rmsd_dna [measure rmsd $sel_dna $frame0_dna]
set rmsd_dna_bb [measure rmsd $sel_dna_bb $frame0_dna_bb]

set rmsd_alpha_core_ca [measure rmsd $sel_alpha_core_ca $frame0_alpha_core_ca]

set time [expr ( $i - 1 ) * 0.1]
puts $outfile "$time\t$rmsd_all\t$rmsd_helix\t$rmsd_helix_CA\t$rmsd_helix_sidech\t$rmsd_dna\t$rmsd_dna_bb\t$rmsd_alpha_core_ca"
 
}	 
close $outfile 

exit


