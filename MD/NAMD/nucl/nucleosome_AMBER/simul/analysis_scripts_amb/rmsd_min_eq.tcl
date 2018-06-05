# mol load psf ../analysis_data/only_nucl_init.psf
mol load pdb ../analysis_data/only_nucl_init.pdb 
mol addfile ../analysis_data/min_eq_nucl.dcd waitfor all

set outfile [open ../analysis_data/rmsd_min_eq.dat w]	 
set nf [molinfo top get numframes]	 

set frame0_all [atomselect top "all and noh" frame 0]	 
set sel_all [atomselect top "all and noh"]	

set frame0_helix [atomselect top "alpha_helix and noh" frame 0]	 
set sel_helix [atomselect top "alpha_helix and noh"]	

set frame0_helix_CA [atomselect top "alpha_helix and name CA" frame 0]	 
set sel_helix_CA [atomselect top "alpha_helix and name CA"]	

set frame0_helix_sidech [atomselect top "sidechain and alpha_helix and noh" frame 0]	 
set sel_helix_sidech [atomselect top "sidechain and alpha_helix and noh"]	

set frame0_dna [atomselect top "nucleic and noh" frame 0]	 
set sel_dna [atomselect top "nucleic and noh"]	

set frame0_dna_bb [atomselect top "backbone and nucleic and noh" frame 0]	 
set sel_dna_bb [atomselect top "backbone and nucleic and noh"]	 
# rmsd calculation loop

puts $outfile "RMSD of nucleosome during minimization and equilibration"
puts $outfile "Align using helix CA structures"
puts $outfile "Frames"
puts $outfile "RMSD, A"
puts $outfile "Frame\tAll\t$\\alpha$-helices\t$\\alpha$-helix C$\\alpha$\t$\\alpha$-helix sidechains\tDNA\tDNA backbone"
for { set i 1 } { $i<$nf } { incr i } {	 
$sel_all frame $i
$sel_helix frame $i
$sel_helix_CA frame $i
$sel_helix_sidech frame $i
$sel_dna frame $i
$sel_dna_bb frame $i

set rmsd_all [measure rmsd $sel_all $frame0_all]
set rmsd_helix [measure rmsd $sel_helix $frame0_helix]
set rmsd_helix_CA [measure rmsd $sel_helix_CA $frame0_helix_CA]
set rmsd_helix_sidech [measure rmsd $sel_helix_sidech $frame0_helix_sidech]
set rmsd_dna [measure rmsd $sel_dna $frame0_dna]
set rmsd_dna_bb [measure rmsd $sel_dna_bb $frame0_dna_bb]

puts $outfile "$i\t$rmsd_all\t$rmsd_helix\t$rmsd_helix_CA\t$rmsd_helix_sidech\t$rmsd_dna\t$rmsd_dna_bb"
 
}	 
close $outfile 

exit


