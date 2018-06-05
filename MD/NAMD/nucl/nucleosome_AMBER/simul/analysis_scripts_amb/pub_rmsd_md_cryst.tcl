# mol load psf ../analysis_data/only_nucl_init.psf
mol load pdb ../analysis_data/only_nucl_init.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd waitfor all

set outfile [open ../analysis_data/pub_rmsd_md_cryst.dat w]	 
set nf [molinfo top get numframes]	 

source def_struct_param.tcl


set rf 0

set frame0_all [atomselect top "all and noh" frame $rf]	 
set sel_all [atomselect top "all and noh"]	


set frame0_helix_CA [atomselect top $alpha_ext_ca frame $rf]	 
set sel_helix_CA [atomselect top $alpha_ext_ca]	

set frame0_dna_n1_9 [atomselect top "nucleic and name N1 N9" frame $rf]	 
set sel_dna_n1_9 [atomselect top "nucleic and name N1 N9"]

set frame0_alpha_hfolds_ca [atomselect top $alpha_hfolds_ca frame $rf]	 
set sel_alpha_hfolds_ca [atomselect top $alpha_hfolds_ca]	 

set frame0_core_ca [atomselect top $core_ca frame $rf]
set sel_core_ca [atomselect top $core_ca]
# rmsd calculation loop

puts $outfile "RMSD of nucleosome during MD simulations relative to crystal stucture"
# puts $outfile "RMSD of nucleosome during MD simulations relative to 500 ns frame"
puts $outfile "Align using helix CA structures"
puts $outfile "Time, ns"
puts $outfile "RMSD, A"
puts $outfile "Time\tHistone fold helices C$\\alpha$\tCore C$\\alpha$\tDNA N1,N9"

# puts $outfile "Time\t$\\alpha$-helix C$\\alpha$\tHistone fold helices C$\\alpha$\tCore C$\\alpha$\tDNA N1"
for { set i 0 } { $i<$nf } { incr i } {	 
# $sel_all frame $i
# $sel_helix frame $i
$sel_helix_CA frame $i
# $sel_helix_sidech frame $i
# $sel_dna frame $i
$sel_dna_n1_9 frame $i
$sel_alpha_hfolds_ca frame $i

$sel_core_ca frame $i

#if not crystal as reference let's realign
# $sel_all move [measure fit $sel_alpha_core_ca $frame0_alpha_core_ca]
# $sel_all move [measure fit $sel_dna_bb $frame0_dna_bb]


# set rmsd_all [measure rmsd $sel_all $frame0_all]
# set rmsd_helix [measure rmsd $sel_helix $frame0_helix]
set rmsd_helix_CA [measure rmsd $sel_helix_CA $frame0_helix_CA]
# set rmsd_helix_sidech [measure rmsd $sel_helix_sidech $frame0_helix_sidech]
# set rmsd_dna [measure rmsd $sel_dna $frame0_dna]
set rmsd_dna_n1_9 [measure rmsd $sel_dna_n1_9 $frame0_dna_n1_9]

set rmsd_alpha_hfolds_ca [measure rmsd $sel_alpha_hfolds_ca $frame0_alpha_hfolds_ca]

set rmsd_core_ca [measure rmsd $sel_core_ca $frame0_core_ca]


set time [expr ( $i - 1 ) * 0.1]
# puts $outfile "$time\t$rmsd_helix_CA\t$rmsd_alpha_core_ca\t$rmsd_core_ca\t$rmsd_dna_n1"
puts $outfile "$time\t$rmsd_alpha_hfolds_ca\t$rmsd_core_ca\t$rmsd_dna_n1_9"
 
}	 
close $outfile 

exit


