mol load psf ../analysis_data/only_nucl_init.psf
# mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd step 10 waitfor all
mol addfile ../analysis_data/only_nucl_average.pdb

set outfile [open ../analysis_data/cluster_md_hfolds.dat w]	 
set nf [molinfo top get numframes]


set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"

set sel_all [atomselect top "all and noh"]


set ff [expr $nf / 2]
set lf [expr $nf - 2]
set rf [expr $nf - 1]

set ff 1



# set av_pos [measure avpos $sel_all first $ff last $lf step 1]
set sel_hfolds [atomselect top "$alpha_core_ca and noh"]
set sel_hfolds_ref [atomselect top "$alpha_core_ca and noh" frame $rf]

# #Let's optionally realign with average structure
for { set i 0 } { $i<$lf } { incr i } {	 
$sel_all frame $i
$sel_hfolds frame $i
$sel_all move [measure fit $sel_hfolds $sel_hfolds_ref]
 
}

set clust_res [measure cluster $sel_hfolds num 5 distfunc rmsd cutoff 0.8 first $ff last $lf step 1 ]
#set clust_res [measure cluster $sel_hfolds num 5 distfunc fitrmsd cutoff 0.5 first $ff last $lf step 1 ]
# it seems there is currently a bug in fitrmsd

set nc [expr [llength $clust_res] - 1]

for { set i 0 } { $i<$nc } { incr i } {	
puts $outfile [lindex $clust_res $i]

}
close $outfile 



set outfile [open ../analysis_data/cluster_md_dnabb.dat w]	 




# set av_pos [measure avpos $sel_all first $ff last $lf step 1]
set sel_dnabb [atomselect top "nucleic and backbone and noh"]
set clust_res [measure cluster $sel_dnabb num 5 distfunc rmsd cutoff 3 first $ff last $lf step 1 ]

set nc [expr [llength $clust_res] - 1]

for { set i 0 } { $i<$nc } { incr i } {	
puts $outfile [lindex $clust_res $i]

}
close $outfile 

# set frame0_all [atomselect top "all and noh" frame $nf]
# $frame0_all set {x y z} $av_pos 
# set sel_all [atomselect top "all and noh"]

# set frame0_helix [atomselect top "alpha_helix and noh" frame $nf] 
# set sel_helix [atomselect top "alpha_helix and noh"]

# set frame0_helix_CA [atomselect top "alpha_helix and name CA" frame $nf]
# set sel_helix_CA [atomselect top "alpha_helix and name CA"]

# set frame0_helix_sidech [atomselect top "sidechain and alpha_helix and noh" frame $nf]
# set sel_helix_sidech [atomselect top "sidechain and alpha_helix and noh"]

# set frame0_dna [atomselect top "nucleic and noh" frame $nf]
# set sel_dna [atomselect top "nucleic and noh"]

# set frame0_dna_bb [atomselect top "backbone and nucleic and noh" frame $nf] 
# set sel_dna_bb [atomselect top "backbone and nucleic and noh"]

# set frame0_alpha_core_ca [atomselect top $alpha_core_ca frame $nf]
# set sel_alpha_core_ca [atomselect top $alpha_core_ca]
# # rmsd calculation loop

# puts $outfile "RMSD of nucleosome during MD simulations relative to average structure"
# puts $outfile "Align using helix CA structures"
# puts $outfile "Time, ns"
# puts $outfile "RMSD, A"
# puts $outfile "Time\tAll\t$\\alpha$-helices\t$\\alpha$-helix C$\\alpha$\t$\\alpha$-helix sch\tDNA\tDNA bb\t$\\alpha$-hHfolds C$\\alpha$"
# for { set i 0 } { $i<$nf } { incr i } {	 
# $sel_all frame $i
# $sel_helix frame $i
# $sel_helix_CA frame $i
# $sel_helix_sidech frame $i
# $sel_dna frame $i
# $sel_dna_bb frame $i
# $sel_alpha_core_ca frame $i

# set rmsd_all [measure rmsd $sel_all $frame0_all]
# set rmsd_helix [measure rmsd $sel_helix $frame0_helix]
# set rmsd_helix_CA [measure rmsd $sel_helix_CA $frame0_helix_CA]
# set rmsd_helix_sidech [measure rmsd $sel_helix_sidech $frame0_helix_sidech]
# set rmsd_dna [measure rmsd $sel_dna $frame0_dna]
# set rmsd_dna_bb [measure rmsd $sel_dna_bb $frame0_dna_bb]

# set rmsd_alpha_core_ca [measure rmsd $sel_alpha_core_ca $frame0_alpha_core_ca]

# set time [expr ( $i ) * 0.1]
# puts $outfile "$time\t$rmsd_all\t$rmsd_helix\t$rmsd_helix_CA\t$rmsd_helix_sidech\t$rmsd_dna\t$rmsd_dna_bb\t$rmsd_alpha_core_ca"
 
# }	 
# close $outfile 

exit


