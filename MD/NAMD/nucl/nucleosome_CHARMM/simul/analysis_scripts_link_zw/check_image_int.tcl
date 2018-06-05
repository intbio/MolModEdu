# This sctipt should analyze minimal distance between nucleosome and its image.

#seems impossible to implement easily

package require pbctools

mol load psf ../input/1kx5nt_ready.psf
mol addfile ../output/md.dcd last 1 waitfor all
mol ssrecalc top
#mol addfile ../analysis_data/md_nucl.dcd waitfor all

set outfile [open ../analysis_data/check_image_int.dat w]	 
set nf [molinfo top get numframes]
puts "Numframes $nf" 

set sel [atomselect top "protein or nucleic and (name CA P)"]

puts "Atom number [$sel num]"
set min 100

for { set i 1 } { $i<$nf } { incr i } {

set pbc [pbc get]
set pbcx [lindex $pbc 0]
set pbcy [lindex $pbc 1]
set pbcz [lindex $pbc 2]
$sel frame $i



foreach v1 [$sel get {x y z}] {

foreach v2 [$sel get {x y z}] {

set d [vecdist $v1 $v2]

if {$min > $d} {
set min $d
}

}
}

}

puts "Minimum is $min"



# rmsd calculation loop

# puts $outfile "RMSD of nucleosome during MD simulations"
# puts $outfile "Align using helix CA structures"
# puts $outfile "Time, ns"
# puts $outfile "RMSD, A"
# puts $outfile "Time\tAll\t$\\alpha$-helices\t$\\alpha$-helix C$\\alpha$\t$\\alpha$-helix sch\tDNA\tDNA bb\t$\\alpha$-hHfolds C$\\alpha$"
# for { set i 1 } { $i<$nf } { incr i } {	 
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

# set time [expr ( $i - 1 ) * 0.1]
# puts $outfile "$time\t$rmsd_all\t$rmsd_helix\t$rmsd_helix_CA\t$rmsd_helix_sidech\t$rmsd_dna\t$rmsd_dna_bb\t$rmsd_alpha_core_ca"
 
# }	 
close $outfile 


#exit


