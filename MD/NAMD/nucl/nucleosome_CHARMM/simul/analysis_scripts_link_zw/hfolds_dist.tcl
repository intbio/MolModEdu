#here we calculate distances between different histone folds

mol load psf ../analysis_data/only_nucl_init.psf
mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd waitfor all

set outfile [open ../analysis_data/hfolds_dist.dat w]	 
set nf [molinfo top get numframes]	 

set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"

set H3H4_1s "((segname CHA and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93)))  and name CA"
set H3H4_2s "((segname CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93)))  and name CA"
set H2AH2B_1s "((segname CHC and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
set H2AH2B_2s "((segname CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"


set rf 0

set H3H4_1 [atomselect top $H3H4_1s]	 
set H3H4_2 [atomselect top $H3H4_2s]
set H2AH2B_1 [atomselect top $H2AH2B_1s]	 
set H2AH2B_2 [atomselect top $H2AH2B_2s]	 
# rmsd calculation loop

puts $outfile "Distance between histone dimers"
# puts $outfile "RMSD of nucleosome during MD simulations relative to 500 ns frame"
puts $outfile "Calculated using centers of mass of CA atoms of histone folds"
puts $outfile "Time, ns"
puts $outfile "Distance, A"
puts $outfile "Time\tH3H4-H3H4\tH2AH2B-H2AH2B\tH2AH2B-H3H4_1\tH2AH2B-H3H4_2"
for { set i 1 } { $i<$nf } { incr i } {	 

$H3H4_1 frame $i
$H3H4_2 frame $i
$H2AH2B_1 frame $i
$H2AH2B_2 frame $i

#if not crystal as reference let's realign
# $sel_all move [measure fit $sel_alpha_core_ca $frame0_alpha_core_ca]
# $sel_all move [measure fit $sel_dna_bb $frame0_dna_bb]


set H3H4_1_p [measure center $H3H4_1 weight mass]
set H3H4_2_p [measure center $H3H4_2 weight mass]
set H2AH2B_1_p [measure center $H2AH2B_1 weight mass]
set H2AH2B_2_p [measure center $H2AH2B_2 weight mass]

set H3H4_H3H4 [veclength [vecsub $H3H4_1_p $H3H4_2_p]]
set H2AH2B_H2AH2B [veclength [vecsub $H2AH2B_1_p $H2AH2B_2_p]]
set H2AH2B_H3H4_1 [veclength [vecsub $H2AH2B_1_p $H3H4_1_p]]
set H2AH2B_H3H4_2 [veclength [vecsub $H2AH2B_2_p $H3H4_2_p]]



set time [expr ( $i - 1 ) * 0.1]
puts $outfile "$time\t$H3H4_H3H4\t$H2AH2B_H2AH2B\t$H2AH2B_H3H4_1\t$H2AH2B_H3H4_2"
 
}	 
close $outfile 

exec python2.7 myplot.py ../analysis_data/hfolds_dist.dat

exit


