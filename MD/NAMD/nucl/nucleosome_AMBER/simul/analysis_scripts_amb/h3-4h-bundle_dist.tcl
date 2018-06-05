#here we calculate distances between different histone folds

mol load psf ../analysis_data/only_nucl_init.psf
mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd waitfor all

set outfile [open ../analysis_data/h3-4h-bundle.dat w]	 
set nf [molinfo top get numframes]	 

set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"

set H3A_4hs "((segname CHA and (resid 104 to 114 or resid 121 to 131)) and name CA)"
set H3E_4hs "((segname CHE and (resid 104 to 114 or resid 121 to 131)) and name CA)"


set rf 0

set H3A_4h [atomselect top $H3A_4hs]	 
set H3E_4h [atomselect top $H3E_4hs]
# rmsd calculation loop

puts $outfile "Distance between H3 4-helical bundle parts (COM of C-alphas)"
# puts $outfile "RMSD of nucleosome during MD simulations relative to 500 ns frame"
puts $outfile "Calculated using centers of mass of CA atoms"
puts $outfile "Time, ns"
puts $outfile "Distance, A"
puts $outfile "Time\tH3A_4h-H3E_4h"
for { set i 1 } { $i<$nf } { incr i } {	 

$H3A_4h frame $i
$H3E_4h frame $i

#if not crystal as reference let's realign
# $sel_all move [measure fit $sel_alpha_core_ca $frame0_alpha_core_ca]
# $sel_all move [measure fit $sel_dna_bb $frame0_dna_bb]


set H3A_4h_p [measure center $H3A_4h weight mass]
set H3E_4h_p [measure center $H3E_4h weight mass]

set dist [veclength [vecsub $H3A_4h_p $H3E_4h_p]]



set time [expr ( $i - 1 ) * 0.1]
puts $outfile "$time\t$dist"
 
}	 
close $outfile 

exec python2.7 myplot.py ../analysis_data/h3-4h-bundle.dat

exit


