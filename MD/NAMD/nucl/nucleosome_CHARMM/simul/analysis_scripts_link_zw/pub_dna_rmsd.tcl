#we need to make a heat plot of RMSD evolution of certain bases
#
# 
#
#



mol load psf ../analysis_data/only_nucl_init.psf
mol addfile ../analysis_data/only_nucl_aligned.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd step 10 first 1 waitfor all


set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile [open ../analysis_data/pub_dna_rmsd.dat w]	


set chi [atomselect top "segname CHI"]	 
set chj [atomselect top "segname CHJ"]	 
# rmsd calculation loop


puts $outfile "DNA RMSD profile"
puts $outfile "Calculated using VMD"
puts $outfile "Time, ps"
puts $outfile "RMSD, A"
puts -nonewline $outfile "Time"
for { set r1 -73 } { $r1<=73 } { incr r1 } {
puts -nonewline $outfile "\t$r1"
}
puts -nonewline $outfile "\n"

for { set i 1 } { $i<=$nframes } { incr i } {
set time [expr 1.0 * $i]
puts -nonewline $outfile [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]

for { set r1 -73 } { $r1<=73 } { incr r1 } {

set r2 [expr $r1 * (-1)]
set sel2 [atomselect top "((segname CHI and resid '$r1' and name N1 N9) or (segname CHJ and resid '$r2' and name N1 N9)) and noh" frame $i]
set sel0 [atomselect top "((segname CHI and resid '$r1' and name N1 N9) or (segname CHJ and resid '$r2' and name N1 N9)) and noh" frame 0]

set rmsd [measure rmsd $sel2 $sel0]


puts -nonewline $outfile "\t$rmsd"

$sel2 delete
$sel0 delete

}
puts -nonewline $outfile "\n"
}

close $outfile 


exit






