#we calculate parts of final RMSD

mol load psf ../analysis_data/only_nucl_init.psf

mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd  waitfor all

# We will see following components
# all 8 histones * (core-fold, alpha, all) * ( CA, all noh )
# dna * ( backbone, sugar, bases )
#
#
#
#
#
#

set alpha_core "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98)))"

set hlist [list "segname CHA" "segname CHE" "segname CHB" "segname CHF" "segname CHC" "segname CHG" "segname CHD" "segname CHH"]
set hname [list "H3_A" "H3_E" "H4_B" "H4_F" "H2A_C" "H2A_G" "H2B_D" "H2B_H"]

set sel_opt [list $alpha_core "helix" "all"]
set sel_opt_name [list "HFold" "helix" "all"]

set sel_opt2 [list "name CA" "all"]
set sel_opt2_name [list "CA" "all"]

set bases "nucleic and not backbone and noh and not name C1' C2' C3' O4'"


set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile [open ../analysis_data/parts_rmsd.dat w]	

puts $outfile "RMSD for parts of nucleosome"
puts $outfile "Calculated using VMD"
puts $outfile "Time, ps"
puts $outfile "RMSD, A"

set ts "Time"

foreach seg $hname {

foreach selopt $sel_opt_name {
foreach selopt2 $sel_opt2_name {
append ts "\t$seg-$selopt-$selopt2"

}}}

append ts "\tDNA-bases\tDNA-backbone"

puts  $outfile $ts

for { set i 1 } { $i<=$nframes } { incr i } {

set ds ""
set time [expr 0.1 * $i]
puts -nonewline $outfile [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]

foreach seg $hlist {

foreach selopt $sel_opt {
foreach selopt2 $sel_opt2 {

set sel0 [atomselect top "$seg and $selopt and $selopt2" frame 0]

set sel1 [atomselect top "$seg and $selopt and $selopt2" frame $i]

set d [measure rmsd $sel0 $sel1]

append ds "\t$d"

$sel0 delete
$sel1 delete

}

}

}



set sel0 [atomselect top "$bases" frame 0]
set sel1 [atomselect top "$bases" frame $i]

set d [measure rmsd $sel0 $sel1]
append ds "\t$d"

$sel0 delete
$sel1 delete

set sel0 [atomselect top "nucleic and backbone" frame 0]
set sel1 [atomselect top "nucleic and backbone" frame $i]

set d [measure rmsd $sel0 $sel1]
append ds "\t$d"

$sel0 delete
$sel1 delete

puts $outfile "$ds"
#puts "$ds"

}
#puts "$rA"
# rmsd calculation loop

#  set nframes [expr  [molinfo top get numframes] - 1 ]

#  set outfile [open ../analysis_data/chain_rmsd.dat w]	
# # #set sel1 [atomselect top "segname CHI and resid '$r1'" frame $i]



 close $outfile 


exit


