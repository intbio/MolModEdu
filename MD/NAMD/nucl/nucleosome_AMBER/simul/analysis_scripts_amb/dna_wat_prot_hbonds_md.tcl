mol load psf ../input/1kx5nt_ready.psf
#mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol addfile ../output/md.dcd  waitfor all

set chi [atomselect top "segname CHI"]	 
set chj [atomselect top "segname CHJ"]	 
# rmsd calculation loop
package require hbonds

set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile [open ../analysis_data/dna_wat_prot_hbonds_md.dat w]	
#set sel1 [atomselect top "segname CHI and resid '$r1'" frame $i]

puts $outfile "Water mediated interactions between DNA and protein"
puts $outfile "Calculated using VMD 3.0 dist, 30 deg ang param"
puts $outfile "Time, ps"
puts $outfile "Num waters"
puts -nonewline $outfile "Time"
for { set r1 -73 } { $r1<=73 } { incr r1 } {
puts -nonewline $outfile "\t$r1"
}
puts -nonewline $outfile "\n"
for { set i 0 } { $i<=$nframes } { incr i } {
set time [expr 0.1 * $i]
puts -nonewline $outfile [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]
for { set r1 -73 } { $r1<=73 } { incr r1 } {
set sel1 [atomselect top "protein" frame $i]
set r2 [expr $r1 * (-1)]
set sel2 [atomselect top "(segname CHI and resid '$r1') or (segname CHJ and resid '$r2')" frame $i]
set sel3 [atomselect top "water" frame $i]
set hbonds1 [lindex [measure hbonds 3.0 30 $sel2 $sel3] 1]
set hbonds2 [lindex [measure hbonds 3.0 30 $sel3 $sel2] 2] 
puts $hbonds1
set sel4 [atomselect top "same residue as index 1000000 $hbonds1 $hbonds2" frame $i]
set hbonds1 [lindex [measure hbonds 3.0 30 $sel4 $sel1] 1]
set hbonds2 [lindex [measure hbonds 3.0 30 $sel1 $sel4] 1]
set nh1 [llength $hbonds1]
set nh2 [llength $hbonds2]
set nh [expr $nh1 + $nh2]
set delta [expr $nh1 - $nh2]
#puts [format "Check: %d" "$delta"]
puts -nonewline $outfile "\t$nh"

$sel2 delete
$sel1 delete
$sel3 delete
$sel4 delete
#hbonds -sel1 $sel1 -sel2 $sel2 -dist 3.0 -ang 30 -frames $i:$i -writefile yes -outfile ../analysis_data/dna_hbonds.dat
#-type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

}
puts -nonewline $outfile "\n"
}

close $outfile 

#hbonds -sel1 $chi -sel2 $chj -dist 3.0 -ang 30 -type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

exit


