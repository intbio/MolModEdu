#we will calculate RMSD for histone alpha segments and for DNA

mol load psf ../analysis_data/only_nucl_init.psf

#mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol addfile ../analysis_data/md_nucl.dcd waitfor all


set h3_1 [atomselect top "segname CHA and alpha_helix and noh"] 
set h3_2 [atomselect top "segname CHE and alpha_helix and noh"] 

set h4_1 [atomselect top "segname CHB and alpha_helix and noh"] 
set h4_2 [atomselect top "segname CHF and alpha_helix and noh"] 

set h2a_1 [atomselect top "segname CHC and alpha_helix and noh"] 
set h2a_2 [atomselect top "segname CHG and alpha_helix and noh"] 

set h2b_1 [atomselect top "segname CHD and alpha_helix and noh"] 
set h2b_2 [atomselect top "segname CHH and alpha_helix and noh"] 

set dnai [atomselect top "segname CHI" frame 0] 
set dnaj [atomselect top "segname CHJ" frame 0]

set h3_1_0 [atomselect top "segname CHA and alpha_helix and noh" frame 0] 
set h3_2_0 [atomselect top "segname CHE and alpha_helix and noh" frame 0] 
set h4_1_0 [atomselect top "segname CHB and alpha_helix and noh" frame 0] 
set h4_2_0 [atomselect top "segname CHF and alpha_helix and noh" frame 0] 

set h2a_1_0 [atomselect top "segname CHC and alpha_helix and noh" frame 0] 
set h2a_2_0 [atomselect top "segname CHG and alpha_helix and noh" frame 0] 
set h2b_1_0 [atomselect top "segname CHD and alpha_helix and noh" frame 0] 
set h2b_2_0 [atomselect top "segname CHH and alpha_helix and noh" frame 0] 

set dnai_0 [atomselect top "segname CHI" frame 0] 
set dnaj_0 [atomselect top "segname CHJ" frame 0] 



#puts "$rA"
# rmsd calculation loop

 set nframes [expr  [molinfo top get numframes] - 1 ]

 set outfile [open ../analysis_data/chain_rmsd.dat w]	
# #set sel1 [atomselect top "segname CHI and resid '$r1'" frame $i]

 puts $outfile "RMSD for histone cores (all_alpha) and DNA"
puts $outfile "Calculated using VMD"
puts $outfile "Time, ps"
puts $outfile "RMSD, A"
puts  $outfile "Time\trH3_1\trH3_2\trH4_1\trH4_2\trH2A_1\trH2A_2\trH2B_1\trH2B_2\trDNA_I\trDNA_J"
# for { set r1 -73 } { $r1<=73 } { incr r1 } {
#puts -nonewline $outfile "\t$r1"
# }
# puts -nonewline $outfile "\n"
for { set i 1 } { $i<=$nframes } { incr i } {
set time [expr 0.1 * $i]
puts -nonewline $outfile [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]
# for { set r1 -73 } { $r1<=73 } { incr r1 } {
# set sel1 [atomselect top "protein and noh" frame $i]
# set r2 [expr $r1 * (-1)]
# set sel2 [atomselect top "((segname CHI and resid '$r1') or (segname CHJ and resid '$r2')) and noh" frame $i]
# set contacts [lindex [measure contacts  3.0 $sel1 $sel2] 1]
$h3_1 frame $i
$h3_2 frame $i
$h4_1 frame $i
$h4_2 frame $i
$h2a_1 frame $i
$h2a_2 frame $i
$h2b_1 frame $i
$h2b_2 frame $i
$dnai frame $i
$dnaj frame $i

set rh3_1 [measure rmsd $h3_1_0 $h3_1]
set rh3_2 [measure rmsd $h3_2_0 $h3_2]
set rh4_1 [measure rmsd $h4_1_0 $h4_1]
set rh4_2 [measure rmsd $h4_2_0 $h4_2]
set rh2a_1 [measure rmsd $h2a_1_0 $h2a_1]
set rh2a_2 [measure rmsd $h2a_2_0 $h2a_2]
set rh2b_1 [measure rmsd $h2b_1_0 $h2b_1]
set rh2b_2 [measure rmsd $h2b_2_0 $h2b_2]
set rdnai [measure rmsd $dnai_0 $dnai]
set rdnaj [measure rmsd $dnaj_0 $dnaj]


# set nc [llength $contacts]
# #puts "$contacts"
# #puts [format "Check: %d" "$delta"]
# puts -nonewline $outfile "\t$nc"

# $sel2 delete
# $sel1 delete
# #hbonds -sel1 $sel1 -sel2 $sel2 -dist 3.0 -ang 30 -frames $i:$i -writefile yes -outfile ../analysis_data/dna_hbonds.dat
# #-type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

# }
puts $outfile "\t$rh3_1\t$rh3_2\t$rh4_1\t$rh4_2\t$rh2a_1\t$rh2a_2\t$rh2b_1\t$rh2b_2\t$rdnai\t$rdnaj"
}

 close $outfile 

#hbonds -sel1 $chi -sel2 $chj -dist 3.0 -ang 30 -type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

exit


