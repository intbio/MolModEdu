#we need to make a DNA protein contact map
#
# should be a 
#
#

set main_flag [lindex $argv 0]
set first [lindex $argv 1]
set last [lindex $argv 2]
set file [lindex $argv 3]



mol load psf ../analysis_data/only_nucl_init.psf
mol addfile ../analysis_data/only_nucl_aligned.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd step 10 first 900 waitfor all


set hlist [list "segname CHA" "segname CHE" "segname CHB" "segname CHF" "segname CHC" "segname CHG" "segname CHD" "segname CHH"]
set hname [list "H3_A" "H3_E" "H4_B" "H4_F" "H2A_C" "H2A_G" "H2B_D" "H2B_H"]




set nframes [expr  [molinfo top get numframes]  ]

if {$main_flag == 1} {
set outfile [open $file w]	



puts $outfile "DNA protein contact map"
puts $outfile "Calculated using VMD"
puts $outfile "Protein"
puts $outfile "DNA"



# Let's output two headerlines

foreach sq $hlist sn $hname {
set sel [atomselect top "name CA and $sq"]
puts -nonewline $outfile "[$sel get resid] "
}

puts -nonewline $outfile "\n"

foreach sq $hlist sn $hname {
set sel [atomselect top "name CA and $sq"]

puts -nonewline $outfile "[$sel get resname] "

}
puts -nonewline $outfile "\n"

} else {
set outfile [open $file w]
}

for { set rd $first } { $rd<=$last } { incr rd } {
	set rd2 [expr $rd * (-1)]
	puts "Starting $rd"
	puts -nonewline $outfile "$rd"
	foreach sq $hlist {
		set selp [atomselect top "name CA and $sq"]
		set ids [$selp get resid]
		$selp delete
		foreach rp $ids {
			set n 0.
			for { set i 0 } { $i< $nframes } { incr i } {
				set sel [atomselect top "($sq and resid $rp and noh) and ( within 4 of ( ((segname CHI and resid '$rd') or (segname CHJ and resid '$rd2')) and noh ) )" frame $i]
				set nc [$sel num]
				set n [expr $n + $nc.]
				$sel delete
				}
			set n [expr $n / $nframes.]
			puts -nonewline $outfile "\t$n"
			}
		}
	puts -nonewline $outfile "\n"
	}
close $outfile 


exit


