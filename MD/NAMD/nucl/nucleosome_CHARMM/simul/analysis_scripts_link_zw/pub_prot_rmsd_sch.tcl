#we need to make a heat plot of RMSD evolution of certain bases
#
# 
#
#


mol load psf ../analysis_data/only_nucl_init.psf
mol addfile ../analysis_data/only_nucl_aligned.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl.dcd step 10 first 1 waitfor all

source def_struct_param.tcl


set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile [open ../analysis_data/pub_prot_rmsd_sch.csv w]	

puts "Frames"
puts $nframes

# set frame0_alpha_all_ca [atomselect top $alpha_all_ca frame 0]	 
# set sel_alpha_all_ca [atomselect top $alpha_all_ca]	
# set atomsid [$frame0_alpha_all_ca get index]

set frame0_all_ca [atomselect top $all_ca frame 0]	 
set sel_all_ca [atomselect top $all_ca]	
set atomsid [$frame0_all_ca get index]

# set chi [atomselect top "segname CHI"]	 
# set chj [atomselect top "segname CHJ"]	 
# rmsd calculation loop


# puts $outfile "Prot RMSD profile"
# puts $outfile "Calculated using VMD"
# puts $outfile "Time, ps"
# puts $outfile "RMSD, A"
puts -nonewline $outfile "Time\tChain\tResid\tRMSD"

puts -nonewline $outfile "\n"

for { set i 1 } { $i<=$nframes } { incr i } {

set time [expr 1.0 * $i]
puts  [format "Time %.3f" "$time"]

foreach a $atomsid {

set sel_ids [atomselect top "index $a" frame 0]

set resname [$sel_ids get resname]

if { $resname eq "GLY" } {
set rmsd 0 
} else {
set sel2 [atomselect top "(same residue as (index $a)) and not name C O N CA and noh" frame $i]
set sel0 [atomselect top "(same residue as (index $a)) and not name C O N CA and noh" frame 0]

set rmsd [measure rmsd $sel2 $sel0]


$sel2 delete
$sel0 delete

}

set chain [$sel_ids get segname]

set resid [$sel_ids get resid]

puts -nonewline $outfile [format "%.3f" "$time"]
puts $outfile "\t$chain\t$resid\t$rmsd"


$sel_ids delete

}
puts "Kuku"
}


close $outfile 


exit






