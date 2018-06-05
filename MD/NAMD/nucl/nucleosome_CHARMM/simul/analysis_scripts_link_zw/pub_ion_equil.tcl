#ion equilibraion
mol load psf ../input/1kx5_ready.psf
mol addfile ../analysis_data/system_aligned.pdb waitfor all
mol addfile ../output/min_eq.dcd waitfor all

# mol load psf ../input/1kx5_ready.psf
#mol addfile ../analysis_data/nucl_orient_sol.pdb waitfor all
#mol ssrecalc top
# mol addfile ../analysis_data/md_nucl_solv.dcd first 0 waitfor all


package require pbctools
set nframes [expr  [molinfo top get numframes] - 1 ]
#set ions [atomselect top "resname SOD"]

### The following must be simplified in a case where nucleosome is fixed!!!
### we can delete water, and set center to origin!!!
pbc wrap -all -center com -compound res  -centersel "protein" -sel "resname SOD or resname CLA or resname TIP3"

set outfile [open ../analysis_data/ion_equil.dat w]	 
set nf [molinfo top get numframes]	 

#set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"



set nucl [atomselect top "(protein or nucleic) and noh"]	

set dna_p [atomselect top "name P"]	 
set sod [atomselect top "resname SOD"]
set cla [atomselect top "resname CLA"]



puts $outfile "Number of ions within 1 nm from nucleosomal DNA with time"
puts $outfile "Calculated in VMD using $nf frames"
puts $outfile "Time, steps"
puts $outfile "Number of ions"
puts $outfile "Time\tSodium\tChloride"

set s1 [atomselect top "resname SOD"]
set s2 [atomselect top "resname CLA"]
set nuclcharge [expr [$s1 num] - [$s2 num]]


for {set i 0} { $i < $nf} {incr i} {

#$cla frame $i
#puts "Processing frame $i"
set sel [atomselect top "resname SOD and within 10 of (nucleic)" frame $i]
set pc [$sel num]
$sel delete
set sel [atomselect top "resname CLA and within 10 of (nucleic)" frame $i]
set nc [$sel num]
$sel delete

puts $outfile "$i\t$pc\t$nc"
puts "$i\t$pc\t$nc"
}

#nucleosome charge is this system is -202

# set sod_rdf [measure gofr $dna_p $sod rmax 15 usepbc 1 first 500 last 999]
# set cla_rdf [measure gofr $dna_p $cla rmax 15 usepbc 1 first 500 last 999]

# set sod_rdf_10 [measure gofr $dna_p $sod rmax 15 usepbc 1 first 0 last 100]
# set cla_rdf_10 [measure gofr $dna_p $cla rmax 15 usepbc 1 first 0 last 100]

# set na_sod_rdf [measure gofr $sod $cla rmax 15 usepbc 1 first 0 last 100]
# # rmsd calculation loop
# #puts $sod_rdf

# foreach r [lindex $sod_rdf 0] g_sod [lindex $sod_rdf 1] g_cla [lindex $cla_rdf 1] g_sod_10 [lindex $sod_rdf_10 1] g_cla_10 [lindex $cla_rdf_10 1]  g_na_sod [lindex $na_sod_rdf 1] {	 

# puts $outfile "$r\t$g_sod\t$g_cla\t$g_sod_10\t$g_cla_10\t$g_na_sod"
 
# }	 
close $outfile 

exec python2.7 myplot.py ../analysis_data/ion_equil.dat 

exit


