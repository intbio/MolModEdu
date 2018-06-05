# we calculate ion rdfs


mol load psf ../input/1kx5_ready.psf
#mol addfile ../analysis_data/nucl_orient_sol.pdb waitfor all
#mol ssrecalc top
mol addfile ../output/md.dcd last 10000 waitfor all

set outfile [open ../analysis_data/ion_rdf.dat w]	 
#set nf [molinfo top get numframes]	 

#set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"


set dna_p [atomselect top "name P"]	 
set sod [atomselect top "resname SOD"]
set cla [atomselect top "resname CLA"]

set sod_rdf [measure gofr $dna_p $sod rmax 15 usepbc 1 first 5000 last 9999]
set cla_rdf [measure gofr $dna_p $cla rmax 15 usepbc 1 first 5000 last 9999]

set sod_rdf_10 [measure gofr $dna_p $sod rmax 15 usepbc 1 first 0 last 100]
set cla_rdf_10 [measure gofr $dna_p $cla rmax 15 usepbc 1 first 0 last 100]

set cla_sod_rdf [measure gofr $sod $cla rmax 15 usepbc 1 first 5000 last 9999]
# rmsd calculation loop
#puts $sod_rdf

puts $outfile "RDF of ions in nucleosome simulation"
puts $outfile "Calculated in VMD"
puts $outfile "Distance, A"
puts $outfile "RDF"
puts $outfile "Distance\tP-SOD (last 500ns)\tP-CLA (last 500ns)\tP-SOD (first 10ns)\tP-CLA (first 10)\tCLA-SOD (last 500ns)"
foreach r [lindex $sod_rdf 0] g_sod [lindex $sod_rdf 1] g_cla [lindex $cla_rdf 1] g_sod_10 [lindex $sod_rdf_10 1] g_cla_10 [lindex $cla_rdf_10 1]  g_cla_sod [lindex $cla_sod_rdf 1] {	 

puts $outfile "$r\t$g_sod\t$g_cla\t$g_sod_10\t$g_cla_10\t$g_cla_sod"
 
}	 
close $outfile 

exec python2.7 myplot.py ../analysis_data/ion_rdf.dat

exit


