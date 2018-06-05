# Proceedure to generate BFACTORs from a trajectory within VMD
# VERSION 1.0
#
# All selected frames are superimposed on the reference frame (frame 0).
# The coords of the selected atoms in frame 0 are 
# set to the average atomic positions prior to writing
# PDB file. The B column of the PDB file is set to the calculated
# isotropic B-factor.
# 
# By Aaron Oakley 
# Research School of Chemistry
# Australian National University
# Email: oakley@rsc.anu.edu.au
# 
# Usage:
#     sel = atoms used for calculation e.g. "segname A and noh"
#   start = number of first frame for calculation
#     end = number of last frame for calculation
# outfile = pdb file containing average x,y,z coords

proc bfactor {sel start end outfile scaling} {

set reference [atomselect top  "$sel" frame 0]
set compare [atomselect top  "$sel"]
set all [atomselect top all]
set atindex [[atomselect top "$sel"] get index] 
set pi 3.14159265
set numframes [expr $end - $start + 1]
# NB: Frames 10 to 20 has 11 frames, but 20 - 10 = 10!

# The following variables will contain the mean coords x,y,z,
# B-factor:

foreach r $atindex {
        set px($r) 0
        set py($r) 0
        set pz($r) 0
        set U($r) 0
        set B($r) 0
}

# loop over frames in the trajectory
# to calculate sum of x,y,z positions for each atom.

puts "Aligning selected frames to frame 0."
puts "Calculating mean positions of selected atoms $sel ..."

for {set frame $start} {$frame <= $end} {incr frame} {
	puts -nonewline "\rFrame: $frame"
        flush stdout
	# get the correct frame
	$compare frame $frame
	# compute the transformation
	set trans_mat [measure fit $compare $reference]
	# do the alignment
        set all_current [atomselect top all frame $frame]
	$all_current  move $trans_mat
        $all_current delete
	# loop through all atoms
        set getframe [atomselect top "index $atindex" frame $frame]
        set veclist [$getframe get {x y z}]
	foreach r $atindex vec $veclist {
                 set ix [lindex $vec 0]
                 set iy [lindex $vec 1]
                 set iz [lindex $vec 2] 
		set px($r) [expr $px($r) +  $ix ]
		set py($r) [expr $py($r) +  $iy ]
		set pz($r) [expr $pz($r) +  $iz ]
	}
        $getframe delete

}

puts ""
puts "...Done!"

# Now divide by $numframes to give the average coordinate for each atom

foreach r $atindex {
        set px($r) [expr $px($r) / $numframes]
        set py($r) [expr $py($r) / $numframes]
        set pz($r) [expr $pz($r) / $numframes]
}


# Now calculate bfactors for each atom:

puts "Calculating bfactors..."

set getframe [atomselect top "index $atindex"]
for {set frame $start} {$frame <= $end} {incr frame} {
        #loop through all atoms
	puts -nonewline "\rFrame: $frame"
        flush stdout
        $getframe frame $frame
        set veclist [$getframe get {x y z}]
        foreach r $atindex vec $veclist {
                set ix [lindex $vec 0]
                set iy [lindex $vec 1]
                set iz [lindex $vec 2]
                set dx [expr $ix - $px($r)]
                set dy [expr $iy - $py($r)]
                set dz [expr $iz - $pz($r)]
                set Uxx [expr $dx * $dx]
                set Uyy [expr $dy * $dy]
		set Uzz [expr $dz * $dz]
                set U($r) [expr $U($r) + $Uxx + $Uyy + $Uzz]
        }
}
$getframe delete

# Divide by $numframes to give average:

foreach r $atindex {
	set U($r) [expr $U($r) / $numframes]
}

# Conversion for B = 8/3 Pi**2 <Uxx + Uyy + Uzz>.


foreach r $atindex {
        set B($r)   [expr 8.0 / 3 * $pi * $pi * $U($r)]
}

puts ""
puts "...Done!"


# Set coordinates in frame 0 to "average" positions.
# Set beta to calculated B-factor.

foreach r $atindex {
        set getindex [atomselect top "index $r" frame 0]
        $getindex set x  $px($r)
        $getindex set y  $py($r)
        $getindex set z  $pz($r)
        $getindex set beta [expr $B($r) / $scaling]
        $getindex delete
}

# Write pdb file with ANISOU records:

puts "Writing pdb file $outfile ..."

$reference writepdb $outfile

puts "...Done!"
}
