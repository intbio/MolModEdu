# the idea is to create a carefully aligned pdb strucutre of initial nucleosome
# at center and oriented along principal axes.
#######################
#Let's center and orient nucleosome, using our fancy orientation and centering algorithm
# We find the geometrical center of DNA helix, and orient wisely.




mol load pdb 1kx5_notails.pdb


mol ssrecalc top

# now we need to adjust ss only to histone fold domains
# and unify helix representation only according
# to the first set of chains (there are slight variations in initial structure)

#seems the easiest way is to do it by hand, this should work for all models,
#since resids are taken with respect to sequence

#set alpha_core [atomselect top "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and backbone"]
set alpha_core [atomselect top "((chain A E and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (chain B F and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (chain C G and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (chain D H and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and backbone"]



set all [atomselect top "all"]
set nucl [atomselect top "protein or nucleic"]



##########let's center nucleosome at origin at first####this is optional
#since we will recenter anyway
proc center_of_mass {selection} {
        # some error checking
        if {[$selection num] <= 0} {
                error "center_of_mass: needs a selection with atoms"
        }
        # set the center of mass to 0
        set com [veczero]
        # set the total mass to 0
        set mass 0
        # [$selection get {x y z}] returns the coordinates {x y z} 
        # [$selection get {mass}] returns the masses
        # so the following says "for each pair of {coordinates} and masses,
  #  do the computation ..."
        foreach coord [$selection get {x y z}] m [$selection get mass] {
           # sum of the masses
           set mass [expr $mass + $m]
           # sum up the product of mass and coordinate
           set com [vecadd $com [vecscale $m $coord]]
        }
        # and scale by the inverse of the number of atoms
        if {$mass == 0} {
                error "center_of_mass: total mass is zero"
        }
        # The "1.0" can't be "1", since otherwise integer division is done
        return [vecscale [expr 1.0/$mass] $com]
}


set vec [center_of_mass $alpha_core]
set vec2 [vecscale $vec -1]
$all moveby $vec2



#########
# ### now we need to orient
 lappend auto_path /Users/alexeyshaytan/soft/vmd_scripts/la1.0
 lappend auto_path /Users/alexeyshaytan/soft/vmd_scripts/orient
 package require Orient
 namespace import Orient::orient


 #orient the symmetry axis with Y
set I [draw principalaxes $alpha_core]
set A [orient $alpha_core [lindex $I 1] {0 1 0}]
$all move $A
set I [draw principalaxes $alpha_core]

# set dna [atomselect top "nucleic and y < 0"]



# #now let's try to get other axis from DNA

#the border lies between 38 and 39 bp on chain I
# we take +- 10 bp each direction

set selDNAF [atomselect top "(chain I and resid 28 to 49) or (chain J and resid '-49' to '-28') "]
set selDNAB [atomselect top "(chain I and resid '-49' to '-28') or (chain J and resid 28 to 49) "]

# set selDNAB [atomselect top "same residue as (nucleic and y < 0 and ( (segname CHI and resid > 0) or (segname CHJ and resid < 0) ))"]

 set vecF [center_of_mass $selDNAF]
 set vecB [center_of_mass $selDNAB]

 set vecZ [vecsub $vecB $vecF]


# ####
proc vmd_draw_arrow {mol start end} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
    graphics $mol cylinder $start $middle radius 1.55
    graphics $mol cone $middle $end radius 1.65
}

 

####

set A [orient $alpha_core $vecZ {0 0 1}]
$all move $A
set I [draw principalaxes $alpha_core]

 set vecF [center_of_mass $selDNAF]
 set vecB [center_of_mass $selDNAB]

 set vecZ [vecsub $vecB $vecF]

set vecY [lindex $I 1]

set vecX [veccross $vecZ $vecY]

set vecZ [vecscale $vecZ 2 ]
set vecY [vecscale $vecY 2 ]
set vecX [vecscale $vecX 2 ]

# vmd_draw_arrow 0 {0 0 0} $vecZ
# vmd_draw_arrow 0 {0 0 0} $vecX
# vmd_draw_arrow 0 {0 0 0} $vecY

vmd_draw_arrow 0 {0 0 0} {0 0 100}
vmd_draw_arrow 0 {0 0 0} {100 0 0}
vmd_draw_arrow 0 {0 0 0} {0 100 0}
#$all writepdb ../analysis_data/nucl_orient_sol.pdb
#$nucl writepdb ../analysis_data/nucl_orient_sol.pdb


#######################
#Getting other axis as the vector betwee DNA gyres

set n [atomselect top "nucleic and ( (chain I and resid 0 ) or (chain J and resid 0))"]
set vec [center_of_mass $n]

for {set i 1} {$i < 39} {incr i} {
    set n1 [atomselect top "nucleic and ( (chain I and resid $i ) or (chain J and resid '-$i'))"]
    set n2 [atomselect top "nucleic and ( (chain I and resid '-$i' ) or (chain J and resid $i))"]
    set vec1 [center_of_mass $n1]
    set vec2 [center_of_mass $n2]
    set vec [vecadd $vec $vec1]
    set vec [vecadd $vec $vec2]
   # puts $vec
}

#we need to divide by the number of vectors used
set vec [vecscale $vec -0.012987012987]


$all moveby $vec

draw color red
draw sphere {0 0 0} radius 2.0

display projection orthographic

mol modstyle 0 0 NewCartoon
mol modselect 0 0 protein or nucleic

$all writepdb 1kx5_notails_aligned.pdb


#exit


