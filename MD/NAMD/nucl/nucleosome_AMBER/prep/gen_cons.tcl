# This file is a VMD/psfgen script to generate topology
# for nucleosome in CHARMM 36 ff
# 
# Will generate a solvated nucleosome with 0.15M ions,
#distance to box 20 A
#

mol new 1kx5nt_amber.prmtop type parm7
mol addfile 1kx5nt_amber.inpcrd type rst7


#Let's get the box
set everyone [atomselect top water]	 
set box [measure minmax $everyone]
echo "DX:"
expr [lindex [lindex $box 1] 0] - [lindex [lindex $box 0] 0]
echo "DY:"
expr [lindex [lindex $box 1] 1] - [lindex [lindex $box 0] 1]
echo "DZ:"
expr [lindex [lindex $box 1] 2] - [lindex [lindex $box 0] 2]
echo "Center:"
measure center $everyone 



#let's mark fixed atoms with 1 beta factor
set all [atomselect top all]
#$all writepdb dna_cons.pdb
#$all set beta 0

set fixed [atomselect top "fragment 0 1 2 3 4 5 6 7 8 9"]
$fixed set beta 1
$all writepdb 1kx5nt_fixed.pdb

# #Now generate file with constraint constants
$all set beta 0
set cons [atomselect top "name CA P"]
$cons set beta 300
$all writepdb 1kx5nt_cons.pdb



###Here we need to generate file with constraints
# for the removal of nucleosome diffusion as a whole
# I suggest to constrain CA of histone folds of H3
# maximum RMSD of H3 CA in a structure aligned using all CA
# in simulations without restraints was 1.45 A
# What about fluctuations maximum rmsf?
# maximum is 4.5 A - but I thinks it's not important
# The number of atoms here is 55*2=110
# We need a potential so that 2 A deviation of all structure would be a kT
# 
# k/2*(2)^2*110=0.6
# k = 0.003

# This is original NAMD selection, AMBER is complicated though
#set cons [atomselect top "(segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) and name CA"]
# for AMBER
#set cons [atomselect top "(fragment 0 4 and (resid 64-36 to 78-36 or resid 86-36 to 114-36 or resid 121-36 to 131-36)) and name CA"]
# chain E should have 395+ values for chain A

$all set beta 0
set cons [atomselect top "((resid 28 to 42 or resid 50 to 78 or resid 85 to 95) or (resid 423 to 437 or resid 445 to 473 or resid 480 to 490)) and name CA"]
# We have here 55*2=110 atoms
$cons set beta 1
$all writepdb 1kx5nt_h3_cons.pdb

# exec bash -c "rm t_*"
exec cp 1kx5nt_amber.inpcrd ../simul/input/
exec cp 1kx5nt_amber.prmtop ../simul/input/
exec cp 1kx5nt_fixed.pdb ../simul/input/
exec cp 1kx5nt_cons.pdb ../simul/input/
exec cp 1kx5nt_h3_cons.pdb ../simul/input/
exit


# {1.412140965461731 1.6916090250015259 1.6624139547348022} {149.9282989501953 148.25173950195313 107.44857025146484}
# DX:
# 148.51615798473358
# DY:
# 146.5601304769516
# DZ:
# 105.78615629673004
# Center:
# 75.98223114013672 74.9022445678711 54.71516799926758