#! /bin/bash
# This script performes Principal Component Analysis of MD trajectory using GROMACs utilities
# it allows to calculate free energy surface by means of projection of the trajectory onto selected number of principal modes
# final chart gives ideas about relative probabilities of the transitions between various conformational states for examined structural domain 
# NB! it is possible to define the structural domain in ndx file created by gromacs
#
#
# set name of the reference pdb and ndx file as well as project (e.g. name of the protein!)
#
dir=../
project="EvXyn_bb"
# set a path to the trajectory
trajectory=${dir}/netcdf/${project}/replicas/EvXynTS_apo_340K_3rep.xtc
# all trajectories merged together in the ensemble
trajectory_name='EvXynTS_apo_340K_1rep'
# set a reference pdb and index filex (for gromacs) to be used with trajecotry
reference=${dir}/reference/${project}
ndx=${dir}/reference/${project}


# output for xmgrace plots
structure=${dir}/output/structure
# output for snapshots
snapshots=${dir}/temp_snap





# step by step PCA analysis with GROMACS

# 0 - make a merged trajectory from all trajectories to create shared subspace for projections
# we are going to project each MD trajectory onto the shared collective plane calculated for all md trajectories together
gmx trjcat  -f ${dir}/netcdf/${project}/replicas/*.xtc -o ./merged_${project}.xtc -cat


#1- calculate covariance matrix for the MD trajectory
echo 0 0 | gmx covar -s ${reference} -n ${ndx} -f ./merged_${project}.xtc -o eigenval.xvg -av averaged_structure.pdb -v eigenvect.trr

#2- calculate PC 1 et PC2 to obtain projections of the md trajectory onto these modes
echo 0 0 | gmx anaeig -v eigenvect.trr -s ${reference} -f ${trajectory} -n ${ndx} -first 1 -last 2 -2d 2d_projection_pc1_2.xvg -filt Pc1_Pc2_motions.pdb

#3- make 2D Free "energy landscape" for this 2D projection

gmx sham -f 2d_projection_pc1_2.xvg -ls gibbs.xpm -notime -ngrid 80 -nlevels 80
gmx xpm2ps -f gibbs.xpm -o pc1_2_landscape.eps -rainbow red -nice 90

echo "PCA analysis is finished!"


