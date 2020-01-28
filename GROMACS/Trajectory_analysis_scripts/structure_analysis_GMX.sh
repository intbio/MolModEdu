#! /bin/bash
# set name of the reference pdb and ndx file as well as project (e.g. name of the protein!)
#
dir=$(pwd)
# indicate your gender
sex='sir' # change to madame if you are a lady
# set name of the project
project="EvXyn_bb"
# set a path to the trajectory
trajectory=${dir}/netcdf/${project}/replicas/EvXynTS_apo_340K_1rep.xtc
trajectory_name='EvXynTS_apo_340K_1rep'
# set a reference pdb and index filex (gromacs) to be used with trajecotry
reference=${dir}/reference/${project}
ndx=${dir}/reference/${project}



# output for xmgrace plots
structure=${dir}/output/structure
# output for snapshots
snapshots=${dir}/temp_snap



# indicate structural fragment that are going to be analysed
# E.g. here we are going to analyse two different loops o protein, thus indicating range of the amino-acid residues composed those segments
thumb_loop='res_com of resnr 56 plus res_com of resnr 148' #NpXyn, thumb loop
second_loop='res_com of resnr 56 plus res_com of resnr 107' #NpXyn, 2nd loop

# options for analysis with GROMACS
# note: use ECHO to pipe a number of index group provided in ndx file
# extract and save snapshots from trajectory
echo 0 0 | gmx trjconv -f ${trajectory} -s ${reference} -dt 2000 -fit rot+trans -o ${snapshots}/snapshots.${trajectory_name}.pdb
#measure distance for thumb loop and the second loop
gmx distance -f ${trajectory} -s ${reference} -n ${ndx} -select "${thumb_loop}" -oall ${structure}/Dist.Thumb_${trajectory_name} -oh ${structure}/Hist.Thumb_${trajectory_name} -tu ns -len 2.5 -binw 0.05
gmx distance -f ${trajectory} -s ${reference} -n ${ndx} -select "${second_loop}" -oall ${structure}/Dist.Second_${trajectory_name} -oh ${structure}/Hist.Second_${trajectory_name} -tu ns -len 2.5 -binw 0.05
# calculate SASA for theloops
gmx sasa -f ${trajectory} -s ${reference} -n ${ndx} -surface 'resnr 149 to 155' -o ${structure}/SASA.Thumb_${trajectory_name}
gmx sasa -f ${trajectory} -s ${reference} -n ${ndx} -surface 'resnr 102 to 110' -o ${structure}/SASA.2nd loop_${trajectory_name}
# calculate RMSD during MD trajectory comparing to the 1st snapshot
echo 10 10 | gmx rms -f ${trajectory} -s ${reference} -n ${ndx} -o ${structure}/rmsd_${trajectory_name} -tu ns # -len 2.5 -binw 0.05
# calculate B-factors (RMSF)
echo 0 0 | gmx rmsf -f ${trajectory} -s ${reference} -n ${ndx} -o ${structure}/rmsf_${trajectory_name} -res -fit -oq ${structure}/B-factors_${trajectory_name}

echo "Full analysis for $trajectory_name has been finished! It is time to change a trajectory, ${sex}!"

