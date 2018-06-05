#!/bin/bash

# Global conformational dynamics analysis
#In this section we are interested in describing the overall conformation of the system, how it changes, how it differs from the structure observed in the crystal.
#In this section we will use several approaches to reveal global conformation dynamics: (1) RMSD for the whole structure and its parts (histones, DNA, backbone, side chains) relative to 
#the crystal structure and relative to the average structure, (2) analyze the average structure and see how it differs from crystal, (3) do cluster analysis of the trajectory, (4) analyze the root-mean-square fluctuations of atoms and reveal more and less dynamically flexible elements, (5) check the symmetry of the two halves of nucleosome. 
#


#1. RMSD
# nucleosome is aligned in trj after preparation.sh by minimizing RMSD between
# CA atoms of of all histone fold helices
#
# note: vmd does not automatically align structures when doing RMSD
# so RMSD values for the parts actually would include the component due to 
# movment of the part as a whole in nucleosome structure if we did not align them during
# preparation of trajectory
#
# What components are interesting to calculate as evolution graphs
# 1) RMSD of all, DNA (all and backbone), helices (all and CA), histone folds (CA)
# important note: helices are determined by STRIDE using crystal structure.
# there is slight +-1 amino acid difference in helix assignment in symmetric histones
# but this should not influence the results
./rmsd_min_eq.sh #RMSD relative to crystal structure for equilibration - just technical thing to check

./rmsd_md.sh


# Here is an attempt to check for RMSD values for different parts of the system, a-la Biswas
# as final averaged values (let say, averaged during last 50 ns of run,
# or after whatever time we define as relaxation period) we can compare:
# (sidechains, backbone, CA, bases, phosphates, sugars) * (all, dna, alpha helix,
# histone fold, each histone pair separately x (all, alpha helix, histone fold),)
#
# We just effectively here split the RMSD for all atoms of nucleosome into parts
#

./parts_rmsd.sh

#
# 3) finally we can analyze symmetry 
# let's do RMSD of halves and DNA strands self RMSD
# with time.
#
#
./symmetry_rmsd.sh

####Average structure analysis
#####
#Let's get the average structure - see the exact time frame in files
vmdt -e av_struct.tcl #generate an average strucure for the whole run.

#Now we need to see it's difference from crystal.
vmdt -e dif_av_cryst.tcl
vmd -e visual_dif_av_cryst.tcl
vmdpy -e visual_dif_av_cryst_overl.tcl
python2.7 myplot_rmsf.py ../analysis_data/dif_av_cryst_chains.dat  --legend "upper center"  --max_y 3 3 3 3 12
python2.7 myplot_rmsf.py ../analysis_data/dif_av_cryst_chains_sch.dat





#Let't do cluster analysis
vmdt -e cluster_md.tcl
python2.7 myplot_cluster.py ../analysis_data/rmsd_md_av.dat ../analysis_data/cluster_md_dnabb.dat --x_r 0 1000
python2.7 myplot_cluster.py ../analysis_data/rmsd_md_av.dat ../analysis_data/cluster_md_hfolds.dat --x_r 0 1000



#
# Now we porceed with RMSF
# I think RMSF is related to B-factor.
# But RMSF has more direct meaning.
# So we will stick to it.
# What we need here is to get a pdb-with b-factor set to RMSF
# we need a script to open interactive visualization
# we need a script to make movie with a bar
#
# we will take the last 50 ns of trajectory - change it in scripts if needed!

./get_rmsf.sh # this will calculate and plot RMSF
#vmd -e visual_rmsf.tcl # Adjust scale minmax values here to get best view

#view_rmsf.sh - is available to view RMSF

#you can also try to use prody_analysis.py to extract rydius of giration,
# and as a stub for other prody analysis functions.





#
#
#
# PCA analysis
#
# We need to use Prody to visualize main collective movement modes 
# of protein core and DNA
#
#
./pca.sh
vmd -e visual_eda_dna.tcl
vmd -e visual_eda_hfolds.tcl






