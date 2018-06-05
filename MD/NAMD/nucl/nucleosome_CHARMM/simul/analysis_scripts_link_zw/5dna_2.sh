#!/bin/bash

#This is a reboot on DNA analysis.
#We aim to make python libraries, that would interface to 3DNA and Curves+
#directly

# We need here to characterize the dna conformation
# The main programs to do this are 3DNA and Curves+
# 3DNA uses a local helical parameters, while Curves+
# uses global, it allows to calculate groove width, depth,
# as well as curvature


#We start by studying the geometry of DNA in terms of a polygon approximation.

vmdt -e dna_bl_geom.tcl
/usr/bin/R --vanilla --slave < dna_bl_analyze.r



#Let's plot the evolution of RMSD of individual basepairs
vmdpy -e dna_rmsd.tcl
/usr/bin/R --vanilla --slave < dna_rmsd.r


#Let's run VMD script that will collect all information on DNA
#conformation and put it into data frame

vmdtpy  -python -e dna2_get_df_cryst.vmdpy
vmdtpy  -python -e  dna2_get_df_md.vmdpy


# Let's start with DNA base-pair parameters from x3DNA
# Here we compute parameter distribution and profiles
# check that all BP are recognized!!! see comment in dna_params.tcl
vmdt -e dna_params_3dna.tcl # this makes data frames and panel for analysis

# let's stick to R and do statistical analysis, no dynamics
/usr/bin/R --vanilla --slave < dna_par_3dna_analyze.r

/usr/bin/R --vanilla --slave < dna_par_3dna_analyze_hist.r

#TODO: evolution of parameters with dynamics, heat plots (?)

# Now let's go to Curves+
# Curves+ get major - minor goove width,
# B1-B2 transition
vmd -e dna_params_cur.tcl
/usr/bin/R --vanilla --slave < dna_bb_analyze.r
/usr/bin/R --vanilla --slave < dna_gr_analyze.r

#TODO: check also Sumr!




