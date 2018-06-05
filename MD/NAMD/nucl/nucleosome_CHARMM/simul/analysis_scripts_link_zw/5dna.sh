#!/bin/bash


# We need here to characterize the dna conformation
# The main programs to do this are 3DNA and Curves+
# 3DNA uses a local helical parameters, while Curves+
# uses global, it allows to calculate groove width, depth,
# as well as curvature



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

#Let's plot the evolution of RMSD of individual basepairs
vmdpy -e dna_rmsd.tcl

# TODO: characterize the path of DNA globaly, and how DNA conformation changes with time
# I guess the best is to trace the change in COM of basepair positions
# Curves+ can compare structures!!!


