# Covariance analysis of the MD trajectory using prody package
# http://prody.csb.pitt.edu/tutorials/trajectory_analysis/eda.html
# Warning! both trajectory.dcd and reference.pdb structure should be of the same number of atoms asuming perfect matching between both

from prody import *



structure=  parsePDB('reference.pdb')
#Unlock below string if you'd like select specified atom subset from the reference structure
#structure = structure.select("resnum 5 to 147 and calpha")
structure = structure.select("calpha") # select only c-alpha sub-sets from the reference structure

# load trajectory in DCD format
trajectory = parseDCD('./trajectory.dcd')
# fit reference atoms to the trajectory
trajectory.setAtoms(structure.calpha)
# Snapshot superimposition
trajectory.superpose()

# Make ensemble from the MD snapshots
md_ensemble = EDA('Name of the trajectory')
# Build covariance matrix
md_ensemble.buildCovariance( trajectory )
# Calc 100 modes
md_ensemble.calcModes(n_modes=100)

# Include all snapshots in the trajectory
trajectory = trajectory[:-1]
# Make superimposition of that snapshots
trajectory.superpose()



# VISUALISATION! create NMD files which could be visualized in the Normal Mode Wizard plugin of VMD
writeNMD('MD.nmd', md_ensemble[:30], structure.select('calpha'))

# Plotting
import matplotlib.pyplot as plt
# Calculate fraction of variance of individual modes
plt.figure(figsize=(9,8))
showFractVars(md_ensemble) 
showCumulFractVars(md_ensemble) 
plt.savefig( "./frac_of_variance.png", dpi=300 )

# project trajectory onto first 2 lowest frequency modes ( might be changed to any modes)
plt.figure(figsize=(9,8))
showProjection(trajectory, md_ensemble[:2], color='green')
showProjection(trajectory[0], md_ensemble[:2], color='green', marker='o', ms=12)
showProjection(trajectory[-1], md_ensemble[:2], color='green', marker='s', ms=12)
plt.savefig( "./MD_projection.png", dpi=300 )

#calculate fluctuations for each C-alpha atom along 20 lowest modes
plt.figure(figsize=(9,8))
showNormedSqFlucts(md_ensemble[:20])
plt.legend(prop={'size': 10})
plt.savefig( "./Fluctuations_along_20modes.png", dpi=300 )

#calculate fluctuations for each C-alpha atom along 1 lowest mode
plt.figure(figsize=(9,8))
showNormedSqFlucts(md_ensemble[0])
plt.legend(prop={'size': 10})
plt.savefig( "./Fluctuations_along_1st_mode.png", dpi=300 )

# calculate cross-correlation map along 20 lowest modes. It could be usefull for the analysis of the oligomer proteins ( fluctuations in the distal monomers).
plt.figure(figsize=(9,8))
showCrossCorr(md_ensemble[:20])
plt.savefig( "./cross-correlations.png", dpi=300 )
