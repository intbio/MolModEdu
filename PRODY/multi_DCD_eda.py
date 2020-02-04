# Covariance analysis (PCA, EDA) of the MD trajectory using Prody python package
# http://prody.csb.pitt.edu/tutorials/nmwiz_tutorial/pca.html

from prody import *
import os
# path to all trajectoris in DCD that will be combined together
path="/mydir/dcd"

# reference pdf
structure =  parsePDB('/mydir/NpXyn_bb.pdb')
# define region of interest to be analysed
sele='resnum 147 to 157'


# 1 - make a shared trajectory
dirs=os.listdir(path)
traj_shared = Trajectory(os.path.join(path, dirs[0]))
traj_shared.setAtoms(structure.select(sele))
for i in range(1, len(dirs)):
    traj_shared.addFile(os.path.join(path, dirs[i]))

writeDCD('merged_npxyn.dcd',traj_shared, step=1,start=1,stop=-1)
writePDB('selection.pdb',structure.select(sele))

# 2 - calculate modes from the shared trajectory
trajectory_big = parseDCD('merged_npxyn.dcd')
trajectory_big.setAtoms(structure.select(sele))
trajectory_big.superpose()
eda = EDA('merged_NEW_hel.dcd')
eda.buildCovariance( trajectory_big )
eda.calcModes()

writeDCD('merged_npxyn_super.dcd',trajectory_big)


# 3 - make an array of all the trajectories
tr= []
modes=[]
for traj in dirs:
 trajectory = parseDCD(os.path.join(path, traj),step=1)
 trajectory.setAtoms(structure.select(sele))
 trajectory=trajectory[:-1]
 trajectory.superpose()
 tr.append(trajectory)
 eda_tr = EDA(os.path.join(path, traj))
 eda_tr.buildCovariance( trajectory )
 eda_tr.calcModes()
 modes.append(eda_tr)




# 4 - plot the projection

from pylab import *
import matplotlib.pyplot as plt
plt.close('all')
plt.figure(figsize=(8,7))
axis([-1.5, 6.0, -2.0, 2.0]);
#legend();
#showProjection(tr[3], eda[0,1], color='red')
#showProjection(tr[3][0], eda[0,1], color='red',marker='o', ms=15)
#showProjection(tr[3][-1], eda[0,1], color='red',marker='s', ms=12)
showProjection(tr[1], eda[0,1], color='orange')
showProjection(tr[1][0], eda[0,1], color='orange',marker='o', ms=15)
showProjection(tr[1][-1], eda[0,1], color='orange',marker='s', ms=12)
#showProjection(tr[2], eda[0,1], color='purple')
#showProjection(tr[2][0], eda[0,1], color='purple',marker='o', ms=15)
#showProjection(tr[2][0], eda[0,1], color='purple',marker='s', ms=12)
showProjection(tr[0], eda[0,1], color='cyan')
showProjection(tr[0][0], eda[0,1], color='cyan',marker='o', ms=15)
showProjection(tr[0][-1], eda[0,1], color='cyan',marker='s', ms=12)
plt.savefig('all_PC2-3_MD_NEW.png',dpi=100)
