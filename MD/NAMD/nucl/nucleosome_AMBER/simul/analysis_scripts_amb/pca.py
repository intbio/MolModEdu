#!/usr/bin/env python

"""
Do PCA analysis of nucleosome trajectory

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv
from pylab import *
from prody import *
#import prody as pr
#import scipy


def main():
	

	
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['ps.useafm'] = True
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['pdf.fonttype'] = 42
	alpha_core='((segment CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segment CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segment CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segment CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA'
	
	
	dcd = DCDFile('../analysis_data/md_nucl.dcd')
	pdb=parsePDB('../analysis_data/only_nucl_aligned.pdb')

	dcd.link(pdb)

	dcd.setAtoms(pdb.select(alpha_core))
	
	edaHF = EDA('Histone fold analysis')
	edaHF.buildCovariance( dcd[2500:] )
	edaHF.calcModes()
	HF_proj=calcProjection(dcd[0:], edaHF[:2]);
	print calcFractVariance(edaHF).round(2)

	writeNMD('../analysis_data/eda_hfolds.nmd', edaHF[:5], pdb.select(alpha_core))

	dcd.setAtoms(pdb.select("name P"))
	edaDNA = EDA('DNA')
	edaDNA.buildCovariance( dcd[2500:] )
	edaDNA.calcModes()
	print calcFractVariance(edaDNA).round(4)
	print calcCovariance(edaDNA).round(4)
	DNA_proj=calcProjection(dcd[0:], edaDNA[:2]);
	showProjection(dcd[0:],edaDNA[:3])
	#print DNA_proj

	figure(figsize={12,5})

	

	subplot(121)
	plot( DNA_proj[:,0], DNA_proj[:,1],'bo', label="DNA P")
	plot(HF_proj[:,0],HF_proj[:,1],'ro', label="Histone fold CA")
	
	xlabel("Principle mode 1, A")
	ylabel("Principle mode 2, A")
	legend(loc="upper right")
	legend(numpoints=1)
	title("Principle components of motion in nucleosome")

	subplot(122)

	plot(np.arange(0, float(len(HF_proj[:,0]))/10,0.1), DNA_proj[:,0],'b-', label="DNA P, mode 1")
	plot(np.arange(0, float(len(HF_proj[:,0]))/10,0.1), DNA_proj[:,1],'g-',linewidth=2, label="DNA P, mode 2")
	plot(np.arange(0, float(len(HF_proj[:,0]))/10,0.1), HF_proj[:,0],'r-', label="Histone fold CA, mode 1")
	plot(np.arange(0, float(len(HF_proj[:,0]))/10,0.1), HF_proj[:,1],'y-',linewidth=2, label="Histone fold CA, mode 2")
	
	#plot( DNA_proj[:,0], DNA_proj[:,1],'bo', label="DNA P")

	xlabel("Time, ns")
	ylabel("Principle mode value, A")
	legend(loc="upper right")
	legend(numpoints=1)
	title("Principle components of motion in nucleosome")

	tight_layout()
	#show()
	savefig('../analysis_data/pca_plot.png',dpi=(200))
	#eda_ensemble
	#for mode in eda[:5]:
		#print calcFractVariance(mode).round(2)

	# mdm2ca_sim1 = dcd[:500]
	# mdm2ca_sim1.superpose()
	# mdm2ca_sim2 = dcd[500:]
	# mdm2ca_sim2.superpose()

 # We project independent trajectories in different color
	# showProjection(mdm2ca_sim1, eda[:3], color='red', marker='.');
	# showProjection(mdm2ca_sim2, eda[:3], color='blue', marker='.');
 # # Now let's mark the beginning of the trajectory with a circle
	# showProjection(mdm2ca_sim1[0], eda[:3], color='red', marker='o', ms=12);
	# showProjection(mdm2ca_sim2[0], eda[:3], color='blue', marker='o', ms=12);
	#  # Now let's mark the end of the trajectory with a square
	# showProjection(mdm2ca_sim1[-1], eda[:3], color='red', marker='s', ms=12); 
	# showProjection(mdm2ca_sim2[-1], eda[:3], color='blue', marker='s', ms=12);

	
	writeNMD('../analysis_data/eda_DNA.nmd', edaDNA[:5], pdb.select("name P"))
	
	#show()
	# for i, frame in enumerate(traj):
	# 	rgyr[i]=calcGyradius(protein)

	# plt.figure()
 
	# plt.title("Gyration radius of alpha helix C-alpha")
	# plt.ylabel("Energy, kcal/mol")        
	# plt.plot(rgyr,'b-',linewidth=2)
	# plt.xlabel('Time, ns')
	# plt.ylabel('Radius of gyrations')
	# plt.savefig("../analysis_data/gyration_ahelix.png",dpi=(600))

	
	#plt.show()
	# ensemble.setAtoms(pdb.calpha)
	# repr(ensemble)
	# rmsd=ensemble.getRMSDs()
	# rmsf=ensemble.getRMSFs()
	# ensemble.superpose()
	# print rmsf[:10]





if __name__ == '__main__':
	main()
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 