#!/usr/bin/env python

"""
Do simple prody analysis of nucleosome trajectory
1. Gyration radius

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv
from prody import *
#import prody as pr
#import scipy


def main():
	"Do prody"

	
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['ps.useafm'] = True
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['pdf.fonttype'] = 42
	traj = Trajectory('../analysis_data/md_nucl.dcd')
	pdb=parsePDB('../analysis_data/only_nucl_aligned.pdb')
	traj.link(pdb)
	traj.setCoords(pdb)
	rgyr=np.zeros(traj.numFrames())
	alpha_core='((segment CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segment CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segment CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segment CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA'
	
	

	protein=pdb.select(alpha_core)
	


	for i, frame in enumerate(traj):
		rgyr[i]=calcGyradius(protein)

	plt.figure()
 
	plt.title("Gyration radius of alpha helix C-alpha")
	plt.ylabel("Energy, kcal/mol")        
	plt.plot(rgyr,'b-',linewidth=2)
	plt.xlabel('Time, ns')
	plt.ylabel('Radius of gyration, A')
	plt.savefig("../analysis_data/gyration_ahelix.png",dpi=(600))

	
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