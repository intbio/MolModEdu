#!/usr/bin/env python

"""
Plot var-covar for DNA backbone

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv
from pylab import *
from prody import *
import pickle
#import prody as pr
#import scipy
EVNUM=40

def main():
	

	
	matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':14})
	plt.rcParams['ps.useafm'] = True
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['pdf.fonttype'] = 82
	# matplotlib.rc('font', **font)
	alpha_core='((segment CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segment CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segment CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segment CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA C O N'
	#											15				  29					11 										11 					28 				11 											11 					29 				10 											12 				29 					11



		# sp.plot(np.arange(0, float(edaHF.numAtoms()),1), calcSqFlucts(edaHF[i]), label="%d"%(i+1),linewidth=2)
		# sp.axvspan(HFind['H3A'][0],HFind['H3A'][-1] , color='blue', alpha=0.2)
		# sp.axvline(x=HFind['H3A'][1],color='k')
		# sp.axvline(x=HFind['H3A'][2],color='k')
		# sp.axvspan(HFind['H4B'][0],HFind['H4B'][-1] , color='green', alpha=0.2)
		# sp.axvline(x=HFind['H4B'][1],color='k')
		# sp.axvline(x=HFind['H4B'][2],color='k')
		# sp.axvspan(HFind['H2AC'][0],HFind['H2AC'][-1] , color='yellow', alpha=0.2)
		# sp.axvline(x=HFind['H2AC'][1],color='k')
		# sp.axvline(x=HFind['H2AC'][2],color='k')
		# sp.axvspan(HFind['H2BD'][0],HFind['H2BD'][-1] , color='red', alpha=0.2)
		# sp.axvline(x=HFind['H2BD'][1],color='k')
		# sp.axvline(x=HFind['H2BD'][2],color='k')

		# sp.axvspan(HFind['H3E'][0],HFind['H3E'][-1] , color='blue', alpha=0.2)
		# sp.axvline(x=HFind['H3E'][1],color='k')
		# sp.axvline(x=HFind['H3E'][2],color='k')
		# sp.axvspan(HFind['H4F'][0],HFind['H4F'][-1] , color='green', alpha=0.2)
		# sp.axvline(x=HFind['H4F'][1],color='k')
		# sp.axvline(x=HFind['H4F'][2],color='k')
		# sp.axvspan(HFind['H2AG'][0],HFind['H2AG'][-1] , color='yellow', alpha=0.2)
		# sp.axvline(x=HFind['H2AG'][1],color='k')
		# sp.axvline(x=HFind['H2AG'][2],color='k')
		# sp.axvspan(HFind['H2BH'][0],HFind['H2BH'][-1] , color='red', alpha=0.2)
		# sp.axvline(x=HFind['H2BH'][1],color='k')
		# sp.axvline(x=HFind['H2BH'][2],color='k')


#Now let's do covariance
	#let's dump it
	#pickle.dump(calcCrossCorr(edaHF), open( "../analysis_data/hfolds_var_covar.p", "wb" ) )
	cross_corr = pickle.load( open( "../analysis_data/dnabb_cross_corr.p", "rb" ) )
	# arange = np.arange(1656)
	# cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
	# cross_correlations[arange[0]+1:,arange[0]+1:] = cross_corr
	cross_corr=cross_corr[0:(cross_corr.shape[0]/2),:]
	print cross_corr.shape
	fig4=plt.figure(figsize={15,7})
	sp=fig4.add_subplot(111)
	cax=sp.imshow(cross_corr,interpolation='none')
	fig4.colorbar(cax,shrink=0.75)
	#sp.set_axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
	sp.set_title('Cross-correlation')
	sp.set_xlabel('Super helix location (SHL)')
	sp.set_ylabel('Super helix location (SHL)')
	sp.set_title("Cross-correlation matrix for DNA backbone.")

	# sp.set_xlim(0,100)
#SHL location#0 0.5 1.0 and so on
	SHL_loc=np.array([2,12,22,32,42,52,62,72,82,92,102,112,122,132,142])
	SHL_loc_ticks=np.array([-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7])

	SHL_loc_ticksx=np.concatenate((SHL_loc_ticks,SHL_loc_ticks))
	SHL_locx=np.concatenate((SHL_loc,SHL_loc+146))

	sp.set_xticks(SHL_locx)
	sp.set_xticklabels(SHL_loc_ticksx)

	sp.set_yticks(SHL_loc)
	sp.set_yticklabels(SHL_loc_ticks)
	sp.tick_params( labelright=True)

	sp.axvline(x=145.5,color='k')
	# sp.axvline(x=146,color='k')

	sp.annotate('Chain I', xy=(.2, .06), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top',fontsize=20)
	sp.annotate('Chain J', xy=(.59, .06), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top',fontsize=20)
	sp.annotate('Chain I', xy=(.817, .55), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top',fontsize=20, rotation=270)

	fig4.tight_layout()
	fig4.savefig("../analysis_data/pca_dnabb_cross_corr.png",dpi=(200))

	plt.show()
	# savefig('../analysis_data/pca_plot.png',dpi=(200))

	
	# writeNMD('../analysis_data/eda_DNA.nmd', edaDNA[:5], pdb.select("name P"))
	

	# print "Mode covariance"
	# print calcCovariance(edaHF[0:100])
	# print calcCovariance(edaHF[0:2]).shape
	# print "Orig covar"
	# print edaHF.getCovariance() # returns covariance matrix
	# print edaHF.getCovariance().shape
	# edaHF2 = EDA('Histone fold analysis')
	# #We want to build covariance - since this is dcd object, the frames will be superimposed on crystal using 
	# edaHF2.buildCovariance( dcd[250:],aligned=True )
	# edaHF2.calcModes(n_modes=EVNUM)
	# print "covar overlap"
	# print calcCovOverlap(edaHF,edaHF2)
	# print printOverlapTable(edaHF,edaHF2)


if __name__ == '__main__':
	main()
	



# 