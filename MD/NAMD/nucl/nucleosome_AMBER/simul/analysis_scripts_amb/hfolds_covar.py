#!/usr/bin/env python

"""
Plot var-covar for histone folds

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

	#Let's make correspondence table of histone borders
	HFind=dict()
	HFind['H3A']=np.array([0,15,15+28,15+28+11])
	HFind['H3E']=np.array([207,207+15,207+15+28,207+15+28+11])
	HFind['H4B']=np.array([55,66,66+27,66+27+11])
	HFind['H4F']=np.array([262,262+11,262+11+27,262+11+27+11])
	HFind['H2AC']=np.array([105,116,116+28,116+28+10])
	HFind['H2AG']=np.array([312,312+11,312+11+28,312+11+28+10])
	HFind['H2BD']=np.array([155,167,167+28,167+28+11])
	HFind['H2BH']=np.array([362,362+12,362+12+28,362+12+28+11])

	HFind['H3A']=HFind['H3A']*4+np.array([0,0,3,3])
	HFind['H4B']=HFind['H4B']*4+np.array([0,0,3,3])
	HFind['H2AC']=HFind['H2AC']*4+np.array([0,0,3,3])
	HFind['H2BD']=HFind['H2BD']*4+np.array([0,0,3,3])
	HFind['H3E']=HFind['H3E']*4+np.array([0,0,3,3])
	HFind['H4F']=HFind['H4F']*4+np.array([0,0,3,3])
	HFind['H2AG']=HFind['H2AG']*4+np.array([0,0,3,3])
	HFind['H2BH']=HFind['H2BH']*4+np.array([0,0,3,3])

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
	cross_corr = pickle.load( open( "../analysis_data/hfolds_covar.p", "rb" ) )
	# arange = np.arange(1656)
	# cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
	# cross_correlations[arange[0]+1:,arange[0]+1:] = cross_corr
	cross_corr=cross_corr[:,0:(cross_corr.shape[0]/2)]
	print cross_corr.shape
	fig4=plt.figure(figsize={7.5,11})
	sp=fig4.add_subplot(111)
	cax=sp.imshow(cross_corr,interpolation='none')
	cb=fig4.colorbar(cax,shrink=0.75)
	cb.set_label("$A^2$")

	#sp.set_axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
	sp.set_title('Cross-correlations')
	sp.set_xlabel('Atom indices')
	sp.set_ylabel('Atom indices')
	sp.set_title("Variance-covariance matrix for histone folds backbone.")
	# sp.set_xlim(0,100)

	# sp.set_xticks([HFind['H3A'].mean(),HFind['H4B'].mean(),HFind['H2AC'].mean(),HFind['H2BD'].mean(),HFind['H3E'].mean(),HFind['H4F'].mean(),HFind['H2AG'].mean(),HFind['H2BH'].mean()])
	# sp.set_xticklabels(['H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'])
	
	sp.set_xticks([HFind['H3A'].mean(),HFind['H4B'].mean(),HFind['H2AC'].mean(),HFind['H2BD'].mean()])
	sp.set_xticklabels(['H3 A','H4 B','H2A C','H2B D'])
	
	sp.set_yticks([HFind['H3A'].mean(),HFind['H4B'].mean(),HFind['H2AC'].mean(),HFind['H2BD'].mean(),HFind['H3E'].mean(),HFind['H4F'].mean(),HFind['H2AG'].mean(),HFind['H2BH'].mean()])
	sp.set_yticklabels(['H3 A','H4 B','H2A C','H2B D','H3 E','H4 F','H2A G','H2B H'])

	sp.tick_params(direction='out',length=6, bottom='off',top='off',left='off',right='off')
	
	sp.axvline(x=HFind['H3A'][0], color='yellow')
	sp.axvline(x=HFind['H3A'][-1], color='yellow')
	sp.axvline(x=HFind['H3A'][1],color='k')
	sp.axvline(x=HFind['H3A'][2],color='k')
	sp.axvline(x=HFind['H4B'][0], color='yellow')
	sp.axvline(x=HFind['H4B'][-1], color='yellow')
	sp.axvline(x=HFind['H4B'][1],color='k')
	sp.axvline(x=HFind['H4B'][2],color='k')
	sp.axvline(x=HFind['H2AC'][0] , color='yellow')
	sp.axvline(x=HFind['H2AC'][-1] , color='yellow')
	sp.axvline(x=HFind['H2AC'][1],color='k')
	sp.axvline(x=HFind['H2AC'][2],color='k')
	sp.axvline(x=HFind['H2BD'][0], color='yellow')
	sp.axvline(x=HFind['H2BD'][-1] , color='yellow')
	sp.axvline(x=HFind['H2BD'][1],color='k')
	sp.axvline(x=HFind['H2BD'][2],color='k')

	# sp.axvline(x=HFind['H3E'][0], color='yellow')
	# sp.axvline(x=HFind['H3E'][-1], color='yellow')
	# sp.axvline(x=HFind['H3E'][1],color='k')
	# sp.axvline(x=HFind['H3E'][2],color='k')
	# sp.axvline(x=HFind['H4F'][0], color='yellow')
	# sp.axvline(x=HFind['H4F'][-1], color='yellow')
	# sp.axvline(x=HFind['H4F'][1],color='k')
	# sp.axvline(x=HFind['H4F'][2],color='k')
	# sp.axvline(x=HFind['H2AG'][0] , color='yellow')
	# sp.axvline(x=HFind['H2AG'][-1] , color='yellow')
	# sp.axvline(x=HFind['H2AG'][1],color='k')
	# sp.axvline(x=HFind['H2AG'][2],color='k')
	# sp.axvline(x=HFind['H2BH'][0], color='yellow')
	# sp.axvline(x=HFind['H2BH'][-1] , color='yellow')
	# sp.axvline(x=HFind['H2BH'][1],color='k')
	# sp.axvline(x=HFind['H2BH'][2],color='k')

	sp.axhline(y=HFind['H3A'][0], color='yellow')
	sp.axhline(y=HFind['H3A'][-1], color='yellow')
	sp.axhline(y=HFind['H3A'][1],color='k')
	sp.axhline(y=HFind['H3A'][2],color='k')
	sp.axhline(y=HFind['H4B'][0], color='yellow')
	sp.axhline(y=HFind['H4B'][-1], color='yellow')
	sp.axhline(y=HFind['H4B'][1],color='k')
	sp.axhline(y=HFind['H4B'][2],color='k')
	sp.axhline(y=HFind['H2AC'][0] , color='yellow')
	sp.axhline(y=HFind['H2AC'][-1] , color='yellow')
	sp.axhline(y=HFind['H2AC'][1],color='k')
	sp.axhline(y=HFind['H2AC'][2],color='k')
	sp.axhline(y=HFind['H2BD'][0], color='yellow')
	sp.axhline(y=HFind['H2BD'][-1] , color='yellow')
	sp.axhline(y=HFind['H2BD'][1],color='k')
	sp.axhline(y=HFind['H2BD'][2],color='k')

	sp.axhline(y=HFind['H3E'][0], color='yellow')
	sp.axhline(y=HFind['H3E'][-1], color='yellow')
	sp.axhline(y=HFind['H3E'][1],color='k')
	sp.axhline(y=HFind['H3E'][2],color='k')
	sp.axhline(y=HFind['H4F'][0], color='yellow')
	sp.axhline(y=HFind['H4F'][-1], color='yellow')
	sp.axhline(y=HFind['H4F'][1],color='k')
	sp.axhline(y=HFind['H4F'][2],color='k')
	sp.axhline(y=HFind['H2AG'][0] , color='yellow')
	sp.axhline(y=HFind['H2AG'][-1] , color='yellow')
	sp.axhline(y=HFind['H2AG'][1],color='k')
	sp.axhline(y=HFind['H2AG'][2],color='k')
	sp.axhline(y=HFind['H2BH'][0], color='yellow')
	sp.axhline(y=HFind['H2BH'][-1] , color='yellow')
	sp.axhline(y=HFind['H2BH'][1],color='k')
	sp.axhline(y=HFind['H2BH'][2],color='k')

	fig4.tight_layout()
	fig4.savefig("../analysis_data/pca_hfolds_covar.png",dpi=(200))

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