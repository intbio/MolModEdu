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
import pickle
import hfolds_cross_corr
import hfolds_covar
from myPrody import mycalcCovariance

#import prody as pr
#import scipy
EVNUM=40

def main():
	

	
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['ps.useafm'] = True
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['pdf.fonttype'] = 42
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

	dcd=parseDCD('../analysis_data/md_nucl.dcd',step=1) # get DCD as ensemble - load at once
	psf=parsePSF('../analysis_data/only_nucl_init.psf') # get PSF - to be able to get atomic names
	pdb_init=parsePDB('../analysis_data/only_nucl_init.pdb')
	pdb_aver=parsePDB('../analysis_data/only_nucl_average.pdb')
	dcd.setAtoms(psf.select(alpha_core)) # this will determine what set to use for superposition and covariance
	dcd.setCoords(pdb_init) # we set reference coordinates to crystal structure - they wil be used to superimpose frames
	dcd.superpose()

#Histone folds
	# print "Indices of selection"
	# print psf.select(alpha_core).getIndices()
	print "Num atoms:", len(psf.select(alpha_core).getIndices())

	edaHF = EDA('Histone fold analysis')
	#We want to build covariance - since this is dcd object, the frames will be superimposed on crystal using 
	edaHF.buildCovariance( dcd[0:],aligned=True )
	edaHF.calcModes(n_modes=EVNUM)
	print "Covariance matrix"
	print edaHF.getCovariance() # returns covariance matrix
	print "Eigenvalues in A^2"
	print edaHF.getEigvals() # returns eigenvalues in A^2
	print "Eigenvalues in A"
	print edaHF.getEigvals()**0.5 # returns eigenvalues in A^2
	print "Fraction of variance"
	print calcFractVariance(edaHF) # returns fraction of variance, PRODY knows the full trace
	print "Mode covariance"
	print calcCovariance(edaHF[0])
#Study eigenvectors

	#Let's make the eigenvector plot
	fig1=plt.figure(figsize={9,10})

	sp=fig1.add_subplot(411)
	sp.bar(range(EVNUM),edaHF.getEigvals(), linewidth=2, label="$A^2$")

	sp.set_ylabel("Eigenvalue magnitude, $A^2$")
	sp.set_title("Eigenvalues in $A^2$")
	sp.set_xlabel("Index")
	sp.legend(loc="upper right", ncol=4)

	sp2=fig1.add_subplot(412)
	sp2.bar(range(EVNUM),edaHF.getEigvals()**0.5, linewidth=2, label="$A$")

	sp2.set_ylabel("Eigenvalue magnitude, $A$")
	sp2.set_title("Eigenvalues in $A$")
	sp2.set_xlabel("Index")
	sp2.legend(loc="upper right", ncol=4)

	sp3=fig1.add_subplot(413)
	sp3.bar(range(EVNUM),calcFractVariance(edaHF), linewidth=2, label="Fraction")
	sp3.plot(range(EVNUM),np.cumsum(calcFractVariance(edaHF)), linewidth=2, color='r', label="Cumulative")

	sp3.set_ylabel("Fraction, %")
	sp3.set_title("Fraction of total variance")
	sp3.set_xlabel("Index")
	sp3.legend(loc="lower right", ncol=4)

	#minimax plot
	HF_proj=calcProjection(dcd[0:], edaHF,rmsd=False)

	pmax=np.amax(HF_proj,axis=0)
	pmin=np.amin(HF_proj,axis=0)

	sp4=fig1.add_subplot(414)
	sp4.bar(range(EVNUM),pmax, linewidth=2, label="Max",width=0.4,color='r')
	sp4.bar(np.arange(EVNUM)+0.4,-pmin, linewidth=2, label="-Min",width=0.4,color='b')
	sp4.bar(np.arange(EVNUM)+0.8,pmax-pmin, linewidth=2, label="Max-min",width=0.2,color='g')

	sp4.set_ylabel("Min max, $A$")
	sp4.set_title("Min and max projection on the mode")
	sp4.set_xlabel("Index")
	sp4.legend(loc="upper right", ncol=4)

	fig1.tight_layout()

	plt.subplots_adjust(top=0.9)
	fig1.suptitle("PCA analysis of histone folds movement: eigenvectors", fontsize=14)

	# plt.show()
	fig1.savefig("../analysis_data/pca_hfolds_eigvec.png",dpi=(200))


#Study projections

	HF_proj=calcProjection(dcd[0:], edaHF[0:5],rmsd=False)
	print "Proj max", np.amax(HF_proj,axis=0)
	print "Proj min", np.amin(HF_proj,axis=0)
	print "Proj aver", np.average(HF_proj,axis=0)


	print "Collect", calcCollectivity(edaHF[0]),calcCollectivity(edaHF[1]),calcCollectivity(edaHF[2]),calcCollectivity(edaHF[3]),calcCollectivity(edaHF[4]),calcCollectivity(edaHF[5]),calcCollectivity(edaHF[6]),calcCollectivity(edaHF[7]),calcCollectivity(edaHF[8]),calcCollectivity(edaHF[9]),calcCollectivity(edaHF[10]),calcCollectivity(edaHF[11]),calcCollectivity(edaHF[12]),calcCollectivity(edaHF[13]),calcCollectivity(edaHF[14]),calcCollectivity(edaHF[15]),calcCollectivity(edaHF[16])
#Let's make the projection plot

	fig2=plt.figure(figsize={9,10})

	sp=fig2.add_subplot(211)
	#we do some smoothing here with convolution
	kernel=np.array([0.01]*100)
	for i in range(2):
		sp.plot(np.arange(0, float(len(HF_proj[:,i]))/10,0.1), np.convolve(HF_proj[:,i],kernel,mode='same'), label="%d"%(i+1),linewidth=2)


	sp.set_ylabel("Projection, $A$")
	sp.set_title("Projections onto eigenvectors (10 ns smoothing)")
	sp.set_xlabel("Time, ns")
	sp.legend(loc="upper right", ncol=5)

	sp2=fig2.add_subplot(212)
	#we do some smoothing here with convolution
	for i in range(2,4):
		sp2.plot(np.arange(0, float(len(HF_proj[:,i]))/10,0.1), np.convolve(HF_proj[:,i],kernel,mode='same'), label="%d"%(i+1),linewidth=2)


	sp2.set_ylabel("Projection, $A$")
	sp2.set_title("Projections onto eigenvectors (10 ns smoothing)")
	sp2.set_xlabel("Time, ns")
	sp2.legend(loc="upper right", ncol=5)

	fig2.tight_layout()

	plt.subplots_adjust(top=0.9)
	fig2.suptitle("PCA analysis of histone folds movement: projections", fontsize=14)

	# plt.show()
	fig2.savefig("../analysis_data/pca_hfolds_project.png",dpi=(200))


#Distribution of modes along atoms
#(edaHF[0].getArrayNx3()**2).sum(axis=1)*edaHF[0].getVariance()
# this is what function calcSqFlucts - calculates - the variance in the mode contributed by individual atoms (as a sum of x,y and z components)

	fig3=plt.figure(figsize={9,10})

	sp=fig3.add_subplot(211)
	#we do some smoothing here with convolution
	#sp.set_color_cycle(['c', 'm', 'y', 'k'])
	for i in range(2):
		sp.plot(np.arange(0, float(edaHF.numAtoms()),1), calcSqFlucts(edaHF[i]), label="%d"%(i+1),linewidth=2)
		sp.axvspan(HFind['H3A'][0],HFind['H3A'][-1] , color='blue', alpha=0.2)
		sp.axvline(x=HFind['H3A'][1],color='k')
		sp.axvline(x=HFind['H3A'][2],color='k')
		sp.axvspan(HFind['H4B'][0],HFind['H4B'][-1] , color='green', alpha=0.2)
		sp.axvline(x=HFind['H4B'][1],color='k')
		sp.axvline(x=HFind['H4B'][2],color='k')
		sp.axvspan(HFind['H2AC'][0],HFind['H2AC'][-1] , color='yellow', alpha=0.2)
		sp.axvline(x=HFind['H2AC'][1],color='k')
		sp.axvline(x=HFind['H2AC'][2],color='k')
		sp.axvspan(HFind['H2BD'][0],HFind['H2BD'][-1] , color='red', alpha=0.2)
		sp.axvline(x=HFind['H2BD'][1],color='k')
		sp.axvline(x=HFind['H2BD'][2],color='k')

		sp.axvspan(HFind['H3E'][0],HFind['H3E'][-1] , color='blue', alpha=0.2)
		sp.axvline(x=HFind['H3E'][1],color='k')
		sp.axvline(x=HFind['H3E'][2],color='k')
		sp.axvspan(HFind['H4F'][0],HFind['H4F'][-1] , color='green', alpha=0.2)
		sp.axvline(x=HFind['H4F'][1],color='k')
		sp.axvline(x=HFind['H4F'][2],color='k')
		sp.axvspan(HFind['H2AG'][0],HFind['H2AG'][-1] , color='yellow', alpha=0.2)
		sp.axvline(x=HFind['H2AG'][1],color='k')
		sp.axvline(x=HFind['H2AG'][2],color='k')
		sp.axvspan(HFind['H2BH'][0],HFind['H2BH'][-1] , color='red', alpha=0.2)
		sp.axvline(x=HFind['H2BH'][1],color='k')
		sp.axvline(x=HFind['H2BH'][2],color='k')


	sp.set_ylabel("Atomic fluct, $A^2$")
	sp.set_title("Mode decomposition on atomic fluctuations")
	sp.set_xlabel("Atom index")
	sp.legend(loc="upper right", ncol=5)
	sp.set_xlim((0,edaHF.numAtoms()))


	sp2=fig3.add_subplot(212)
	#sp2.set_color_cycle(['c', 'm', 'y', 'k'])
	#we do some smoothing here with convolution
	for i in range(2,4):
		sp2.plot(np.arange(0, float(edaHF.numAtoms()),1), calcSqFlucts(edaHF[i]), label="%d"%(i+1),linewidth=2)
		sp2.axvspan(HFind['H3A'][0],HFind['H3A'][-1] , color='blue', alpha=0.2)
		sp2.axvline(x=HFind['H3A'][1],color='k')
		sp2.axvline(x=HFind['H3A'][2],color='k')
		sp2.axvspan(HFind['H4B'][0],HFind['H4B'][-1] , color='green', alpha=0.2)
		sp2.axvline(x=HFind['H4B'][1],color='k')
		sp2.axvline(x=HFind['H4B'][2],color='k')
		sp2.axvspan(HFind['H2AC'][0],HFind['H2AC'][-1] , color='yellow', alpha=0.2)
		sp2.axvline(x=HFind['H2AC'][1],color='k')
		sp2.axvline(x=HFind['H2AC'][2],color='k')
		sp2.axvspan(HFind['H2BD'][0],HFind['H2BD'][-1] , color='red', alpha=0.2)
		sp2.axvline(x=HFind['H2BD'][1],color='k')
		sp2.axvline(x=HFind['H2BD'][2],color='k')

		sp2.axvspan(HFind['H3E'][0],HFind['H3E'][-1] , color='blue', alpha=0.2)
		sp2.axvline(x=HFind['H3E'][1],color='k')
		sp2.axvline(x=HFind['H3E'][2],color='k')
		sp2.axvspan(HFind['H4F'][0],HFind['H4F'][-1] , color='green', alpha=0.2)
		sp2.axvline(x=HFind['H4F'][1],color='k')
		sp2.axvline(x=HFind['H4F'][2],color='k')
		sp2.axvspan(HFind['H2AG'][0],HFind['H2AG'][-1] , color='yellow', alpha=0.2)
		sp2.axvline(x=HFind['H2AG'][1],color='k')
		sp2.axvline(x=HFind['H2AG'][2],color='k')
		sp2.axvspan(HFind['H2BH'][0],HFind['H2BH'][-1] , color='red', alpha=0.2)
		sp2.axvline(x=HFind['H2BH'][1],color='k')
		sp2.axvline(x=HFind['H2BH'][2],color='k')

	sp2.set_ylabel("Atomic fluct, $A^2$")
	sp2.set_title("Mode decomposition on atomic fluctuations")
	sp2.set_xlabel("Atom index")
	sp2.legend(loc="upper right", ncol=5)
	sp2.set_xlim((0,edaHF.numAtoms()))

	fig3.tight_layout()

	plt.subplots_adjust(top=0.9)
	fig3.suptitle("PCA analysis of histone folds movement: mode decomposition", fontsize=14)

	# plt.show()
	fig3.savefig("../analysis_data/pca_hfolds_modes.png",dpi=(200))

	# print calcFractVariance(edaHF).round(2)

	writeNMD('../analysis_data/eda_hfolds.nmd', edaHF[:5], pdb_init.select(alpha_core))
	
	rmsd0=((np.amax(np.abs(HF_proj[:,0])))**2/edaHF.numAtoms())**0.5
	print "RMSD1:", rmsd0,"A"
	trj0=traverseMode(edaHF[0],pdb_aver.select(alpha_core),rmsd=rmsd0,n_steps=40)
	writeDCD('../analysis_data/eda_hfolds_m1.dcd',trj0)

	rmsd1=((np.amax(np.abs(HF_proj[:,1])))**2/edaHF.numAtoms())**0.5
	print "RMSD2:", rmsd1,"A"
	trj1=traverseMode(edaHF[1],pdb_aver.select(alpha_core),rmsd=rmsd1,n_steps=40)
	writeDCD('../analysis_data/eda_hfolds_m2.dcd',trj1)

	rmsd2=((np.amax(np.abs(HF_proj[:,2])))**2/edaHF.numAtoms())**0.5
	print "RMSD3:", rmsd2,"A"
	trj2=traverseMode(edaHF[2],pdb_aver.select(alpha_core),rmsd=rmsd2,n_steps=40)
	writeDCD('../analysis_data/eda_hfolds_m3.dcd',trj2)

	rmsd3=((np.amax(np.abs(HF_proj[:,3])))**2/edaHF.numAtoms())**0.5
	print "RMSD4:", rmsd3,"A"
	trj3=traverseMode(edaHF[3],pdb_aver.select(alpha_core),rmsd=rmsd3,n_steps=40)
	writeDCD('../analysis_data/eda_hfolds_m4.dcd',trj3)

#Now let's do covariance
	#let's dump it
	pickle.dump(calcCrossCorr(edaHF[0]), open( "../analysis_data/hfolds_cross_corr_m1.p", "wb" ) )
	pickle.dump(calcCrossCorr(edaHF[1]), open( "../analysis_data/hfolds_cross_corr_m2.p", "wb" ) )
	pickle.dump(calcCrossCorr(edaHF), open( "../analysis_data/hfolds_cross_corr.p", "wb" ) )
	
	pickle.dump(mycalcCovariance(edaHF[0]), open( "../analysis_data/hfolds_covar_m1.p", "wb" ) )
	pickle.dump(mycalcCovariance(edaHF[1]), open( "../analysis_data/hfolds_covar_m2.p", "wb" ) )
	pickle.dump(mycalcCovariance(edaHF), open( "../analysis_data/hfolds_covar.p", "wb" ) )
	

	hfolds_cross_corr.main()
	hfolds_covar.main()

	# arange = np.arange(edaHF.numAtoms())
	# cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
	# cross_correlations[arange[0]+1:,arange[0]+1:] = calcCrossCorr(edaHF)
	# fig4=plt.figure(figsize={9,10})
	# sp=fig4.add_subplot(111)
	# cax=sp.imshow(cross_correlations)
	# fig4.colorbar(cax)
	# #sp.set_axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
	# sp.set_title('Cross-correlations')
	# sp.set_xlabel('Atom indices')
	# sp.set_ylabel('Atom indices')
	# sp.set_title("Variance covariance matrix for histone folds backbone")
	# fig4.savefig("../analysis_data/pca_hfolds_var-covar.png",dpi=(200))

	# dcd.setAtoms(pdb.select("name P"))
	# edaDNA = EDA('DNA')
	# edaDNA.buildCovariance( dcd[2500:] )
	# edaDNA.calcModes()
	# print calcFractVariance(edaDNA).round(4)
	# print calcCovariance(edaDNA).round(4)
	# DNA_proj=calcProjection(dcd[0:], edaDNA[:2]);
	# showProjection(dcd[0:],edaDNA[:3])

	# figure(figsize={12,5})

	

	# subplot(121)
	# plot( DNA_proj[:,0], DNA_proj[:,1],'bo', label="DNA P")
	# plot(HF_proj[:,0],HF_proj[:,1],'ro', label="Histone fold CA")
	
	# xlabel("Principle mode 1, A")
	# ylabel("Principle mode 2, A")
	# legend(loc="upper right")
	# legend(numpoints=1)
	# title("Principle components of motion in nucleosome")

	# subplot(122)

	# plot(np.arange(0, float(len(HF_proj[:,0]))/10,0.1), DNA_proj[:,0],'b-', label="DNA P, mode 1")
	# plot(np.arange(0, float(len(HF_proj[:,0]))/10,0.1), DNA_proj[:,1],'g-',linewidth=2, label="DNA P, mode 2")
	# plot(np.arange(0, float(len(HF_proj[:,0]))/10,0.1), HF_proj[:,0],'r-', label="Histone fold CA, mode 1")
	# plot(np.arange(0, float(len(HF_proj[:,0]))/10,0.1), HF_proj[:,1],'y-',linewidth=2, label="Histone fold CA, mode 2")
	

	# xlabel("Time, ns")
	# ylabel("Principle mode value, A")
	# legend(loc="upper right")
	# legend(numpoints=1)
	# title("Principle components of motion in nucleosome")

	# tight_layout()
	# #show()
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