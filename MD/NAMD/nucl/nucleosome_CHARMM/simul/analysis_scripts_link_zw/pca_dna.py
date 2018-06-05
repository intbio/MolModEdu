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
import dnabb_cross_corr
import dnabb_covar
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
	dna_backbone='name P and resid `-73 to 73`'

	#Let's make correspondence table of histone borders


	dcd=parseDCD('../analysis_data/md_nucl.dcd',step=1) # get DCD as ensemble - load at once
	psf=parsePSF('../analysis_data/only_nucl_init.psf') # get PSF - to be able to get atomic names
	pdb_init=parsePDB('../analysis_data/only_nucl_init.pdb')
	pdb_aver=parsePDB('../analysis_data/only_nucl_average.pdb')
	dcd.setAtoms(psf.select(dna_backbone)) # this will determine what set to use for superposition and covariance
	dcd.setCoords(pdb_init) # we set reference coordinates to crystal structure - they wil be used to superimpose frames
	#dcd.superpose() - avoid DNA superposing

#Histone folds
	# print "Indices of selection"
	# print psf.select(alpha_core).getIndices()
	print "Num atoms:", len(psf.select(dna_backbone).getIndices())

	edaHF = EDA('DNA backbone analysis')
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
	fig1.suptitle("PCA analysis of DNA backbone movement: eigenvectors", fontsize=14)

	# plt.show()
	fig1.savefig("../analysis_data/pca_dnabb_eigvec.png",dpi=(200))


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
	fig2.suptitle("PCA analysis of DNA backbone movement: projections", fontsize=14)

	# plt.show()
	fig2.savefig("../analysis_data/pca_dnabb_project.png",dpi=(200))


#Distribution of modes along atoms
#(edaHF[0].getArrayNx3()**2).sum(axis=1)*edaHF[0].getVariance()
# this is what function calcSqFlucts - calculates - the variance in the mode contributed by individual atoms (as a sum of x,y and z components)

	fig3=plt.figure(figsize={9,10})

	sp=fig3.add_subplot(211)
	#we do some smoothing here with convolution
	#sp.set_color_cycle(['c', 'm', 'y', 'k'])
	for i in range(2):
		sp.plot(np.arange(0, float(edaHF.numAtoms()),1), calcSqFlucts(edaHF[i]), label="%d"%(i+1),linewidth=2)


		#SHL location#0 0.5 1.0 and so on
	SHL_loc=np.array([2,12,22,32,42,52,62,72,82,92,102,112,122,132,142])
	SHL_loc_ticks=np.array([-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7])

	SHL_loc_ticksx=np.concatenate((SHL_loc_ticks,SHL_loc_ticks))
	SHL_locx=np.concatenate((SHL_loc,SHL_loc+146))

	sp.set_xticks(SHL_locx)
	sp.set_xticklabels(SHL_loc_ticksx)

	sp.axvline(x=145.5,color='k')
	# sp.axvline(x=146,color='k')

	sp.annotate('Chain I', xy=(0.25, .5), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top',fontsize=20)
	sp.annotate('Chain J', xy=(.75, .5), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top',fontsize=20)


	sp.set_ylabel("Atomic fluct, $A^2$")
	sp.set_title("Mode decomposition on atomic fluctuations")
	sp.set_xlabel("SHL")
	sp.legend(loc="upper right", ncol=5)
	sp.set_xlim((0,edaHF.numAtoms()))


	sp2=fig3.add_subplot(212)
	#sp2.set_color_cycle(['c', 'm', 'y', 'k'])
	#we do some smoothing here with convolution
	for i in range(2,4):
		sp2.plot(np.arange(0, float(edaHF.numAtoms()),1), calcSqFlucts(edaHF[i]), label="%d"%(i+1),linewidth=2)


	sp2.set_ylabel("Atomic fluct, $A^2$")
	sp2.set_title("Mode decomposition on atomic fluctuations")
	sp2.set_xlabel("SHL")
	sp2.legend(loc="upper right", ncol=5)
	sp2.set_xlim((0,edaHF.numAtoms()))

	sp2.set_xticks(SHL_locx)
	sp2.set_xticklabels(SHL_loc_ticksx)

	sp2.axvline(x=145.5,color='k')
	# sp.axvline(x=146,color='k')

	sp2.annotate('Chain I', xy=(0.25, .0), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top',fontsize=20)
	sp2.annotate('Chain J', xy=(.75, .0), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top',fontsize=20)


	fig3.tight_layout()

	plt.subplots_adjust(top=0.9)
	fig3.suptitle("PCA analysis of DNA backbone movement: mode decomposition", fontsize=14)

	# plt.show()
	fig3.savefig("../analysis_data/pca_dnabb_modes.png",dpi=(200))

	# print calcFractVariance(edaHF).round(2)

	writeNMD('../analysis_data/eda_dnabb.nmd', edaHF[:5], pdb_init.select(dna_backbone))
	
	rmsd0=((np.amax(np.abs(HF_proj[:,0])))**2/edaHF.numAtoms())**0.5
	print "RMSD1:", rmsd0,"A"
	trj0=traverseMode(edaHF[0],pdb_aver.select(dna_backbone),rmsd=rmsd0,n_steps=40)
	writeDCD('../analysis_data/eda_dnabb_m1.dcd',trj0)

	rmsd1=((np.amax(np.abs(HF_proj[:,1])))**2/edaHF.numAtoms())**0.5
	print "RMSD2:", rmsd1,"A"
	trj1=traverseMode(edaHF[1],pdb_aver.select(dna_backbone),rmsd=rmsd1,n_steps=40)
	writeDCD('../analysis_data/eda_dnabb_m2.dcd',trj1)

	rmsd2=((np.amax(np.abs(HF_proj[:,2])))**2/edaHF.numAtoms())**0.5
	print "RMSD3:", rmsd2,"A"
	trj2=traverseMode(edaHF[2],pdb_aver.select(dna_backbone),rmsd=rmsd2,n_steps=40)
	writeDCD('../analysis_data/eda_dnabb_m3.dcd',trj2)

	rmsd3=((np.amax(np.abs(HF_proj[:,3])))**2/edaHF.numAtoms())**0.5
	print "RMSD4:", rmsd3,"A"
	trj3=traverseMode(edaHF[3],pdb_aver.select(dna_backbone),rmsd=rmsd3,n_steps=40)
	writeDCD('../analysis_data/eda_dnabb_m4.dcd',trj3)

#Now let's do covariance
	#let's dump it
	pickle.dump(calcCrossCorr(edaHF[0]), open( "../analysis_data/dnabb_cross_corr_m1.p", "wb" ) )
	pickle.dump(calcCrossCorr(edaHF[1]), open( "../analysis_data/dnabb_cross_corr_m2.p", "wb" ) )
	pickle.dump(calcCrossCorr(edaHF), open( "../analysis_data/dnabb_cross_corr.p", "wb" ) )

	pickle.dump(mycalcCovariance(edaHF[0]), open( "../analysis_data/dnabb_covar_m1.p", "wb" ) )
	pickle.dump(mycalcCovariance(edaHF[1]), open( "../analysis_data/dnabb_covar_m2.p", "wb" ) )
	pickle.dump(mycalcCovariance(edaHF), open( "../analysis_data/dnabb_covar.p", "wb" ) )
	
	dnabb_cross_corr.main()
	dnabb_covar.main()

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