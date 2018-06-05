#!/usr/bin/env python

"""
Program for plotting dat files with legends and names

Special version to plot RMSF in nucleosome 1kx5

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import csv


def myplot(file, max_y,leg_loc):
	"Plot dat file"

	
	font = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 10}

	mpl.rc('font', **font)

#Let's load data
	with open(file,'r') as f:
		#f=open(file,'r')
		r=csv.reader(f,delimiter='\t')
		title=r.next().pop() # we need pop to get an element from a list
		description=r.next().pop()
		xlabel=r.next().pop()
		ylabel=r.next().pop()
		legends=r.next() # here we have a list
		#print legends
		first_row=r.next()

		ox=np.array(first_row[0])
		y=np.array(first_row[1:])
		for row in r:
				ox=np.append(ox,row[0])
				y=np.vstack((y,row[1:]))
 		y=y.astype(float) # to make it of float numpy type
 		ox=ox.astype(float)
 		

#reminder on histone fold
#set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"



 	ym= np.ma.masked_equal(y,0)
	fig=plt.figure(figsize={12,7})
 	sp=fig.add_subplot(321)

 	for i in [0,1]:
		sp.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

#set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
	sp.axvspan(64, 78, color='red', alpha=0.5)
	sp.axvspan(86, 114, color='red', alpha=0.5)
	sp.axvspan(121, 131, color='red', alpha=0.5)
	
	sp.set_ylabel(ylabel)
	sp.set_title("H3 histones")
	sp.set_xlabel(xlabel)
	sp.legend(loc=leg_loc,ncol=4)
	
	if(max_y!='auto'):
		sp.set_ylim((0,max_y[0]))

	sp2=fig.add_subplot(322)

 	for i in [2,3]:
		sp2.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

#set alpha_core_ca " (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
	sp2.axvspan(31, 41, color='red', alpha=0.5)
	sp2.axvspan(49, 76, color='red', alpha=0.5)
	sp2.axvspan(83, 93, color='red', alpha=0.5)

	sp2.set_ylabel(ylabel)
	sp2.set_title("H4 histones")
	sp2.set_xlabel(xlabel)
	sp2.legend(loc=leg_loc,ncol=4)

	if(max_y!='auto'):
		sp2.set_ylim((0,max_y[1]))

	sp3=fig.add_subplot(323)

 	for i in [4,5]:
		sp3.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

#set alpha_core_ca " (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
	sp3.axvspan(27, 37, color='red', alpha=0.5)
	sp3.axvspan(45, 73, color='red', alpha=0.5)
	sp3.axvspan(80, 89, color='red', alpha=0.5)

	sp3.set_ylabel(ylabel)
	sp3.set_title("H2A histones")
	sp3.set_xlabel(xlabel)
	sp3.legend(loc=leg_loc,ncol=4)

	if(max_y!='auto'):
		sp3.set_ylim((0,max_y[2]))

	sp4=fig.add_subplot(324)

 	for i in [6,7]:
		sp4.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])
#set alpha_core_ca "  or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
	sp4.axvspan(34, 45, color='red', alpha=0.5)
	sp4.axvspan(53, 81, color='red', alpha=0.5)
	sp4.axvspan(88, 98, color='red', alpha=0.5)

	#sp4.vlines([26,27],0,np.amax(y[:,6:7]),'g')
	sp4.set_ylabel(ylabel)
	sp4.set_title("H2B histones")
	sp4.set_xlabel(xlabel)
	sp4.legend(loc=leg_loc,ncol=4)

	if(max_y!='auto'):
		sp4.set_ylim((0,max_y[3]))

	sp5=fig.add_subplot(325)
	for i in [8,9]:
		sp5.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

	sp5.set_ylabel(ylabel)
	sp5.set_title("DNA")
	sp5.set_xlabel("Nucleotide number")
	sp5.legend(loc=leg_loc,ncol=4)

	if(max_y!='auto'):
		sp5.set_ylim((0,max_y[4]))

	sp6=fig.add_subplot(326)
	sp6.set_title("Boxplot")
	sp6.boxplot(y)
	sp6.set_ylabel(ylabel)
	#sp5.set_xticklabels([20,40],['a','b'])
	sp6.set_xticklabels(legends[1:])

	fig.tight_layout()

	#plt.boxplot(y)
	#plt.show()
	plt.subplots_adjust(top=0.9)
	fig.suptitle(title, fontsize=14)
	fig.savefig(file.replace('.dat','.png'),dpi=(200))





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('--max_y',dest='max_y',type=float, nargs=5, help='Maximum y for plots, an array of 5 numbers')
	parser.add_argument('--legend',dest='leg_loc',type=str, nargs=1, help='Legen location (def: upper left)')
	
	#myplot("../analysis_data/rmsf_chains.dat")
	args = parser.parse_args()
	if(not args.leg_loc):
		leg_loc='upper left'
	else:
		leg_loc=args.leg_loc[0]
	#print args.max_y
	if(args.max_y):
		myplot(args.file_name,args.max_y,leg_loc)
	else:
		myplot(args.file_name,'auto',leg_loc)

	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 