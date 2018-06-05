#!/usr/bin/env python

"""
Program for plotting dat files with legends and names

Special version to plot RMSF in nucleosome 1kx5 only for DNA

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import csv
import re
from collections import defaultdict


def myplot(file_bb,file_sch,pdbfile,max_y,leg_loc):
	"Plot dat file"

	
	font = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 10}

	mpl.rc('font', **font)

#Let's load data

	with open(file_bb,'r') as f:
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

	with open(file_sch,'r') as f:
		#f=open(file,'r')
		r=csv.reader(f,delimiter='\t')
		title_sch=r.next().pop() # we need pop to get an element from a list
		description_sch=r.next().pop()
		xlabel_sch=r.next().pop()
		ylabel_sch=r.next().pop()
		legends_sch=r.next() # here we have a list
		#print legends
		first_row=r.next()

		ox_sch=np.array(first_row[0])
		y_sch=np.array(first_row[1:])
		for row in r:
				ox_sch=np.append(ox_sch,row[0])
				y_sch=np.vstack((y_sch,row[1:]))
		y_sch=y_sch.astype(float) # to make it of float numpy type
		ox_sch=ox_sch.astype(float)
 
	exp_r_sch=defaultdict(lambda: [[],[]])
	with open(pdbfile,'r') as f:
		#f=open(file,'r')
#ATOM  13309  CB  LYS H  82      30.557  67.488   3.180  1.00 44.63           C
		
####!!!!Correct PDB so that there is a space before B-factor
		r=re.compile('ATOM\s+\d+\s+(\S+)\s+\S+\s+(\S)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)')
		resid_old=100000 #fake number, which will not match ever
		mean=0
		numat=0
		for line in f:  #let's do averaging here
			if(r.match(line)):
				atomname=r.match(line).group(1)
				if(not re.match('^CA$|^C$|^O$|^N$|^P$|^OP1$|^OP2$|^C1\'$|^C2\'$|^C3\'$|^C4\'$|^C5\'$|^O3\'$|^O4\'$|^O5\'$',atomname)):
					# print "match"
					resid=r.match(line).group(3)
					if((resid==resid_old) or (mean==0)):
						# print r.match(line).group(4)
						mean=mean+float(r.match(line).group(4))**0.5
						numat=numat+1
						chainid=r.match(line).group(2)
					else:
						mean=mean/numat
						exp_r_sch[chainid][0].append(resid_old)
						exp_r_sch[chainid][1].append(mean)
						mean=float(r.match(line).group(4))**0.5
						chainid=r.match(line).group(2)
						numat=1

				# exp_r[chainid][0].append(r.match(line).group(2))
				# exp_r[chainid][1].append(r.match(line).group(3))
					resid_old=resid
				#print r.match(line).group(1), r.match(line).group(2), r.match(line).group(3)
	
	exp_rn_sch=defaultdict(list)
	for key in exp_r_sch:
		exp_rn_sch[key]=np.array([exp_r_sch[key][0],exp_r_sch[key][1]])
		exp_rn_sch[key]=exp_rn_sch[key].astype(float)
		#what we have is average of square roots of B-factors
		#since B=8*pi*pi*RMSF^2 , we have averaged the quantity sqrt(8)*pi*RMSF
		exp_rn_sch[key][1,:]=(exp_rn_sch[key][1,:]/(8**0.5)/3.14)
		# exp_rn[key][1,:]=(exp_rn[key][1,:]/8/3.14/3.14)**0.5

	print exp_rn_sch


	exp_r=defaultdict(lambda: [[],[]])
	with open(pdbfile,'r') as f:
		#f=open(file,'r')
#ATOM  13309  CB  LYS H  82      30.557  67.488   3.180  1.00 44.63           C
		r=re.compile('ATOM\s+\d+\s+(?:CA|P)\s+\S+\s+(\S)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)')
		for line in f:
			if(r.match(line)):
				exp_r[r.match(line).group(1)][0].append(r.match(line).group(2))
				exp_r[r.match(line).group(1)][1].append(r.match(line).group(3))
				#print r.match(line).group(1), r.match(line).group(2), r.match(line).group(3)
	
	exp_rn=defaultdict(list)
	for key in exp_r:
		exp_rn[key]=np.array([exp_r[key][0],exp_r[key][1]])
		exp_rn[key]=exp_rn[key].astype(float)
		exp_rn[key][1,:]=(exp_rn[key][1,:]/8/3.14/3.14)**0.5
	






 	ym=np.ma.masked_equal(y,0)
 	ym_sch=np.ma.masked_equal(y_sch,0)

	fig=plt.figure(figsize={12,7})
 	sp=fig.add_subplot(211)

	for i in [8,9]:
		sp.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

	chain='I'
	sp.plot(exp_rn[chain][0,:],exp_rn[chain][1,:]*1.2, linewidth=2, label="Exp "+chain)
	chain='J'
	sp.plot(exp_rn[chain][0,:],exp_rn[chain][1,:]*1.2, linewidth=2, label="Exp "+chain)



	sp.set_ylabel(ylabel)
	sp.set_title("DNA")
	sp.set_xlabel(xlabel)
	sp.legend(loc=leg_loc,ncol=4)

	if(max_y!='auto'):
		sp.set_ylim((0.5,max_y[0]))
		sp.set_xlim((-80,80))




	sp5=fig.add_subplot(212)
	for i in [8,9]:
		sp5.plot(ox_sch,ym_sch[:,i], linewidth=2, label=legends_sch[i+1])

	chain='I'
	sp5.plot(exp_rn_sch[chain][0,:],exp_rn_sch[chain][1,:]*1.2, linewidth=2, label="Exp "+chain)
	chain='J'
	sp5.plot(exp_rn_sch[chain][0,:],exp_rn_sch[chain][1,:]*1.2, linewidth=2, label="Exp "+chain)



	sp5.set_ylabel(ylabel_sch)
	sp5.set_title("DNA")
	sp5.set_xlabel(xlabel_sch)
	sp5.legend(loc=leg_loc,ncol=4)

	if(max_y!='auto'):
		sp5.set_ylim((0.5,max_y[1]))
		sp5.set_xlim((-80,80))



	fig.tight_layout()

	plt.subplots_adjust(top=0.93)
	fig.suptitle("Comparison of calculated RMSF with B-factor derived RMSF for DNA, B-factors are scaled by 1.2", fontsize=14)

	# plt.show()
	fig.savefig(file_bb.replace('.dat','_dna_exp.png'),dpi=(200))





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file for backbone ')
	parser.add_argument('file_name2', help='Name of dat file for side chians')
	parser.add_argument('--max_y',dest='max_y',type=float, nargs=2, help='Maximum y for plots, an array of 2 numbers')
	parser.add_argument('--legend',dest='leg_loc',type=str,  help='Legen location (def: upper left)')
	
	args = parser.parse_args()
	# myplot(args.file_name,"../analysis_data/1kx5.pdb")
	if(not hasattr(args, 'leg_loc')):
		leg_loc='upper left'
	else:
		leg_loc=args.leg_loc
	if(args.max_y):
		myplot(args.file_name,args.file_name2,"../analysis_data/1kx5.pdb",args.max_y,leg_loc)
	else:
		myplot(args.file_name,args.file_name2,"../analysis_data/1kx5.pdb",'auto',leg_loc)
	#myplot(args.file_name)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 