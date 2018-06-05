#!/usr/bin/env python

"""
Program for plotting DNA protein contact profiles

Copyright 2013 (c) Alexey Shaytan, NCBI/NLM/NIH

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import csv
from collections import defaultdict




def myplot(file, interactive):
	"Plot dat file"

	
	font = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 14}

	mpl.rc('font', **font)

#Let's load data
	with open(file,'r') as f:
		#f=open(file,'r')
		r=csv.reader(f,delimiter='\t')
		title=r.next().pop() # we need pop to get an element from a list
		description=r.next().pop()
		xlabel=r.next().pop()
		ylabel=r.next().pop()
		resids=np.array(r.next().pop().split()) # here we have a list
		resnames=r.next().pop().split() # here we have a list
		#print legends
		first_row=r.next()

		ox=np.array(first_row[0])
		y=np.array(first_row[1:])
		for row in r:
				ox=np.append(ox,row[0])
				y=np.vstack((y,row[1:]))
 		y=y.astype(float) # to make it of float numpy type
 		ox=ox.astype(int)
 		resids=resids.astype(int)
 		
 		# print resids


# Let's now process the data, determine histone borders
#
#
#
	hist_list=['H3_A','H3_E','H4_B','H4_F','H2A_C','H2A_G','H2B_D','H2B_H']
	hist_border = defaultdict(list)

	ic=0
	hist_border[hist_list[ic]].append(0)
	for i in range(1,len(resids)):
		if(resids[i]<resids[i-1]):
			hist_border[hist_list[ic]].append(i-1)
			ic+=1
			hist_border[hist_list[ic]].append(i)
	hist_border[hist_list[ic]].append(i)
	print hist_border

#Let's do some data rearrangement

#Let's compress twice summing for each histone and separate by histone
	#DHmap=np.empty((y.shape[0],y.shape[1]/2))
	#DHresids=np.empty(y.shape[1]/2)

	print y.shape
	# print comp.shape
	DHmap=list()
	DHresids=list()

	# y=y > 0

	for p in range(4):
		i=p*2
		DHmap.append(np.empty((len(y[:,0]),hist_border[hist_list[i]][1]-hist_border[hist_list[i]][0]+1)))
		DHresids.append(list())
		cn=0
		for k in range(hist_border[hist_list[i]][0],hist_border[hist_list[i]][1]+1):
			DHmap[p][:,cn]=y[:,k]+y[:,k-hist_border[hist_list[i]][0]+hist_border[hist_list[i+1]][0]]
			DHresids[p].append(resids[k])
			cn+=1
	# print DHresids
	# y=np.transpose(y)

	#Let's now collapse data to get profiles
	Pprof=list()
	for p in range(4):
		Pprof.append(np.sum(DHmap[p],axis=0))
	Dprof=np.sum(y,axis=1)

	fig=plt.figure(figsize={10,5})

	# cmap = plt.get_cmap('cool')
	# cmap = plt.get_cmap('jet')

	h3=fig.add_subplot(111)
	h3.plot(ox,Dprof,linewidth=2)
	h3.set_title("DNA prtoien contact profile")
	h3.set_xlabel("DNA basepair")
	h3.set_ylabel("Number of contacts")

	fig.tight_layout()


# 	#plt.boxplot(y)
	if(interactive):
		plt.show()

	fig.savefig(file.replace('.dat','_dna_prof.png'),dpi=(200))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('--int', default=0, help='Interactive')
	myplot("../analysis_data/dna_prot_contact_map.dat",1)
	#args = parser.parse_args()
	#myplot(args.file_name,args.int)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 