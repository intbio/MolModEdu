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

	fig=plt.figure(figsize={8,10})

	# cmap = plt.get_cmap('cool')
	# cmap = plt.get_cmap('jet')

	h3=fig.add_subplot(411)
	h3.plot(DHresids[0],Pprof[0],linewidth=2)
	h3.set_title("Histone H3")
	h3.set_xlabel("Histone resid")
	h3.set_ylabel("Number of contacts")
	h3.axvspan(64, 78, color='red', alpha=0.5)
	h3.axvspan(86, 114, color='red', alpha=0.5)
	h3.axvspan(121, 131, color='red', alpha=0.5)

	h4=fig.add_subplot(412)
	h4.plot(DHresids[1],Pprof[1],linewidth=2)
	h4.set_title("Histone H4")
	h4.set_xlabel("Histone resid")
	h4.set_ylabel("Number of contacts")
	h4.axvspan(31, 41, color='red', alpha=0.5)
	h4.axvspan(49, 76, color='red', alpha=0.5)
	h4.axvspan(83, 93, color='red', alpha=0.5)

	h2a=fig.add_subplot(413)
	h2a.plot(DHresids[2],Pprof[2],linewidth=2)
	h2a.set_title("Histone H2A")
	h2a.set_xlabel("Histone resid")
	h2a.set_ylabel("Number of contacts")
	h2a.axvspan(27, 37, color='red', alpha=0.5)
	h2a.axvspan(45, 73, color='red', alpha=0.5)
	h2a.axvspan(80, 89, color='red', alpha=0.5)

	h2b=fig.add_subplot(414)
	h2b.plot(DHresids[3],Pprof[3],linewidth=2)
	h2b.set_title("Histone H2B")
	h2b.set_xlabel("Histone resid")
	h2b.set_ylabel("Number of contacts")
	h2b.axvspan(34, 45, color='red', alpha=0.5)
	h2b.axvspan(53, 81, color='red', alpha=0.5)
	h2b.axvspan(88, 98, color='red', alpha=0.5)

	
	fig.suptitle("Protein DNA contact profiles", fontsize=18)

	fig.tight_layout()
	fig.subplots_adjust(top=0.90)


# 	#plt.boxplot(y)
	if(interactive):
		plt.show()

	fig.savefig(file.replace('.dat','_prot_prof.png'),dpi=(200))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('--int', default=0, help='Interactive')
	myplot("../analysis_data/dna_prot_contact_map.dat",1)
	#args = parser.parse_args()
	#myplot(args.file_name,args.int)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 