#!/usr/bin/env python

"""
Program for plotting DNA protein map

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

	y=y > 0

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
	fig=plt.figure(figsize={17,7})

	cmap = plt.get_cmap('cool')
	# cmap = plt.get_cmap('jet')

	h3=fig.add_subplot(141)
	pl=h3.imshow(DHmap[0],interpolation='nearest',cmap=cmap, extent=[DHresids[0][0],DHresids[0][-1],ox[0],ox[-1]],aspect="auto")
	h3.set_title("Histone H3")
	h3.set_xlabel("Histone resid")
	h3.set_ylabel("DNA resid")

	h4=fig.add_subplot(142)
	pl=h4.imshow(DHmap[1],interpolation='nearest',cmap=cmap,extent=[DHresids[1][0],DHresids[1][-1],ox[0],ox[-1]],aspect="auto")
	h4.set_title("Histone H4")
	h4.set_xlabel("Histone resid")
	h4.set_yticklabels([])
	# h3.set_ylabel("DNA resid")

	h2a=fig.add_subplot(143)
	pl=h2a.imshow(DHmap[2],interpolation='nearest',cmap=cmap,extent=[DHresids[2][0],DHresids[2][-1],ox[0],ox[-1]],aspect="auto")
	h2a.set_title("Histone H2A")
	h2a.set_xlabel("Histone resid")
	h2a.set_yticklabels([])
	# h3.set_ylabel("DNA resid")

	h2b=fig.add_subplot(144)
	pl=h2b.imshow(DHmap[3],interpolation='nearest',cmap=cmap,extent=[DHresids[3][0],DHresids[3][-1],ox[0],ox[-1]],aspect="auto")
	fig.colorbar(pl)
	h2b.set_yticklabels([])
	h2b.set_title("Histone H2B")
	h2b.set_xlabel("Histone resid")

	fig.suptitle("Presence of contacts between residues and DNA dinucl pairs", fontsize=18)
	# h3.set_ylabel("DNA resid")
	# plt.xticks(range(len(compresids)),compresids)

#  	ym= np.ma.masked_equal(y,0)
# 	fig=plt.figure(figsize={20,5})
#  	sp=fig.add_subplot(231)

#  	for i in [0,1]:
# 		sp.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

# #set alpha_core_ca "((segname CHA CHE and (resid 64 to 78 or resid 86 to 114 or resid 121 to 131)) or (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
# 	sp.axvspan(64, 78, color='red', alpha=0.5)
# 	sp.axvspan(86, 114, color='red', alpha=0.5)
# 	sp.axvspan(121, 131, color='red', alpha=0.5)
	
# 	sp.set_ylabel(ylabel)
# 	sp.set_title("H3 histones")
# 	sp.set_xlabel(xlabel)
# 	sp.legend(loc="upper left",ncol=4)
	
# 	sp2=fig.add_subplot(232)

#  	for i in [2,3]:
# 		sp2.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

# #set alpha_core_ca " (segname CHB CHF and (resid 31 to 41 or resid 49 to 76 or resid 83 to 93))  or (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
# 	sp2.axvspan(31, 41, color='red', alpha=0.5)
# 	sp2.axvspan(49, 76, color='red', alpha=0.5)
# 	sp2.axvspan(83, 93, color='red', alpha=0.5)

# 	sp2.set_ylabel(ylabel)
# 	sp2.set_title("H4 histones")
# 	sp2.set_xlabel(xlabel)
# 	sp2.legend(loc="upper left",ncol=4)

# 	sp3=fig.add_subplot(233)

#  	for i in [4,5]:
# 		sp3.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

# #set alpha_core_ca " (segname CHC CHG and (resid 27 to 37 or resid 45 to 73 or resid 80 to 89)) or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
# 	sp3.axvspan(27, 37, color='red', alpha=0.5)
# 	sp3.axvspan(45, 73, color='red', alpha=0.5)
# 	sp3.axvspan(80, 89, color='red', alpha=0.5)

# 	sp3.set_ylabel(ylabel)
# 	sp3.set_title("H2A histones")
# 	sp3.set_xlabel(xlabel)
# 	sp3.legend(loc="upper left",ncol=4)

# 	sp4=fig.add_subplot(234)

#  	for i in [6,7]:
# 		sp4.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])
# #set alpha_core_ca "  or (segname CHD CHH and (resid 34 to 45 or resid 53 to 81 or resid 88 to 98))) and name CA"
# 	sp4.axvspan(34, 45, color='red', alpha=0.5)
# 	sp4.axvspan(53, 81, color='red', alpha=0.5)
# 	sp4.axvspan(88, 98, color='red', alpha=0.5)

# 	sp4.vlines([26,27],0,np.amax(y[:,6:7]),'g')
# 	sp4.set_ylabel(ylabel)
# 	sp4.set_title("H2B histones")
# 	sp4.set_xlabel(xlabel)
# 	sp4.legend(loc="upper left",ncol=4)

# 	sp5=fig.add_subplot(235)
# 	for i in [8,9]:
# 		sp5.plot(ox,ym[:,i], linewidth=2, label=legends[i+1])

# 	sp5.set_ylabel(ylabel)
# 	sp5.set_title("DNA")
# 	sp5.set_xlabel(xlabel)
# 	sp5.legend(loc="upper left",ncol=4)

# 	sp6=fig.add_subplot(236)
# 	sp6.set_title("Boxplot")
# 	sp6.boxplot(y)
# 	sp6.set_ylabel(ylabel)
# 	#sp5.set_xticklabels([20,40],['a','b'])
# 	sp6.set_xticklabels(legends[1:])

	fig.tight_layout()
	fig.subplots_adjust(top=0.85)


# 	#plt.boxplot(y)
	if(interactive):
		plt.show()

	fig.savefig(file.replace('.dat','_has_cont.png'),dpi=(200))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('--int', default=0, help='Interactive')
	myplot("../analysis_data/dna_prot_contact_map.dat",1)
	#args = parser.parse_args()
	#myplot(args.file_name,args.int)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 