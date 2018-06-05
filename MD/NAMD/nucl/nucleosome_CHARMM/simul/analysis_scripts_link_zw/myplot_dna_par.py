#!/usr/bin/env python

"""
Program for plotting x3dna param files
IMPORTANT: since termninal BP fray, the data often omits terminal BP,
this results in reference shift. So we have to correct for it!
Check this!

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv


def myplot(file, interactive):
	"Plot dat file"

	
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['ps.useafm'] = True
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['pdf.fonttype'] = 42

	ox=[]
	y=np.array([[1,1],[1,1]])

	# let's strip white spaces at first and replace with tabs
	tmpfile=file.replace('.dat','.tmp')
	with open(file,'r') as f:
		with open(tmpfile,'w') as tf:
			for line in f:
				tf.write('\t'.join(line.split())+'\n')
	with open(tmpfile,'r') as f:
		#f=open(file,'r')
		r=csv.reader(f,delimiter='\t')
		title=r.next().pop()
		description=r.next().pop()
		#xlabel=r.next().pop()
		#ylabel=r.next().pop()
		legends=r.next()
		#print legends
		first_row=r.next()
		ox.append(first_row[0])
		y=np.array(first_row[1:])
		for row in r:
				ox.append(row[0])
				y=np.vstack((y,row[1:]))
 		y=y.astype(float) # to make it of float numpy type
 		print ox
 		#ox=ox.astype(float)
 	plt.figure(figsize={9,10})
	plt.subplot(4,1,1)
 	
 	for i in [0,1,2,]:
		plt.plot(np.arange(len(ox)),y[:,i], linewidth=2, label=legends[i+1])

	plt.ylabel("Value")
	# plt.title("DNA parameters")
	plt.xlabel("BP number")
	plt.legend(loc="upper left",ncol=3)
	plt.ylim((-2,2))
	plt.xlim((0,147))

	plt.subplot(4,1,2)
	for i in [6,7,8]:
		plt.plot(np.arange(len(ox)),y[:,i], linewidth=2, label=legends[i+1])
	
	plt.ylabel("Value")
	# plt.title("DNA parameters")
	plt.xlabel("BP number")
	plt.legend(loc="upper left",ncol=3)
	plt.ylim((-2,2))
	plt.xlim((0,147))

	plt.subplot(4,1,3)
	for i in [3,4,5]:
		plt.plot(np.arange(len(ox)),y[:,i], linewidth=2, label=legends[i+1])

	plt.ylabel("Value")
	# plt.title("DNA parameters")
	plt.xlabel("BP number")
	plt.legend(loc="upper left",ncol=3)
	plt.ylim((-30,30))
	plt.xlim((0,147))

	plt.subplot(4,1,4)
	for i in [9,10,11]:
		plt.plot(np.arange(len(ox)),y[:,i], linewidth=2, label=legends[i+1])

	plt.ylabel("Value")
	# plt.title("DNA parameters")
	plt.xlabel("BP number")
	plt.legend(loc="upper left",ncol=3)
	plt.ylim((-30,30))
	plt.xlim((0,147))

	plt.tight_layout()
	plt.savefig(file.replace('.dat','.png'),dpi=(200))
	
	if(interactive):
		plt.show()




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('--int', default=0, help='Interactive')
	# myplot("dna_params/step_0.dat",1)
	args = parser.parse_args()
	myplot(args.file_name,args.int)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 