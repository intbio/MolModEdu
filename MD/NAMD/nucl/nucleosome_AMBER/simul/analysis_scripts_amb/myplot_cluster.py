#!/usr/bin/env python

"""
Program for plotting dat files with legends and names

this version tuned to display vertical lines supplied by data points from second file
Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv


def myplot(file, file2, x_r, interactive):
	"Plot dat file"

	
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['ps.useafm'] = True
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['pdf.fonttype'] = 42

	clustp=list()
	with open(file2,'r') as f2:
		r=csv.reader(f2,delimiter=' ')
		for row in r:
			clustp.append(np.array(row,dtype=float)/1) #2000 to recalulate frame numbers to time!




	ox=[]
	y=np.array([[1,1],[1,1]])
	with open(file,'r') as f:
		#f=open(file,'r')
		r=csv.reader(f,delimiter='\t')
		title=r.next().pop()
		description=r.next().pop()
		xlabel=r.next().pop()
		ylabel=r.next().pop()
		legends=r.next()
		#print legends
		first_row=r.next()
		ox.append(first_row[0])
		y=np.array(first_row[1:])
		for row in r:
				ox.append(row[0])
				y=np.vstack((y,row[1:]))
		y=y.astype(float) # to make it of float numpy type
 		

 	plt.figure(figsize={12,6})
 	for i in range(len(legends)-1):
		plt.plot(ox,y[:,i], linewidth=2, label=legends[i+1])

	color=['#00ff00','#ffff00','#00ffff','#ff00ff','#009999','#990099']
	for clust,i in zip(clustp,range(len(clustp))):
		plt.vlines(clust,0,np.amax(y),colors=color[i])

	if(x_r!='auto'):
		plt.xlim((x_r[0],x_r[1]))

	plt.ylabel(ylabel)
	plt.title(title)
	plt.xlabel(xlabel)
	plt.legend(loc="upper left",ncol=3)
	
	plt.tight_layout()
	plt.savefig(file.replace('.dat','_')+file2.rsplit('/',1).pop().replace('.dat','.png'),dpi=(200))
	
	if(interactive):
		plt.show()




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('file_name2', help='Name of dat file with cluster points')
	parser.add_argument('--int', default=0, help='Interactive')
	parser.add_argument('--x_r',dest='x_r',type=float, nargs=2, help='X-axis range, two numbers')
	#myplot("rmsd.dat")
	args = parser.parse_args()
	x_r='auto'
	if(args.x_r):
		x_r=args.x_r
	myplot(args.file_name,args.file_name2,x_r,args.int)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 