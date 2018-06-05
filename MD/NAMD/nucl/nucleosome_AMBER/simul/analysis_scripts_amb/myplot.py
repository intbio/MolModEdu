#!/usr/bin/env python

"""
Program for plotting dat files with legends and names

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


def myplot(file,  max_y, x_r, interactive,leg_loc,bar,grid=False):
	"Plot dat file"

	
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['ps.useafm'] = True
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	plt.rcParams['pdf.fonttype'] = 42

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
				ox=np.append(ox,row[0])
				y=np.vstack((y,row[1:]))
 		y=y.astype(float) # to make it of float numpy type
 		ox=ox.astype(float)

 	plt.figure(figsize={12,6})
 	for i in range(len(legends)-1):
		if(bar):
			plt.bar(ox,y[:,i], width=0.5, label=legends[i+1],align="center")
		else:
			plt.plot(ox,y[:,i], linewidth=2, label=legends[i+1])

	if(grid):
		minorLocator   = MultipleLocator(0.5)
		# plt.gca().xaxis.set_minor_locator(minorLocator)
		plt.gca().yaxis.set_minor_locator(minorLocator)
		plt.gca().tick_params(axis='x',which='minor',bottom='on')
		plt.grid(which='minor')
		plt.grid()

	plt.ylabel(ylabel)
	plt.title(title)
	plt.xlabel(xlabel)
	plt.legend(loc="upper left",ncol=3)
	if(max_y!='auto'):
		plt.ylim((0,max_y[0]))
	if(x_r!='auto'):
		plt.xlim((x_r[0],x_r[1]))
	plt.legend(loc=leg_loc,ncol=4)
	plt.tight_layout()
	plt.savefig(file.replace('.dat','.png'),dpi=(200))
	
	if(interactive):
		plt.show()




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('--max_y',dest='max_y',type=float, nargs=1, help='Maximum y for plots')
	parser.add_argument('--x_r',dest='x_r',type=float, nargs=2, help='X-axis range, two numbers')
	parser.add_argument('--int', default=0, help='Interactive')
	parser.add_argument('--legend',dest='leg_loc',type=str, nargs=1, help='Legen location (def: upper left)')
	parser.add_argument('--bar',dest='bar',type=str, nargs=1, help='Display as bar')
	parser.add_argument('--grid',dest='grid',type=str, nargs=1, help='Display grid')

	#myplot("rmsd.dat")
	args = parser.parse_args()
	if(not args.leg_loc):
		leg_loc='upper left'
	else:
		leg_loc=args.leg_loc[0]
	if(not args.bar):
		bar=0
	else:
		bar=1
	x_r='auto'
	max_y='auto'
	if(args.max_y):
		max_y=args.max_y
	if(args.x_r):
		x_r=args.x_r

	if(args.grid):
		grid=True
	else:
		grid=False

	myplot(args.file_name,max_y,x_r,args.int,leg_loc,bar,grid)


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 