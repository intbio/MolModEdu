#!/usr/bin/env python

"""
Program for parsing 3DNA ref_frame files for DNA broken line analysis and making a single file for the whole trajectory.
The ideas are to use Pandas data frame and R
Computes angles on the fly.

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
#import matplotlib.pyplot as plt
#from matplotlib import rc
import argparse
import csv
import os
import glob
import pandas as pd
import re

def main():
	"Parse 3DNA ref_frame files"

	#Let's get a list of files
	flist=glob.glob("../analysis_data/dna_params_bl/step*.dat")
	panel_dict={}
	for file in flist:
		print "Processing ", file
		#let's get the data from file
		# let's strip white spaces at first and replace with tabs
		tmpfile=file.replace('.dat','.tmp')
		with open(file,'r') as f:
			with open(tmpfile,'w') as tf:
				f.next()
				f.next()
				for line in f:
					tf.write('\t'.join(line.partition('#')[0].split())+'\n')
					next(f,'')
					next(f,'')
					next(f,'')
					next(f,'')


		#now we will read data frames, merege them to Panel, and output with pickle
		ind=re.search('\d+',file).group(0)
		df=pd.read_csv(tmpfile, sep=u'\t',skiprows=0,header=None)
		df.columns=['x','y','z']


		panel_dict[ind]=df


	p=pd.Panel(panel_dict)
#workaround to set types - does not help, run it in the script after data loading
# df=p['1']
# print p['1'].dtypes
# for k in p.items:
# 	for n in df.columns[1:]:
# 		p[k][n]=p[k][n].astype(float)
# 	p[k][df.columns[0]]=p[k][df.columns[0]].astype(str) # seems it will still be of object type
# print p['1'].dtypes
	p.to_pickle('../analysis_data/dna_bl_pandas_panel.dat')

# Let's recombine all into a single data frame with number of nucleotide pair being a factor
# and adding time variable
	bdf = pd.DataFrame()
	for key in panel_dict:
		value=panel_dict[key]
		value['Time']=key
		bdf=bdf.append(value)
	bdf.to_csv('../analysis_data/dna_bl_big_df.csv')

#Let's save also a dataframe for crystal structure
####
####

	file="../analysis_data/dna_params_bl/cryst.dat"
	print "Processing ", file
	#let's get the data from file
	# let's strip white spaces at first and replace with tabs
	tmpfile=file.replace('.dat','.tmp')
	with open(file,'r') as f:
		with open(tmpfile,'w') as tf:
			f.next()
			f.next()
			for line in f:
				tf.write('\t'.join(line.partition('#')[0].split())+'\n')
				next(f,'')
				next(f,'')
				next(f,'')
				next(f,'')


	df=pd.read_csv(tmpfile, sep=u'\t',skiprows=0,header=None)
	# df.index=range(147)
	df.columns=['x','y','z']
	df.to_csv('../analysis_data/dna_bl_cryst.csv')


	
	# if(interactive):
	# 	plt.show()




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('--int', default=0, help='Interactive')
	main()
	# args = parser.parse_args()
	# myplot(args.file_name,args.int)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 