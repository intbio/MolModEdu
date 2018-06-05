#!/usr/bin/env python

"""
Program for parsing 3DNA files and making a single file for the whole trajectory.
The ideas are to use Pandas data frame and R

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
	"Parse 3DNA files"

	#Let's get a list of files
	flist=glob.glob("../analysis_data/dna_params/step*.dat")
	panel_dict={}
	for file in flist:
		print "Processing ", file
		#let's get the data from file
		# let's strip white spaces at first and replace with tabs
		tmpfile=file.replace('.dat','.tmp')
		with open(file,'r') as f:
			with open(tmpfile,'w') as tf:
				for line in f:
					tf.write('\t'.join(line.split())+'\n')
#Now we need to see what BP are not determined
#!!!! This depends on the seq. numbering adjust it !!!!!!
#   1  294  0 #    1 | ....>C:.-73_:[ADE]A-----T[THY]:..73_:C<....  0.74  0.34  5.88  8.43 -1.29
		bp_list=list()
		inpfile=file.replace('.dat','.inp')
		with open(inpfile,'r') as inf:
			for line in inf:
				if(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line)):
					bp_list.append(int(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line).group(1))+73)
		print bp_list
		#now we will read data frames, merege them to Panel, and output with pickle
		ind=re.search('\d+',file).group(0)
		df=pd.read_csv(tmpfile, sep=u'\t',skiprows=2)
		#The first column is of the type A-T or A+T, the latter indicating that we do not have normal basepairing
		# We need to converti it to 0 (not normal), 1 (normal base pairing)
		df=df.rename(columns={'#':'BP'})
		for i in range(len(df)):
			if(i in bp_list):
				df.iloc[i,0]=1
				#print  i
			else:
				df.iloc[i,0]=0
# we need to set the first params of base-steps to NaN
		df.iloc[0,7:]=np.nan

#Now let's add the torsion parameters to data frame.
# We need to slice stepN.tor file into two files for angles and sugar puckers, which we will read as csv
		torfile=file.replace('.dat','.tor')
		angfile=file.replace('.dat','.bba')
		puckfile=file.replace('.dat','.pck')
		with open(torfile,'r') as f:
			f.readline()
			with open(angfile,'w') as tf:
				for line in f:
					if('****' not in line):
						plist=line.split()
						if(len(plist)>5):
							if((len(plist)!=12) and (plist[0]!='base')):
								plist.insert(3,'no')
						tf.write('\t'.join(plist)+'\n')
					else:
						break
			for line in f:
				if('*****' in line): break
			with open(puckfile,'w') as tf:
				for line in f:
					tf.write('\t'.join(line.split())+'\n')

		dfa=pd.read_csv(angfile, sep=u'\t',skiprows=19)
		dfp=pd.read_csv(puckfile, sep=u'\t',skiprows=18)
		#Now the idea is to slice this data frames into strands
		dfa1=dfa.iloc[0:len(dfa)/2,1:]
		dfa2=dfa.iloc[len(dfa)/2:,1:]
		n=dfa1.columns
		n1=[]
		n2=[]
		for i in n:
			n1.append(i+'_1')
			n2.append(i+'_2')
		dfa1.columns=n1
		dfa2.columns=n2
		dfa1.index=range(len(dfa1))
		dfa2.index=range(len(dfa2)-1,-1,-1)
		dfa_new=pd.concat([dfa1,dfa2],axis=1)

		dfp1=dfp.iloc[0:len(dfp)/2,1:]
		dfp2=dfp.iloc[len(dfp)/2:,1:]
		n=dfp1.columns
		n1=[]
		n2=[]
		for i in n:
			n1.append(i+'_1')
			n2.append(i+'_2')
		dfp1.columns=n1
		dfp2.columns=n2
		dfp1.index=range(len(dfp1))
		dfp2.index=range(len(dfp2)-1,-1,-1)
		dfp_new=pd.concat([dfp1,dfp2],axis=1)

		df_new=pd.concat([df,dfa_new,dfp_new],axis=1)
		panel_dict[ind]=df_new


	p=pd.Panel(panel_dict)
#workaround to set types - does not help, run it in the script after data loading
# df=p['1']
# print p['1'].dtypes
# for k in p.items:
# 	for n in df.columns[1:]:
# 		p[k][n]=p[k][n].astype(float)
# 	p[k][df.columns[0]]=p[k][df.columns[0]].astype(str) # seems it will still be of object type
# print p['1'].dtypes
	p.to_pickle('../analysis_data/dna_param_pandas_panel.dat')

# Let's recombine all into a single data frame with number of nucleotide pair being a factor
# and adding time variable
	bdf = pd.DataFrame()
	for key in panel_dict:
		value=panel_dict[key]
		value['Time']=key
		bdf=bdf.append(value)
	bdf.to_csv('../analysis_data/dna_param_big_df.csv')

#Let's save also a dataframe for crystal structure
####
####

	file="../analysis_data/dna_params/cryst.dat"
	print "Processing ", file
	#let's get the data from file
	# let's strip white spaces at first and replace with tabs
	tmpfile=file.replace('.dat','.tmp')
	with open(file,'r') as f:
		with open(tmpfile,'w') as tf:
			for line in f:
				tf.write('\t'.join(line.split())+'\n')

	bp_list=list()
	inpfile=file.replace('.dat','.inp')
	with open(inpfile,'r') as inf:
		for line in inf:
			if(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line)):
				bp_list.append(int(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line).group(1))+73)
	print bp_list

	df=pd.read_csv(tmpfile, sep=u'\t',skiprows=2)
	#The first column is of the type A-T or A+T, the latter indicating that we do not have normal basepairing
	# We need to converti it to 0 (not normal), 1 (normal base pairing)
	df=df.rename(columns={'#':'BP'})
	for i in range(len(df)):
		if(i in bp_list):
			df.iloc[i,0]=1
		else:
			df.iloc[i,0]=0
	df.iloc[0,7:]=np.nan


#Now let's add the torsion parameters to data frame.
# We need to slice stepN.tor file into two files for angles and sugar puckers, which we will read as csv
	torfile='../analysis_data/dna_params/cryst.tor'
	angfile='../analysis_data/dna_params/cryst.bba'
	puckfile='../analysis_data/dna_params/cryst.pck'
	with open(torfile,'r') as f:
		f.readline()
		with open(angfile,'w') as tf:
			for line in f:
				if('****' not in line):
						plist=line.split()
						if(len(plist)>5):
							if((len(plist)!=12) and (plist[0]!='base')):
								plist.insert(3,'no')
						tf.write('\t'.join(plist)+'\n')
				else:
					break
		for line in f:
			if('*****' in line): break
		with open(puckfile,'w') as tf:
			for line in f:
				tf.write('\t'.join(line.split())+'\n')

	dfa=pd.read_csv(angfile, sep=u'\t',skiprows=19)
	dfp=pd.read_csv(puckfile, sep=u'\t',skiprows=18)
	#Now the idea is to slice this data frames into strands
	dfa1=dfa.iloc[0:len(dfa)/2,1:]
	dfa2=dfa.iloc[len(dfa)/2:,1:]
	n=dfa1.columns
	n1=[]
	n2=[]
	for i in n:
		n1.append(i+'_1')
		n2.append(i+'_2')
	dfa1.columns=n1
	dfa2.columns=n2
	dfa1.index=range(len(dfa1))
	dfa2.index=range(len(dfa2)-1,-1,-1)
	dfa_new=pd.concat([dfa1,dfa2],axis=1)

	dfp1=dfp.iloc[0:len(dfp)/2,1:]
	dfp2=dfp.iloc[len(dfp)/2:,1:]
	n=dfp1.columns
	n1=[]
	n2=[]
	for i in n:
		n1.append(i+'_1')
		n2.append(i+'_2')
	dfp1.columns=n1
	dfp2.columns=n2
	dfp1.index=range(len(dfp1))
	dfp2.index=range(len(dfp2)-1,-1,-1)
	dfp_new=pd.concat([dfp1,dfp2],axis=1)

	df_new=pd.concat([df,dfa_new,dfp_new],axis=1)


	df_new.to_csv('../analysis_data/dna_param_cryst.csv')


	
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