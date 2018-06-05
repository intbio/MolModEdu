#!/usr/bin/env python

"""
Program for parsing Curves+ files and making a single file for the whole trajectory.
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
	"Parse Curve+ files"

	#Let's get a list of files
	flist=glob.glob("../analysis_data/dna_params_cur/step*.lis")
	panel_bb_dict={}
	panel_gr_dict={}

	for file in flist:
		print "Processing ", file
		#What we need to do here is to parse the Curves+ file
		#and get info from it
		# Right now we are only interested in DNA groove data and backbone data
		# with backbone we will be able to study BI-BII conformation

		# Probably later we can add parsing the parameters of base pair steps.
		# Also there are valies of H-Twi and H-Rise - am not sure how to interpret them
		# and whether we can get overall twist from them
		#now we will read data frames, merege them to Panel, and output with pickle
		ind=re.search('\d+',file).group(0)

		tmpfile=file.replace('.lis','.tmpb')
		tmpfile2=file.replace('.lis','.tmpg')

		#Let's parse the backbone
		with open(file,'r') as f: 
			with open(tmpfile,'w') as t:
				with open(tmpfile2,'w') as t2:
					#loop till we find
					#  (D) Backbone Parameters
					for line in f:
						if(re.search("\(D\)",line)): break
					# Strand 1     Alpha  Beta   Gamma  Delta  Epsil  Zeta   Chi    Phase  Ampli  Puckr  
					for line in f:
						if(re.search("Strand",line)):
							t.write('Strand\t'+'Resid\t'+'\t'.join(line.split()[2:])+'\n')
							f.next() # skip one blank line
							break
					for line in f:
						if(re.search("\d+\)\s+[ATGC]",line)):
							t.write('1'+'\t'+'\t'.join(line.split()[2:]).replace('----','')+'\n')
						else: break
					for line in f:
						if(re.search("Strand",line)):
							f.next() # skip one blank line
							break
					for line in f:
						if(re.search("\d+\)\s+[ATGC]",line)):
							t.write('2'+'\t'+'\t'.join(line.split()[2:]).replace('----','')+'\n')
						else: break
					for line in f:
						if(re.search("\(E\)",line)): break
					for line in f:
						if(re.search("Level",line)):
							t2.write('\t'.join(line.split())+'\n')
							f.next() # skip one blank line
							break
					for line in f:
						if(re.search("\s+\d+",line)):
							t2.write(line[0:8].strip()+'\t'+line[16:22].strip()+'\t'+line[23:30].strip()+'\t'+line[31:38].strip()+'\t'+line[39:46].strip()+'\n')
						else: break
					


					#print line	

		df_bb=pd.read_csv(tmpfile, sep=u'\t')
		panel_bb_dict[ind]=df_bb

		df_gr=pd.read_csv(tmpfile2, sep=u'\t')
		panel_gr_dict[ind]=df_gr

	p_bb=pd.Panel(panel_bb_dict)
	p_bb.to_pickle('../analysis_data/dna_bb_param_pandas_panel.dat')


	p_gr=pd.Panel(panel_gr_dict)
	p_gr.to_pickle('../analysis_data/dna_gr_param_pandas_panel.dat')


# Let's recombine all into a single data frame with number of nucleotide pair being a factor
# and adding time variable
	bdf = pd.DataFrame()
	for key in panel_bb_dict:
		value=panel_bb_dict[key]
		value['Time']=key
		bdf=bdf.append(value)
	bdf.to_csv('../analysis_data/dna_bb_param_big_df.csv')

	bdf = pd.DataFrame()
	for key in panel_gr_dict:
		value=panel_gr_dict[key]
		value['Time']=key
		bdf=bdf.append(value)
	bdf.to_csv('../analysis_data/dna_gr_param_big_df.csv')


################
####
#### Crystal
####

	file="../analysis_data/dna_params_cur/cryst.lis"
	tmpfile=file.replace('.lis','.tmpb')
	tmpfile2=file.replace('.lis','.tmpg')

	#Let's parse the backbone
	with open(file,'r') as f: 
		with open(tmpfile,'w') as t:
			with open(tmpfile2,'w') as t2:
				#loop till we find
				#  (D) Backbone Parameters
				for line in f:
					if(re.search("\(D\)",line)): break
				# Strand 1     Alpha  Beta   Gamma  Delta  Epsil  Zeta   Chi    Phase  Ampli  Puckr  
				for line in f:
					if(re.search("Strand",line)):
						t.write('Strand\t'+'Resid\t'+'\t'.join(line.split()[2:])+'\n')
						f.next() # skip one blank line
						break
				for line in f:
					if(re.search("\d+\)\s+[ATGC]",line)):
						t.write('1'+'\t'+'\t'.join(line.split()[2:]).replace('----','')+'\n')
					else: break
				for line in f:
					if(re.search("Strand",line)):
						f.next() # skip one blank line
						break
				for line in f:
					if(re.search("\d+\)\s+[ATGC]",line)):
						t.write('2'+'\t'+'\t'.join(line.split()[2:]).replace('----','')+'\n')
					else: break
				for line in f:
					if(re.search("\(E\)",line)): break
				for line in f:
					if(re.search("Level",line)):
						t2.write('\t'.join(line.split())+'\n')
						f.next() # skip one blank line
						break
				for line in f:
					if(re.search("\s+\d+",line)):
						t2.write(line[0:8].strip()+'\t'+line[16:22].strip()+'\t'+line[23:30].strip()+'\t'+line[31:38].strip()+'\t'+line[39:46].strip()+'\n')
					else: break
					

	df_bb=pd.read_csv(tmpfile, sep=u'\t')

	df_gr=pd.read_csv(tmpfile2, sep=u'\t')

	df_bb.to_csv('../analysis_data/dna_bb_param_cryst.csv')

	df_gr.to_csv('../analysis_data/dna_gr_param_cryst.csv')




# #Let's save also a dataframe for crystal structure
####
# ####

# 	file="../analysis_data/dna_params/cryst.dat"
# 	print "Processing ", file
# 	#let's get the data from file
# 	# let's strip white spaces at first and replace with tabs
# 	tmpfile=file.replace('.dat','.tmp')
# 	with open(file,'r') as f:
# 		with open(tmpfile,'w') as tf:
# 			for line in f:
# 				tf.write('\t'.join(line.split())+'\n')

# 	bp_list=list()
# 	inpfile=file.replace('.dat','.inp')
# 	with open(inpfile,'r') as inf:
# 		for line in inf:
# 			if(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line)):
# 				bp_list.append(int(re.search('\.\.\.\.>C:\.*(-?\d+)_:',line).group(1))+73)
# 	print bp_list

# 	df=pd.read_csv(tmpfile, sep=u'\t',skiprows=2)
# 	#The first column is of the type A-T or A+T, the latter indicating that we do not have normal basepairing
# 	# We need to converti it to 0 (not normal), 1 (normal base pairing)
# 	df=df.rename(columns={'#':'BP'})
# 	for i in range(len(df)):
# 		if(i in bp_list):
# 			df.iloc[i,0]=1
# 		else:
# 			df.iloc[i,0]=0
# 	df.iloc[0,7:]=np.nan
# 	df.to_csv('../analysis_data/dna_param_cryst.csv')


	
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