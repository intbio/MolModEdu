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
import re


def namdplot(file_in,file_out,step):
	"""
	Parses NAMD log file
	"""
	with open(file_in,'r') as f:
		for line in f:
        		if(re.match('ETITLE',line)):
        			print line
        			titles=line
        			break

	tlist=titles.split()[1:]
	print tlist
	p=re.compile('ENERGY')
	out=open(file_out,'w')
	out.write("NAMD energy terms from file %s \n" % file_in)
	out.write("parsed with python script \n")
	out.write("Time, ps\n")
	out.write("Energy \ pressure terms\n")
	out.write('\t'.join(tlist)+'\n')
	print 'Reading with step ', step
	with open(file_in,'r') as f:
		s=0
		for line in f:
        		if(p.match(line)):
        			s=s+1;
        			if(s==step):
        				s=0
        				terms=line.split()[1:]
        				out.write('%f\t'%(float(terms[0])*0.002)+'\t'.join(terms[1:])+'\n')
        			#print terms
        			

	out.close()






if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Extract data from NAMD log file')
	parser.add_argument('file_name', help='Name of dat file to analyze')
	parser.add_argument('-o', dest='output_f_name', help='Output file name')
	parser.add_argument('-s', dest='step', type=int, default=1, help='Load data every Nth step')
	args = parser.parse_args()
	namdplot(args.file_name,args.output_f_name,args.step)
	
	
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 