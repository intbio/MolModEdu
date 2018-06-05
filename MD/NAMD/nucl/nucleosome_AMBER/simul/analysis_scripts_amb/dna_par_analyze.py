#!/usr/bin/env python

"""
Using pandas to analyze our DNA params data set

Copyright 2013 (c) Alexey Shaytan

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
import csv
import os
import glob
import pandas as pd
import re

def main():
	"Analyze DNA params"

	p=pd.read_pickle('../analysis_data/dna_param_pandas_panel.dat')
	for i in p.items:
		p[i]=p[i].convert_objects()
	a=p.minor_xs('Roll')
	a=a.transpose()
	b=a.unstack()
	b.hist(bins=1000)
	plt.show()
	plt.savefig('../analysis_data/test.png',dpi=(200))




if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('--int', default=0, help='Interactive')
	main()
	# args = parser.parse_args()
	# myplot(args.file_name,args.int)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 