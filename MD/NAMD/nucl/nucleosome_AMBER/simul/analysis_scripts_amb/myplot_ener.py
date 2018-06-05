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
#import scipy


def myplot(file,first,width):
	"Plot dat file"

	
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['ps.useafm'] = True
	# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	# plt.rcParams['pdf.fonttype'] = 42
	font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 8}

	rc('font', **font)

	ox=[]
	#y=np.array([[1,1],[1,1]])
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
		
 		y=y.astype(float)
 		

 	plt.figure()
 	plt.title(title)
 	plt.gca()
 	plt.subplot(321)
 	#plt.yscale('log')
 	for i in [9]:
		plt.plot(ox[first:],y[first:,i], linewidth=width, label=legends[i+1])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
	#plt.xlim(xmin=1000)
	#plt.ylim(ymin=np.amin(y[100:,9]),ymax=np.amax(y[100:,9]))
	#print np.amin(y[100:,9])
	#plt.ylim(ymin=0,ymax=100)
	#plt.autoscale(enable=True,axis='y')
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.legend(loc="upper right",ncol=5)

	#let's get autocorrelation time
	data=y[first:,i]
	m=np.mean(data)
	if(len(data)<100):
		raise Exception("Data less than 100!!! quiting!")
	d1=map(lambda x: x-m ,data[:100])
	d2=d1[:50]
	acf=np.correlate(d1,d2,mode='valid')
	acf_max=acf[0]
	acf=map(lambda x: x/acf_max,acf)
	#Let's now find the right step, when acf falls below 1/2.72
	for t in range(0,49):
		if(acf[t]< (1/2.72)):
			ac_step=t;
			break;
	#print acf
	
	stdev=np.std(data)
	unc=stdev/((float(len(data))/float(ac_step))**0.5)
	print legends[i+1],' mean: ',m, 'Std.dev.: ',stdev , ' Uncert.:', unc, ' Cor.steps:', ac_step  


	plt.subplot(322)
	for i in [10]:
		plt.plot(ox[first:],y[first:,i], linewidth=width, label=legends[i+1])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.legend(loc="upper right",ncol=5)
	
#let's get autocorrelation time
	data=y[first:,i]
	m=np.mean(data)
	if(len(data)<100):
		raise Exception("Data less than 100!!! quiting!")
	d1=map(lambda x: x-m ,data[:100])
	d2=d1[:50]
	acf=np.correlate(d1,d2,mode='valid')
	acf_max=acf[0]
	acf=map(lambda x: x/acf_max,acf)
	#Let's now find the right step, when acf falls below 1/2.72
	for t in range(0,49):
		if(acf[t]< (1/2.72)):
			ac_step=t;
			break;
	#print acf
	
	stdev=np.std(data)
	unc=stdev/((float(len(data))/float(ac_step))**0.5)
	print legends[i+1],' mean: ',m, 'Std.dev.: ',stdev , ' Uncert.:', unc, ' Cor.steps:', ac_step  



	plt.subplot(323)
	for i in [15]:
		plt.plot(ox[first:],y[first:,i], linewidth=width, label=legends[i+1])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.legend(loc="upper right",ncol=5)
	
#let's get autocorrelation time
	data=y[first:,i]
	m=np.mean(data)
	if(len(data)<100):
		raise Exception("Data less than 100!!! quiting!")
	d1=map(lambda x: x-m ,data[:100])
	d2=d1[:50]
	acf=np.correlate(d1,d2,mode='valid')
	acf_max=acf[0]
	acf=map(lambda x: x/acf_max,acf)
	#Let's now find the right step, when acf falls below 1/2.72
	for t in range(0,49):
		if(acf[t]< (1/2.72)):
			ac_step=t;
			break;
	#print acf
	
	stdev=np.std(data)
	unc=stdev/((float(len(data))/float(ac_step))**0.5)
	print legends[i+1],' mean: ',m, 'Std.dev.: ',stdev , ' Uncert.:', unc, ' Cor.steps:', ac_step  


	plt.subplot(324)
	for i in [16]:
		plt.plot(ox[first:],y[first:,i], linewidth=width, label=legends[i+1])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.legend(loc="upper right",ncol=5)
	#let's get autocorrelation time
	data=y[first:,i]
	acf=np.correlate(data[:100],data[:100],mode='full')[99:]
	acf_max=acf[0]
	acf=map(lambda x: x/acf_max,acf)
	
#let's get autocorrelation time
	data=y[first:,i]
	m=np.mean(data)
	if(len(data)<100):
		raise Exception("Data less than 100!!! quiting!")
	d1=map(lambda x: x-m ,data[:100])
	d2=d1[:50]
	acf=np.correlate(d1,d2,mode='valid')
	acf_max=acf[0]
	acf=map(lambda x: x/acf_max,acf)
	#Let's now find the right step, when acf falls below 1/2.72
	for t in range(0,49):
		if(acf[t]< (1/2.72)):
			ac_step=t;
			break;
	#print acf
	
	stdev=np.std(data)
	unc=stdev/((float(len(data))/float(ac_step))**0.5)
	print legends[i+1],' mean: ',m, 'Std.dev.: ',stdev , ' Uncert.:', unc, ' Cor.steps:', ac_step  


	plt.subplot(325)
	for i in [4]:
		plt.plot(ox[first:],y[first:,i], linewidth=width, label=legends[i+1])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.legend(loc="upper right",ncol=5)
	
	
#let's get autocorrelation time
	data=y[first:,i]
	m=np.mean(data)
	if(len(data)<100):
		raise Exception("Data less than 100!!! quiting!")
	d1=map(lambda x: x-m ,data[:100])
	d2=d1[:50]
	acf=np.correlate(d1,d2,mode='valid')
	acf_max=acf[0]
	acf=map(lambda x: x/acf_max,acf)
	#Let's now find the right step, when acf falls below 1/2.72
	for t in range(0,49):
		if(acf[t]< (1/2.72)):
			ac_step=t;
			break;
	#print acf
	
	stdev=np.std(data)
	unc=stdev/((float(len(data))/float(ac_step))**0.5)
	print legends[i+1],' mean: ',m, 'Std.dev.: ',stdev , ' Uncert.:', unc, ' Cor.steps:', ac_step  



	plt.subplot(326)
	for i in [5]:
		plt.plot(ox[first:],y[first:,i], linewidth=width, label=legends[i+1])
	plt.ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.legend(loc="upper right",ncol=5)
	
#let's get autocorrelation time
	data=y[first:,i]
	m=np.mean(data)
	if(len(data)<100):
		raise Exception("Data less than 100!!! quiting!")
	d1=map(lambda x: x-m ,data[:100])
	d2=d1[:50]
	acf=np.correlate(d1,d2,mode='valid')
	acf_max=acf[0]
	acf=map(lambda x: x/acf_max,acf)
	#Let's now find the right step, when acf falls below 1/2.72
	for t in range(0,49):
		if(acf[t]< (1/2.72)):
			ac_step=t;
			break;
	#print acf
	
	stdev=np.std(data)
	unc=stdev/((float(len(data))/float(ac_step))**0.5)
	print legends[i+1],' mean: ',m, 'Std.dev.: ',stdev , ' Uncert.:', unc, ' Cor.steps:', ac_step  



	plt.savefig(file.replace('.dat','.png'),dpi=(600))





if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot data files')
	parser.add_argument('file_name', help='Name of dat file to plot')
	parser.add_argument('-first',dest="first", default=0, help='Name of dat file to plot')
	parser.add_argument('-w',dest="width", default=2, help='Name of dat file to plot')
	#myplot("rmsd.dat")
	args = parser.parse_args()
	myplot(args.file_name,int(args.first),args.width)
	


# dd = np.arange(2.5,3.5,0.1)    # the x locations for the groups

# 