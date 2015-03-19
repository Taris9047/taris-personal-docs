'''
Created on Apr. 5, 2011

@author: Taylor Shin

polynomial fitting for TFT transfer curve at saturation region.

'''

from Tkinter import *
import tkFileDialog as tkFD
import tkSimpleDialog as tkSD
import os, sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from scipy.interpolate import interp1d

def poly_fit(p0, x):
	#return p0[0]*x**2 + p0[1]*x + p0[2]
	return p0[0]*(x+p0[1])**2 + p0[2]
	#return p0[0]*x + p0[1]
	#return p0[0]*x**5 + p0[1]*x**4 + p0[2]*x**3 + p0[3]*x**2 + p0[4]*x + p0[5]

def poly_residual(p, y, x):
	"""
	Defining residual for Transfer curve fitting
	"""	
	err = y - poly_fit(p ,x)
	return err
	
# returns whether the input string is numeric format or not.
def isnumeric(value):
	return str(value).replace("e","")\
	.replace("E","").replace(".", "")\
	.replace("-", "").replace("+","").isdigit()


# Reads data from spectroscopy output to python workspace.
def trans_readfile(fn):
	"""
	Reads optical spectra data.

	Input: 
		Filename (including path)

	Output: 
		wl, tm

	"""
	originalFile = open(fn, "r")
	# adjusting filename for further use
	fn = fn.replace(".txt", "")
	voltage = []
	current = [] 
	lineIndex = 0 

	for line in originalFile.readlines():
		# Excluding text field
		if str.find(line, "\"") != -1:
			lineIndex += 1
		else:
			line = line.replace(" ", "\t")
			line = line.replace(",", "\t")
			if line.find("\r") != -1:
				line = line.replace("\r","")
			if line.find("\n") != -1:
				line = line.replace("\n","")
			volt, curr = line.split()
			#print line.split("\t")
			if isnumeric(volt) == False or isnumeric(curr) == False:
				print "Invalid point at ", volt, " V, (empty data)"
			else:
				voltage.append(float(volt))
				current.append(float(curr))
				lineIndex += 1

	print "Total data points: " + str(lineIndex-1) + "\n"

	return np.array(voltage), np.array(current)


# Main Function
def main():
	# Reading the file
	'''
	root = Tk()
	root.withdraw()
	filename = \
		tkFD.askopenfilename( \
			defaultextension="txt", initialdir="./", \
			title="Select an Transfer Characteristic Data.")
	if os.path.isfile(filename) != True:
		raise ValueError, "Invalid file or cannot find designated file!"
	else:
		print "Reading the spectrum data from ..."
		print filename
		# Reading data from file.
		volt, curr = trans_readfile(filename)
	'''
	filename = './SiTFTTRANSDK.txt'
	print "Reading the spectrum data from ..."
	print filename
	# Reading data from file.
	volt, curr = trans_readfile(filename)
	volt = np.array(volt)
	curr = np.array(curr)
	
	# Refining data steps for original data.
	interp_volt = np.linspace(0,20,21)
	interp_curr_function = interp1d(volt, curr)
	curr = interp_curr_function(interp_volt)
	volt = interp_volt

	# Now, fitting the Transfer curve with 2nd order Polynomial function
	p0 = np.array([.404, -5, 1.3])
	#p0 = np.array([1, -1, 1])
	#p0 = np.array([0.404,-5])
	#p0 = np.array([1,1,1,1,1,0])
	
	# selecting valid range
	volt_fit_range = []
	curr_fit_range = []
	index = 0
	threshold = 5
	for volt_val in volt:
		if volt_val >= threshold:
			volt_fit_range.append(volt[index])
			curr_fit_range.append(curr[index])
		index += 1
	# converting both arrays to np.array
	volt_fit_range = np.array(volt_fit_range)
	curr_fit_range = np.array(curr_fit_range)
	
	# Fitting with normalized parameters
	norm_factor = 1e6
	curr_fit_range_norm = curr_fit_range * norm_factor
	
	para_sq, success = \
		leastsq(poly_residual, p0, \
		 	args=(volt_fit_range, curr_fit_range_norm))
	print para_sq, success
	print "Ids = ", para_sq[0], "(Vgs + ", para_sq[1], ")^2"
	
	# Regenerating fitted function into array form
	fitted_curr = []
	volt_fit_graph =  \
		np.linspace(volt_fit_range[0], volt_fit_range[np.size(volt_fit_range)-1])
	for volt_val in volt_fit_graph:
		fitted_curr.append(poly_fit(para_sq, volt_val)/norm_factor)
	fitted_curr = np.array(fitted_curr)
	#print fitted_curr
	
	initial_fit = []
	for volt_val in volt_fit_graph:
		initial_fit.append(poly_fit(p0, volt_val)/norm_factor)
	initial_fit = np.array(initial_fit)
	#print initial_fit

	
	# Checking the fitting compared to original data.
	plt.figure()
	#plt.yscale('log')
	plt.plot(volt, curr*norm_factor, 'r-')
	plt.plot(volt_fit_graph, fitted_curr*norm_factor, 'go')
	plt.plot(volt_fit_graph, initial_fit*norm_factor, 'bs-')
	plt.legend([r'Original Data', r'Fitted Curve', r'Initial Fit'], loc='best')
	plt.title(r'Transfer Curve Fitting Check.', fontsize=16)
	plt.xlabel(r'Gate Voltage (V)', fontsize=16)
	plt.ylabel(r'Drain Current ($\mu$A/$\mu$m)', fontsize=16)
	ax = plt.gca()
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	plt.show()
	
	
	
	
	
# running the program	
if __name__ == "__main__":
	main()