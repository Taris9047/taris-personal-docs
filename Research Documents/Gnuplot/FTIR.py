#!/usr/bin/env python
"""
    FTIR.py

	FTIR plot implementation for gnuplot.
	Also, can be extended for stochiometry extraction.
"""
import numpy as np
import os as os

debug_mode = 1

def isnumeric(value):
	return str(value).replace(".", "").replace("-", "")\
		   .replace("E","").replace("e","").isdigit()

# Smoothing function (just for 1st and 2nd order differentiation smoothing)
def smooth(x, window_len=10, window='hanning'):
	"""
	Smooth the data using a window with requested size.

	input:
		x: the input signal as array
		window_len: the dimension of window from 'flat', 'hanning','hamming','bartlett','blackman' flat window will produce a moving average smooting.

	output:
		the smoothed signal as array

	example:
		t=linspace(-2,2,0.1)
		x=sin(t)+randn(len(t))*.1
		y=smooth(x)

	see also:
		numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve, scipy.signal.lfilter

	"""
	x = np.array(x)
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."

	if x.size < window_len:
		raise ValueError, "Input Vector needs to be larger than window size."

	if window_len < 3:
		return x

	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window need to be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

	s = np.r_[2 * x[0] - x[window_len:1:-1], x, 2 * x[-1] - x[-1:-window_len:-1]]

	if window == 'flat': # moving average
		w = np.ones(window_len, 'd')
	else:
		w = eval('np.' + window + '(window_len)')

	y = np.convolve(w / w.sum(), s, mode='same')
	return y[window_len - 1:-window_len + 1]


def read_data(fn):
	data_file = open(fn, "r")
	
	# Reading the data file contents.
	wavenumber = []
	absorbance = []
	line_index = 0
	for line in data_file.readlines():
		# Excluding title text
		if str.find(line, "#") != -1:
			print "Excluding line:", line, "\n"
			line_index += 1
		else:
			wn, absor = line.split()
			# Checking whether the data points are valid or not.
			if isnumeric(wn) == False or isnumeric(absor) == False:
				print "Invalid data point at ", wn, " /cm."
			else:
				wavenumber.append(float(wn))
				absorbance.append(float(absor))
				line_index += 1

	data_file.close()
	print "\n"
	print "From file "+"\""+str(fn)+"\""+"...\n"
	print "Total data points extracted: "+str(line_index)+"\n"
	return wavenumber, absorbance

def search_value(key_value, given_array, fn_given_array):
	"""
	    search_value

		Finds the nearest value from function type array when
		the input value for an array defined function cannot
		be found easily.

		Input: key_value, x of function array, f(x) of function array
		
		Output: found value in x, f(x) of the found x
	"""
	# Checking data size
	if len(given_array) != len(fn_given_array):
		raise ValueError,\
			  "Array size mismatch! f(x) should show same size as its input."

	# Checking key_value range
	if key_value < min(given_array) or key_value > max(given_array):
		raise ValueError,\
			  "The value to evaluate must be within the range of the defined function!!"

	index_given_array = 0
	for valx_given in given_array:
		if valx_given > key_value:
			if (valx_given - key_value) >\
			   (key_value - given_array[index_given_array-1]):
				return given_array[index_given_array-1],\
					   fn_given_array[index_given_array-1]
			else:
				return given_array[index_given_array],\
					   fn_given_array[index_given_array]
		elif valx_given == key_value:
			return given_array[index_given_array],\
				   fn_given_array[index_given_array]
		else:
			index_given_array += 1

def generate_baseline(bslpt, wavenumber, absorbance):
	"""
	    genrate_baseline

		Generates the baseline for FTIR data.

		Input: x coordinates of baseline (must be obtained by inspection),
		       FTIR data in (wavenumber, absorbance) format.
			   
        Output: baseline data in (wavenumber, baseline_fit) format.
	"""
	# Obtaining the correct location of baseline fit.
	bslpt_raw = []
	bslval_raw = []
	for bslpt_element in bslpt:
		bslpt_raw.append(search_value(bslpt_element, wavenumber, absorbance)[0])
		bslval_raw.append(search_value(bslpt_element, wavenumber, absorbance)[1])		
	print "Baseline points found at:\n", bslpt_raw
	if debug_mode == 1:
		print "\n"
		print "bslpt_raw has, "+str(len(bslpt_raw))+" of elements."
		print "bslval_raw has, "+str(len(bslval_raw))+" of elements."
		print "\n"

	# Generating the linear fit values between baseline points.
	"""
	    To fit the baselinen points to f(x) = Ax+B format,
		we can derive A and B when f(x) is crossing (a1,b1), (a2,b2).
	
     	A = (b2-b1)/(a2-a1) and B = (a2b1 - a1b2)/(a2-a1)
	
	    In this case, we can assume (a2, b2) is the next data point.
	"""
	baseline = []
	for i in range(1, len(bslpt_raw)-1):
		A = (bslval_raw[i]-bslval_raw[i-1])/(bslpt_raw[i]-bslpt_raw[i-1])
		B = (bslpt_raw[i]*bslval_raw[i-1] - bslpt_raw[i-1]*bslval_raw[i])/\
			(bslpt_raw[i]-bslpt_raw[i-1])
		# Setting up wavenumber range for current section.
		wn_local = []
		for wn in wavenumber:
			if wn >= bslpt_raw[i-1] and wn < bslpt_raw[i]:
				wn_local.append(wn)
		# Update the baseline data points.
		for wnl in wn_local:
			baseline.append(A*wnl+B)

	
	#baseline.append(bslval_raw[len(bslval_raw)-1])

	return baseline

def main():
	# File name and path defined here
	data_file = "SiNxHadi.txt"
	baseline_file = "Baseline.txt"
	baseline_output = "Baseline_ary.txt"
	FTIR_output = "FTIR_sub.txt"
	print "Reading data file from, "+data_file
	print "Baseline definition found at "+baseline_file
	print "\n"

	wavenumber, absorbance = read_data(data_file)
	abosrbance = smooth(absorbance)
	bslindex, bslpt = read_data(baseline_file)

	# Generating baseline data plot
	baseline = generate_baseline(bslpt, wavenumber, absorbance)

	# Executing baseline subtraction.
	wn_trunc = []
	absorb_trunc = []
	for i in range(0,len(wavenumber)-1):
		if wavenumber[i] >= min(bslpt) and wavenumber[i] <= max(bslpt):
			wn_trunc.append(wavenumber[i])
			absorb_trunc.append(absorbance[i])

	FTIR_sub = []
	for i in range(0,len(baseline)):
		FTIR_sub.append(\
			absorb_trunc[i]-baseline[i])

	# Save Baseline and Subtracted data.
	bsl_out = open(baseline_output, "w")
	bsl_out_text = ["Generated by FTIR.py. Contains baseline plot data.\n"]
	for i in range(0,len(baseline)-1):
		bsl_out_text.append(str(wn_trunc[i])+"\t"+str(baseline[i])+"\n")
	bsl_out.writelines(bsl_out_text)
	bsl_out.close()
	FTIR_out = open(FTIR_output, "w")
	FTIR_out_text = ["Generated by FTIR.py. Contains subtracted FTIR data.\n"]
	for i in range(0,len(baseline)-1):
		FTIR_out_text.append(str(wn_trunc[i])+"\t"+str(FTIR_sub[i])+"\n")
	FTIR_out.writelines(FTIR_out_text)
	FTIR_out.close()
	
	return 0

# Running the Main Function
if __name__ == "__main__":
	main()
