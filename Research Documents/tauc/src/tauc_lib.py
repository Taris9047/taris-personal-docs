"""
	tauc_lib.py
	
	Library for tauc.py
	The main purpose is to contain various functions in a certain file 
	to avoid complex code in the main function.
"""
import numpy as np
import os as os

# Preprocessor for Debug
# Set, debug_mode = 1 for scientific debug.
debug_mode = 0

# isnumeric
# returns whether the input string is numeric format or not.
def isnumeric(value):
	return str(value).replace(".", "").replace("-", "").isdigit()

# tauc_readfile
# Reads data from spectroscopy output to python workspace.
def tauc_readfile(fn):
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
	wavelength = []
	transmittance = [] 
	lineIndex = 0 
	
	for line in originalFile.readlines():
		# Excluding text field
		if str.find(line, "\"") != -1:
			lineIndex += 1
		else:
			line = line.replace(",", "\t")
			if line.find("\r") != -1:
				line = line.replace("\r","")
			if line.find("\n") != -1:
				line = line.replace("\n","")
			wl, tm = line.split("\t")
			#print line.split("\t")
			if isnumeric(wl) == False or isnumeric(tm) == False:
				print "Invalid point at ", wl, " nm, (empty data)"
			else:
				wavelength.append(float(wl))
				transmittance.append(float(tm))
				lineIndex += 1

	print "Total data points: " + str(lineIndex-1) + "\n"
	
	return np.array(wavelength), np.array(transmittance)

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
	if np.array(x).ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."

	if np.array(x).size < window_len:
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

# findzero
# Find out zero position for array.
def findzero(x, fn, threshold=0.01):
	"""
	Finds out zero crossings from an array defined function.
	
	Input:
		x: x axis of array defined function.
		fn: y value of array defined function.
		threhsold: tolerance of zero value. if y < threshold, then it will be
			considered as zero.
	
	Output:
		zeroc: array of zero crossing points.
	"""
	if np.array(x).ndim != 1 and np.array(fn).ndim != 1:
		raise ValueError, "x,fn must be linear vector."
	if np.array(x).size != np.array(fn).size:
		raise ValueError, "x,fn must have same sizes."
	
	zeroc = []
	x_index = 0
	#zeros_count = 0
	for item in fn:
		if item >= -threshold and item <= threshold:
			# for debug
			#print "Zero cross at " + str(x[x_index]) + " found!"
			#zeros_count = zeros_count + 1
			zeroc.append(x[x_index])
		else:
			# Do nothing
			x_index = x_index + 0
			# for debug
			#print "Zero cross not found at " + str(x[x_index])
		x_index = x_index + 1
		
	#for debug
	#print zeroc
	#print "Total # of zero crossing: " + str(zeros_count)
	return zeroc

# readfunction
# Reads out a value from array defined function using linear interpolation
# method		
def readfunction(x, fnx, fn): 
	"""
	Input:
		x: X value input.
		fnx: X component of function.
		fn: function itself.
	
	Output:
		fny: Y value at x in X of function.
		
	"""
	if np.array(fnx).ndim != 1 or np.array(fn).ndim != 1:
		raise ValueError, "Function should be defined as array."
	if np.array(fnx).size != np.array(fn).size:
		raise ValueError, "X axis data and Y axis data should match in dimension"
	
	index = 0
	for lower in fnx:
		higher = fnx[index + 1]
		lowfn = fn[index]
		highfn = fn[index + 1]
		if x > higher:
			index = index + 1
		else:
			fny = (highfn - lowfn) / (higher - lower) * (x - lower) + lowfn
			index = index + 1
			break
	
	return fny

# readfunction_a (Array version)
# Reads out a value from array defined function using linear interpolation
# method.		
def readfunction_a(x, fnx, fn): 
	"""
	Input:
		x: X value array input.
		fnx: X component of function.
		fn: function itself.
	
	Output:
		fny: Y value array at x in X of function.
		
	"""
	if np.array(fnx).ndim != 1 or np.array(fn).ndim != 1:
		raise ValueError, "Function should be defined as array."
	if np.array(fnx).size != np.array(fn).size:
		raise ValueError, "X axis data and Y axis data should match in dimension"
	
	fny = []
	for value in x:
		index = 0
		for lower in fnx:
			higher = fnx[index + 1]
			lowfn = fn[index]
			highfn = fn[index + 1]
			if value > higher:
				index = index + 1
			else:
				fny.append((highfn - lowfn) / (higher - lower) * (value - lower) + lowfn)
				index = index + 1
				break
	
	#for Debug
	#print "Lower limit: " + str(lower)
	#print "Higher limit: " + str(higher)
	#print "Value: " + str(fny) 
	
	return fny

# findpeak
# Finds out peak value using 1st and 2nd order derivatives.
def findpeak(x, fn, threshold=0.01):
	"""
	Input:
		x: x portion of array defined function.
		fn: the function value itself.

	Output:
		maxx: array of x value where peak values reside at.
	""" 
	if np.array(x).ndim != 1 or np.array(fn).ndim != 1:
		raise ValueError, "Function portions must be Vector form."

	if np.array(x).size != np.array(fn).size:
		raise ValueError, "X, Y need to have same size as array."

	fnd = np.diff(fn, 1)
	fnds = smooth(fnd, 30)
	fndd = np.diff(fnds, 1)
	fndds = smooth(fndd, 100)

	zerosd = findzero(np.linspace(min(x), max(x), len(x) - 1), fnds, threshold)

	maxx = []
	for zero_pos in zerosd:
		if readfunction(zero_pos, np.linspace(min(x), max(x), len(x) - 2), fndds) > 0:
			pass
		else:
			maxx.append(zero_pos)

	return maxx

# findpeak_revised
def findpeak_rev(x, fn, delta=3):
	"""
	Input:
		x: x portion of array defined function.
		fn: the function value itself.
		delta: threshold of determining peak or bottom. 
		       Default is 3 for this application.

	Output:
		maxx: array of x value where peak values reside at.
		maxxfn: array of peak values. (not implemented for this app.)
		minn: array of x value where bottom values reside at.
		minnfn: array of bottom values. (not implemented for this app.)
	""" 
	if np.array(x).ndim != 1 or np.array(fn).ndim != 1:
		raise ValueError, "Function portions must be Vector form."

	if np.array(x).size != np.array(fn).size:
		raise ValueError, "X, Y need to have same size as array."
		
	if delta <= 0:
		raise ValueError, "delta must be positive value!"	
		
	# Getting peaks
	findpeak = 1
	max_local = -float('inf')
	min_local = float('inf')
	index = 0
	maxx = []
	maxxfn = []
	minn = []
	minnfn = []
	for fnval in fn:
		# Updating local max/min points
		if fnval > max_local:
			max_local = fnval
			max_local_pos = x[index]
		
		if fnval < min_local:
			min_local = fnval
			min_local_pos = x[index]
		
		# Updating peak/bottom position depending on delta value
		if findpeak == 1:
			if fnval < max_local - delta:
				maxx.append(max_local_pos)
				maxxfn.append(max_local)
				min_local = fnval
				min_local_pos = x[index]
				findpeak = 0
		else:
			if fnval > min_local + delta:
				minn.append(min_local_pos)
				minnfn.append(min_local)
				max_local = fnval
				max_local_pos = x[index]
				findpeak = 1
			
		index += 1
		
	if debug_mode == 1:
		print "Maxx, minn", maxx, minn
		
	return maxx, maxxfn, minn, minnfn
	
# findbottom
# Finds out peak value using 1st and 2nd order derivatives.
def findbottom(x, fn, threshold=0.01):
	"""
	Input:
		x: x portion of array defined function.
		fn: the function value itself.
	Output:
		maxx: array of x value where peak values reside at.
	""" 
	if np.array(x).ndim != 1 or np.array(fn).ndim != 1:
		raise ValueError, "Function portions must be Vector form."
	
	if np.array(x).size != np.array(fn).size:
		raise ValueError, "X, Y need to have same size as array."
	
	fnd = np.diff(fn, 1)
	fnds = smooth(fnd, 20)
	fndd = np.diff(fnds, 1)
	fndds = smooth(fndd, 30)
	
	zerosd = findzero(np.linspace(min(x), max(x), len(x) - 1), fnds, threshold)
	
	minn = []
	for zero_pos in zerosd:
		if readfunction(zero_pos, np.linspace(min(x), max(x), len(x) - 2), fndds) < 0:
			pass
		else:
			minn.append(zero_pos)
	
	return minn


# Thickness input parser
def thick_input(thick):
	"""
	Converts the thickness input to number.
	
	Input:
		thick (string type w or w/o unit)
		
	Output:
		floating point type of thick(nm unit), unit of the input.
	"""
	
	if type(thick) != str:
		raise ValueError, "Input must be string formatted integer or float!"
	
	# Stripping out space
	thick = thick.replace(' ','')
	
	# If the input is negative?
	if thick[0] == '-':
		raise ValueError, "Negative Thickness cannot be realistic!"

	# Unit detection
	unit = thick[len(thick)-3:].replace(' ','')
	if unit.find('nm') != -1:
		umod = 1 # Default return unit: 1e-9, nm
		thick = thick.replace('nm','')
		unit = 'nm'
	elif unit.find('um') != -1:
		umod = 1e-3
		thick = thick.replace('um','')
		unit = 'um'
	elif unit.find('mm') != -1:
		umod = 1e-6
		thick = thick.replace('mm','')
		unit = 'mm'
	elif unit.find('cm') != -1:
		umod = 1e-7
		thick = thick.replace('cm','')
		unit = 'cm'
	elif unit.find('A') != -1:
		umod = .1
		unit = 'A'
		thick = thick.replace('A','')
	elif unit.isdigit() == True:
		# When the input has no unit, consider it as nm
		umod = 1
		unit = 'nm'
	else:
		raise ValueError, "The input must have valid unit(A,nm,um,mm,cm) or has no unit!"
		
	# Now the unit was stripped off. Converting numeric part to float and return the value.
	return float(thick)/umod, unit 
