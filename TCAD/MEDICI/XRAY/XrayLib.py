"""
	XrayLib.py
	
	Library for XrayConc.py
	Library contains various functions for XrayConc.py:
	From reading files to generating medici input deck.
	
"""
import numpy as np
import os as os

# Preprocessor for Debug
# Set as debug_mode = 1 for complete debug.
debug_mode = 0

def readfile(fn, skip_rows=0):
	"""
		Reads PENELOPE data.
		
		Input:
			data (including path), # of rows to skip
		
		Output:
			depth, X-ray intensity
			
	"""
	return np.loadtxt(open(fn, "r"), \
			skiprows=skip_rows, \
			unpack=True \
			); 
	
# Functions for plt. formatting
def tick_num_format(x, multi=0):
	label_text = str(float(x/10**multi))
	if multi == 3:
		label_text = label_text + 'k'
	elif multi == 6:
		label_text = label_text + 'M'
	elif multi == 9:
		label_text = label_text + 'G'
	elif multi == -2:
		label_text = label_text + 'c'
	elif multi == -3:
		label_text = label_text + 'm'
	elif multi == -6:
		label_text = label_text + '$\mu$'
	elif multi == -9:
		label_text = label_text + 'n'
	elif multi == -12:
		label_text = label_text + 'p'
	elif multi == -15:
		label_text = label_text + 'f'
	else:
		pass

	return label_text

# Smoothing function (just for 1st and 2nd order
# differentiation smoothing)
def smooth(x, window_len=10, window='hanning'):
	"""
	Smooth the data using a window with requested size.

	input:
        x: the input signal as array
		window_len: the dimension of window from 'flat', 
		'hanning','hamming','bartlett','blackman' 
		flat window will produce a moving average smooting.

	output:
	    the smoothed signal as array

    example:
	    t=linspace(-2,2,0.1)
		x=sin(t)+randn(len(t))*.1
		y=smooth(x)

	see also:
	    numpy.hanning, numpy.hamming, numpy.bartlett, 
	numpy.blackman, numpy.convolve, scipy.signal.lfilter

	"""
	if np.array(x).ndim != 1: 
		raise ValueError, \
		"smooth only accepts 1 dimension arrays."

	if np.array(x).size < window_len: 
		raise ValueError, \
		"Input Vector needs to be larger than window size."

	if window_len < 3:
		return x

	if not window in \
		['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, \
			"Window need to be one of 'flat', 'hanning', \
			'hamming', 'bartlett', 'blackman'"

	s = \
		np.r_[2 * x[0] - x[window_len:1:-1], \
		x, 2 * x[-1] - x[-1:-window_len:-1]]

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
# Reads out a value from array defined function 
# using linear interpolation
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
		raise ValueError, \
			"Function should be defined as array."
	if np.array(fnx).size != np.array(fn).size:
		raise ValueError, \
			"X axis data and Y axis data should match in dimension"

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
		raise ValueError, \
			"Function should be defined as array."
	if np.array(fnx).size != np.array(fn).size:
		raise ValueError, \
			"X axis data and Y axis data should match in dimension"

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
				fny.append((highfn - lowfn) / \
					(higher - lower) * (value - lower) + lowfn)
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
	fnds = smooth(fnd, 20)
	fndd = np.diff(fnds, 1)
	fndds = smooth(fndd, 30)

	zerosd = findzero(np.linspace(min(x), max(x), len(x) - 1), \
		fnds, threshold)

	maxx = []
	for zero_pos in zerosd:
		if readfunction(zero_pos, \
			np.linspace(min(x), max(x), len(x) - 2), fndds) > 0:
			pass
		else:
			maxx.append(zero_pos)

	return maxx


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

	zerosd = findzero(np.linspace(min(x), max(x), len(x) - 1), \
		fnds, threshold)

	minn = []
	for zero_pos in zerosd:
		if readfunction(zero_pos, \
			np.linspace(min(x), max(x), len(x) - 2), fndds) < 0:
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
		raise ValueError, \
			"Input must be string formatted integer or float!"

	# Stripping out space
	thick = thick.replace(' ','')

	# If the input is negative?
	if thick[0] == '-':
		raise ValueError, \
			"Negative Thickness cannot be realistic!"

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
		raise ValueError, \
			"The input must have valid unit(A,nm,um,mm,cm) \
			or has no unit!"

	# Now the unit was stripped off. 
	# Converting numeric part to float and return the value.
	return float(thick)/umod, unit	
	
	
# Extracting EHP concentrations.
def getEHP(pEnergy, aCoeff, sThick):
	"""
		Obtains Electron Hole Pair profile in silicon substrate
		according to given attenuation coeff. data.
		
		Input:
			Energy (only works with 20keV atm), Attenuation Coeff., 
			Substrate Thickness
		
		Output:
			EHP concentration in (/cm^3)
		
	"""
	# Using Mammography X-ray profile.
	# mean exposure: 12 mR (=P)
	P = 12
	# exposure range: 0.6 - 240 mR
	# Conversion of the exposure to number of photons/mm2,
	# at 20 keV, Theta = 5.05552e+07 ph/mm2/R
	# TODO: Figure out how this number came out.
	Theta = 5.05552e+07
	# Number of photons incident on the device (/area)
	numPhotonInc = Theta*P
	
	# Electron-hole pair creation energy (Eg is the bandgap of Silicon)
	Eg = 1.15
	W = 2.8*Eg + 0.5
	
	# TODO: Continue on this... when have time...
	
	return depth, EHP
	
# Shink the # of data points under 100	
def shrinkdp(depth, photon):
	depth = np.array(depth, 'f')
	photon = np.array(photon, 'f')	
	if depth.shape[0] != photon.shape[0]:
		raise ValueError, "Wrong insertion on shrinkdp function at XrayLib!!"
	
	if depth.shape[0] <= 100:
		print "No need to shrink the data points."
		return depth, photon
	else:
		print "Shrinking the # of data points under 100."
		print "Using linear approximation."
		depth_new = np.linspace(np.min(depth), np.max(depth), 100)
		photon_new = readfunction_a(depth_new, depth, photon)
		
		return depth_new, photon_new
	
	
# Generate EHP contration coordinate from X-ray data in Medici form.
def Medici_photogen_codegen(depth, Xphoton, xp_filename, \
	dim_x_min, dim_x_max, dim_y_min, dim_y_max, \
	datapoints_x, datapoints_y):
	"""
		Generates EHP concentration on each position to get Medici code.
		
		Input:
			depth, X-ray photon concentration as pair, 
			Filename to be opened in the Medici code,
			dimension x min, dimension x max,
			dimension y min, dimension y max (in um)
			(Simulation mesh area that needs X-ray photons be placed),
			datapoints_x, datapoints_y
			
		Output:
			Medici form of EHP concentration string.
			
	"""
	# Floatifying resolution
	photogen_instruction = []
	for j in np.linspace(dim_x_min,dim_x_max,int(datapoints_x)):
		x_text = \
			"X.START="+str(float(j))+" X.END="+str(float(j))+" "
		y_text = \
			"Y.START="+str(float(dim_y_min))+" Y.END="+str(float(dim_y_max))+" "
		"""
		photogen_fase = \
			"PHOTOGEN "+x_text+y_text+\
			"LETFILE=\""+xp_filename+"\""+" "+\
			"DCHR=0.2 "+\
			"UNIFORM ^CLEAR"+'\n'
		"""
		photogen_fase = \
			"PHOTOGEN "+x_text+y_text+\
			"LETFILE=\""+xp_filename+"\""+" "+\
			"DCHR=0.2 "+\
			"T0=3.0E-12 TC=1.5E-12 GAUSS ^CLEAR"+'\n'
		
		photogen_instruction.append(photogen_fase)		
	
	return photogen_instruction






