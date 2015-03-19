'''
Created on Jun 14, 2009

@author: taris

Curve fitting feature for tauc plot project

'''

from tauc_lib import *
import numpy as np
from scipy.optimize import leastsq

# defining exponential type residual function for leastsq fitting process.
# exponential residual for Transmittance fitting
def exp_residual(p, y, x):
	"""
	Defining residual function for Transmittance fitting,
	y = A*(1 - k*exp((t-x)/alpha))
	"""
	err = y - exp_fit(p,x)
	return err

# 2 exponential associated residual for Transmittance fitting
def expAssoc_residual(p, y, x):
	"""
	Defining residual function for Transmittance fitting,
	y = A1*(1 - k1*exp((t1-x)/alpha1)) + A2*(1 - k2*exp((t2-x)/alpha2))
	"""
	err = y - expAssoc_fit(p,x)
	return err

# Polynomial residual for Transmittance Fitting.
def poly_residual(p, y, x):
	"""
	Defining residual for Transmittance Fitting
	y = Ax^5 + Bx^4 + Cx^3 + Dx^2 + Ex + F
	"""
	# Ensuring x is np data type
	x = np.array(x)
	err = y - poly_fit(p,x)
	return err

def log_residual(p, y, x):
	"""
	Defining residual for logarithmic function.
	y = A + B*log(x-C)
	"""
	x = np.array(x)
	err = y - log_fit(p,x)
	return err
	
# refractance formula residual
def refrac_residual(p, y, lamb):
	"""
	Defining residual function for refractive index model,
	n = A/lambda^2 + B
	"""
	
	err = y - refrac_fit(p, lamb)
	return err

# actual function to get data after curve fitting
def poly_fit(p0, x):
	return p0[0]*x**5 + p0[1]*x**4 + p0[2]*x**3 + p0[3]*x**2 + p0[4]*x + p0[5]
	
def exp_fit(p0,x):	  
	return p0[0] * (1 - p0[1] * np.exp((p0[2] - x) / p0[3]))

def expAssoc_fit(p0,x):
	return p0[0] * (1 - p0[1] * np.exp((p0[2] - x) / p0[3])) \
		+ p0[4] * (1 - p0[5] * np.exp((p0[6] - x) / p0[7]))

def log_fit(p0,x):
	return p0[0] + p0[1] * np.log(x - p0[2])

def peak_fit(p0, x , fit_type):
	if fit_type == 'exp':
		return exp_fit(p0,x)
	elif fit_type == 'expAssoc':
		return expAssoc_fit(p0,x)
	elif fit_type == 'poly':
		return poly_fit(p0,x)
	elif fit_type == 'log':
		return log_fit(p0,x)
	else:
		return x
	
# actual function to get data after curve fitting (refractive index)
def refrac_fit(p0, Lambda):
	return p0[0]/(Lambda**2) + p0[1]
	#return 3e5/(Lambda**2) + p0[1]


# Fitting peaks for TM
# use the array data with wl
def TM_fit(peaks, peaksfn, wl, fixed, fixedfn, p0, fit_type='exp'):
	'''
	TM_fit : 
		Performs curve fitting on extracted peaks and bottom values.
	Input :
		peaks, peaksfn, wl, fixed, fixedfn, p0, fit_type.
	Output :
		A vector of fitted curve which corresponds with wl.
	'''
	# Attaching fixed points to the peaks/bottoms data for curve fitting.			
	peaks_cv = np.append(fixed, peaks)
	peaksfn_cv = np.append(fixedfn, peaksfn)
		
	if fit_type == 'exp':
		para_sq = leastsq(exp_residual, p0, args=(peaksfn_cv, peaks_cv))
	elif fit_type == 'expAssoc':
		para_sq = leastsq(expAssoc_residual, p0, args=(peaksfn_cv, peaks_cv))
	elif fit_type == 'log':
		para_sq = leastsq(log_residual, p0, args=(peaksfn_cv, peaks_cv))
	else:
		para_sq = leastsq(poly_residual, p0, args=(peaksfn_cv, peaks_cv))
	
	#for debug
	if debug_mode == 1:
		print "TM_fit init: ", p0
		print "TM_fit: ", para_sq[0]
	
	# Generating TM
	TM = []
	for wlval in wl:
		if peak_fit(para_sq[0], wlval, fit_type) < 0:
			# Dropping negative 'fitted' Transmittance points
			# and replaces them with zero since they don't make any sense.
			TM.append(0)
		else:
			TM.append(peak_fit(para_sq[0], wlval, fit_type))
	
	return TM