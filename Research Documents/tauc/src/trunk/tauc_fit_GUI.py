#!/usr/bin/env python
'''
Created on Jun 14, 2009

@author: taris

Curve fitting feature for tauc plot project

'''

from tauc_lib_GUI import *
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
    y = A*(1 - k*exp((t-x)/alpha))
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
    return p0[0] * (1 - p0[1] * np.exp((p0[2] - x) / p0[3])) + p0[4] * (1 - p0[5] * np.exp((p0[6] - x) / p0[7]))

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
def refrac_fit(p0, lamb):
    #return p0[0]/(lamb**2) + p0[1]
    return 3e5/(lamb**2) + p0[1]


# Fitting peaks for TM
# use the array data with wl
def TM_fit(peaks, maxfit, wl, p0, fit_type='exp'):
    # reforming maxfit and peaks for curve fitting
    peaks_cv = []
    maxfit_cv = []
    index = 0
    for val in peaks:
        if maxfit[index] < .5:
            peaks_cv_init = val
            maxfit_cv_init = maxfit[index]
        else:
            peaks_cv.append(val)
            maxfit_cv.append(maxfit[index])
        index = index + 1
        
    peaks_cv = np.append(peaks_cv_init, peaks_cv)
    maxfit_cv = np.append(maxfit_cv_init, maxfit_cv)
        
    if fit_type == 'exp':
        para_sq = leastsq(exp_residual, p0, args=(maxfit_cv, peaks_cv))
    elif fit_type == 'expAssoc':
        para_sq = leastsq(expAssoc_residual, p0, args=(maxfit_cv, peaks_cv))
    elif fit_type == 'log':
        para_sq = leastsq(log_residual, p0, args=(maxfit_cv, peaks_cv))
    else:
        para_sq = leastsq(poly_residual, p0, args=(maxfit_cv, peaks_cv))
    
    #for debug
    if debug_mode == 1:
        print "TM_fit init: ", p0
        print "TM_fit: ", para_sq[0]
    
    # Generating TM
    TM = []
    for wlval in wl:
        if peak_fit(para_sq[0], wlval, fit_type) < 0:
            # Dropping zero Transmittance since it doesn't make any sense.
            TM.append(0)
        else:
            TM.append(peak_fit(para_sq[0], wlval, fit_type))
    
    return TM

# Fitting peaks for Tm
# use the array data with wl
def Tm_fit(bottoms, minfit, wl, p0, fit_type='exp'):
    # reforming maxfit and peaks for curve fitting.
    # (Removing useless data from peak data)
    bottoms_cv = []
    minfit_cv = []
    index = 0
    for val in bottoms:
        if minfit[index] < 5:
            bottoms_cv_init = val
            minfit_cv_init = minfit[index]
        else:
            bottoms_cv.append(val)
            minfit_cv.append(minfit[index])
        index = index+1
        
    bottoms_cv = np.append(bottoms_cv_init, bottoms_cv)
    minfit_cv = np.append(minfit_cv_init, minfit_cv)
        
    if fit_type == 'exp':
        para_sq = leastsq(exp_residual, p0, args=(minfit_cv, bottoms_cv))
    elif fit_type == 'expAssoc':
        para_sq = leastsq(expAssoc_residual, p0, args=(minfit_cv, bottoms_cv))
    elif fit_type == 'log':
        para_sq = leastsq(log_residual, p0, args=(minfit_cv, bottoms_cv))
    else:
        para_sq = leastsq(poly_residual, p0, args=(minfit_cv, bottoms_cv))
    
    #for debug
    if debug_mode == 1:
        print "Tm_fit init: ", p0
        print "Tm_fit: ", para_sq[0]
        
    # Generating TM
    Tm = []
    for wlval in wl:
        if peak_fit(para_sq[0], wlval, fit_type) < 0:
            Tm.append(0)
        else:
            Tm.append(peak_fit(para_sq[0], wlval, fit_type))
    
    return Tm

