'''
Created on Jun 14, 2009

@author: taris

Includes acutal functions to draw Tauc plot.

'''

from tauc_lib import *
from tauc_fit import *
import matplotlib.pyplot as plt


# Returning refractive index along wavelength
# Transmittance value is not % in here.
def refrac(wl, tm, threshold=[200,350,800],s = 1.52, fit_type='exp'):
	# Refractive index of Corning 1737
	# s = 1.52
	# tm Fitting to derive TM and Tm
	
	peaks, maxfit, bottoms, minfit = findpeak_rev(wl, tm, delta=.5)
	
	##################################################
	# Initializing fitting parameter for TM/Tm fits. #
	##################################################
	if fit_type == 'exp':
		#p0 for exponential fit
		p0 = [max(maxfit), 1, threshold[0], 50]
	elif fit_type == 'poly':
		#p0 for polynomial fit
		p0 = [1,1,1,1,1,1]
	elif fit_type == 'expAssoc':
		#p0 for ExpAssoc fit
		p0 = [max(maxfit), 1, threshold[0], 50, \
			max(minfit), 1, threshold[0], 50]
	elif fit_type == 'log':
		#p0 for Logarithm fit
		p0 = [0, 1, 200]
	else:
		p0 = [0]
	
	# for Debug, double-checking peaks(bottoms) and maxfit(minfit)
	if debug_mode == 1:
		plt.figure(110)
		plt.plot(wl,tm, 'b-')
		plt.plot(peaks,maxfit,'o')
		plt.plot(bottoms,minfit,'.')
		plt.axis([0, max(wl), 0, 100])
		plt.title("Peaks and Bottoms for Transmittance")
		plt.xlabel("Wavelength (nm)")
		plt.ylabel("Transmittane (%)")
		plt.draw()
	
	# Running the fitting process!
	print "Beginning fitting process..."
	# Looking for initial point.
	min_tm_threshold = 5 # TM 5% as minimal point.
	init_index = 0
	for tmval in tm:
		if tmval > min_tm_threshold:
			init_index += 1
			break
			
	print "Setting up fixed ponint at, ", str(wl[init_index]) , "nm."
	TM = TM_fit(peaks, maxfit, wl, \
	 	 wl[init_index], tm[init_index], p0, fit_type)
	Tm = TM_fit(bottoms, minfit, wl, \
		 wl[init_index], tm[init_index], p0, fit_type)
	
	# for Curve Fitting Debug.
	if debug_mode == 1:
		plt.figure(100)
		plt.plot(wl, TM, 'r-')
		plt.plot(peaks,maxfit,'ro')
		plt.plot(wl, Tm, 'g-')
		plt.plot(bottoms,minfit,'go')
		plt.axis([0, max(wl), 0, 100])
		plt.title("Curve Fitting for Transmittance")
		plt.xlabel("Wavelength (nm)")
		plt.ylabel("Transmittane (%)")
		plt.draw()
	
	# for Transparent Region,
	n_trans = []
	print "Calculating refractive index for Transparent region."
	for wval in wl:
		if wval <= threshold[2]:
			pass
		else:
			Tmval = readfunction(wval,wl,Tm)/100
			M = 2*s/Tmval - (s**2+1)/2
			n_trans.append(np.sqrt(M+np.sqrt(M**2-s**2)))
	
	# for Medium/Weak Absorption Region,
	n_ma = []
	print "Calculating refractive index for Medium/Weak Absorption region."
	for wval in wl:
		if wval <= threshold[1] or wval > threshold[2]:
			pass
		else:
			TMval = readfunction(wval,wl,TM)/100
			Tmval = readfunction(wval,wl,Tm)/100
			N = 2*s*(TMval-Tmval)/(TMval*Tmval) + (s**2+1)/2
			n_ma.append(np.sqrt(N+np.sqrt(N**2-s**2)))
	
	# for Strong Absorption region,
	n_str = []
	p0 = [3e5,3]
	print "Initiating refractive index fitting procedure."
	print "Initial model: y = ", p0[0],"/wl**2 + ", p0[1]
	parasq = leastsq(refrac_residual, \
		p0, args=(n_ma+n_trans, wl[len(wl)-len(n_ma+n_trans):]))
	print "Fitted model: y = ", \
		parasq[0][0],"/wl**2 + ", parasq[0][1]
	print "Estimating refractive index."
	
	for wval in wl:
		if wval > threshold[1]:
			pass
		else:
			n_str.append(refrac_fit(parasq[0],wval))
	
	if debug_mode == 1:
		print "Drawing refractive index."
		plt.figure(101)
		plt.plot(wl, n_str + n_ma + n_trans)
		plt.title("Refractive Index")
		plt.xlabel("Wavelength (nm)")
		plt.ylabel("Refractive Index (A.U.)")
		plt.draw()
	
	return n_str + n_ma + n_trans

# absorption
# obtaining alpha(absorption coeff.) along wavelength
def absorption(wl,tm,threshold=[300,450,950], \
	thickness=360e-9,s=1.52, fit_type='exp'):
	'''
	Extracts absorption coefficient(alpha) for given spectra.
	
	Inputs :
		wavelength, transmission spectra, 
		threshold for tm, thickness, 
		refractive coefficient of substrate, fitting type.
	
	Output :
		[wavelength, absorption coeff.]
	
	'''
	
	print "Setting up fit type: ", fit_type
	print "Setting up threshold range: ", threshold
	
	refractive_index = refrac(wl,tm,threshold,s,fit_type)
	wval_index = 0
	alpha = []
	for wval in wl:
		n = refractive_index[wval_index]
		
		# Printing current refractive index for debug.
		if debug_mode == 1:
			print "n for wl(",wval," nm) is, ", n
		
		# Defining A,B,C,D
		A = 16*(n**2)*s
		B = (n+1)**3*(n+s**2)
		C = 2*(n**2-1)*(n**2-s**2)
		D = ((n-1)**3)*(n-s**2)
		psi = 4*np.pi*n*thickness/(wval*1e-9)
		
		if debug_mode == 1:
			print "A, B, C, D, cos(psi) are, ", A, B, C, D, np.cos(psi)
			print "Thickness used in calcluation : ", thickness
		
		# Solving equation for x using tm.
		# In this case, we have to solve a quadratic equation,
		# x^2 - (C/D*cos(psi)+A/TD)x + (B/D) = 0
		# Therefore,
		a = 1
		b = -((C/D)*np.cos(psi) + A/((tm[wval_index]/100)*D))
		c = B/D
		
		if debug_mode == 1:
			print "a,b,c for wval(",wval,") are, ", a,b,c
		
		# Obtaining determinant for the 2nd order polynomial equation.
		determinant = b**2 - 4*a*c
		
		if debug_mode == 1:
			print "At wavelength ", wval , "nm", ", Determinant : ", determinant
		
		# If the equation has imaginary solutions,
		if determinant < 0:
			print "** Unable to reach Real space solutions. **"
			print "Using 0 as alpha"
			alpha.append(0)
		
		# If the equation has one solution,
		elif determinant == 0:
			x = -b/(2*a)
			alpha_index = -np.log(x[0])/thickness
			
			if debug_mode == 1:
				print "Alpha on ", wval, "nm is,", alpha_index
		
			alpha.append(alpha_index)
		
		# If the equation has two solutions, select one between 0 and 1. 
		else:
			sols = \
				[(-b-np.sqrt(determinant))/(2*a), \
				 (-b+np.sqrt(determinant))/(2*a)]
			
			# One if the solutions must be fall into between 0 and 1.
			x = np.array([])
			for sol in sols:
				if sol > 0.0 and sol < 1.0: 
					 x = np.append(x,sol)
			
			if len(x) == 1:
				alpha_index = -np.log(x[0])/thickness
				if debug_mode == 1:
					print "Alpha on ", wval, "nm is,", alpha_index
				alpha.append(alpha_index)
			else:
				print "At wavelength ", wval, "nm"
				print "Solutions are ", str(sols[0]), str(sols[1])
				print "** Unable to resolve reasonable absorbance value. **"
				print "Appending alpha with ", str(-np.log(1)/thickness)
				alpha.append(-np.log(1)/thickness)

		# increasing wavelength array index
		wval_index += 1
		
	if debug_mode == 1:
		print alpha
		print "Drawing alpha."
		plt.figure(102)
		plt.plot(wl,alpha)
		plt.title(r'$\alpha$ vs. Wavelength')
		plt.xlabel(r'Wavelength (nm)')
		plt.ylabel(r'$\alpha$ (A.U.)')
		plt.draw()
	
	return alpha


def tauc_plot(wl, alpha):
	"""
	Drawing tauc plot using fitted data from above.
	
	Inputs:
		wl, alpha(absorption coeff)
	
	Outputs:
		energy, sqrt(hv*alpha)
	"""
	# Physical Constants (unit:MKS)
	h = 4.13566733e-15
	electron = 1.60217646e-19
	lightspd = 299792458
	
	tauc_alpha = np.array([])
	index = 0
	for wval in wl:
		v = lightspd/(wval*1e-9) # Converting the Wavelength to Frequency
		tauc_alpha = np.append(tauc_alpha, np.sqrt(h*v*alpha[index]))
		index += 1
	
	if debug_mode == 1:
		print alpha
		print tauc_alpha
	
	energy = []
	for wval in wl:
		v = lightspd/(wval*1e-9)
		energy.append(h*v)
	
	return energy, tauc_alpha
