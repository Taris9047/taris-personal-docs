'''
Created on Apr. 6, 2011

@author: Taylor Shin

Threshold voltage processing and plotting for SiDetector Simulation.

'''
# Importing Tk graphics library. Not used for now.
#from Tkinter import *
#import tkFileDialog as tkFD
#import tkSimpleDialog as tkSD
# Sorting purpose.
import operator
# Importing numerical calculation and data plotting library.
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as plt_font
# Lastly, importing own library.
from dynrng_lib import *

# Main Function
def main():
	# Reading the data file
	filename = './380_interest_uR.csv'
	print "Reading the EHP distribution data from ..."
	print filename
	# Reading data from file.
	ronts, depth, conc = ront_readfile(filename)
		
	# Re-arranging data on EHP concentration
	print "Re-arranging data on EHP concentration."
	EHP = dict()
	for index, ront in enumerate(ronts):
		ront = ront.replace('\"','',2)
		ront = ront.replace('\'','',2)
		ront = ront.replace('R','')
		ront = ront.replace(' ','', 2)
		if ront.find('m') != -1:
			mod_unit = 1e-3
			ront = float(ront.replace('m',''))*mod_unit
		elif ront.find('u') != -1:
			mod_unit = 1e-6
			ront = float(ront.replace('u',''))*mod_unit
		else:
			pass
		ront = float(ront)
		EHP[ront] = conc[:,index]
	
	print "Calculating delta_Vth according to EHP concentration data."
	# Unit conversion for depth (m --> cm)
	device_depth = depth * 100
	# Actually, the area of device is L = 15 um and W = 1 um. 
	# Thus, in cm-2,
	Device_Area = 15e-4 * 1e-4
	
	# Number of EHP, (cm-2) 
	EHP_TOTAL = dict()
	for ront, EHP_list in EHP.iteritems():
		EHP_TOTAL[ront] = np.trapz(EHP[ront], device_depth)

	# Evaluating Vth shift
	q = -1.602e-19
	depletion_depth = 80e-4
	eps_r_silicon = 11.69
	eps_0 = 8.854e-14
	Cap_dep = eps_r_silicon*eps_0*(Device_Area/depletion_depth)
	print "Charge of an electron is ", q, " C"
	print "Depletion Depth : ", depletion_depth * 1e4, " um"
	print "Relative Permittivity of Silicon : ", eps_r_silicon
	print "Permittivity : ", eps_0, " F/cm"
	print "Thus, the capacitance is ", Cap_dep, " F"
	
	DELTA_Vth = dict()
	for ront, ehp_total in EHP_TOTAL.iteritems():
		DELTA_Vth[ront] = ehp_total*Device_Area*q/Cap_dep
		print "At ", ront, \
			", Vth varition is,", ehp_total*Device_Area*q/Cap_dep
	
	print "delta_Vth has been processed. Processing Plot data."
	
	# Generating expected EHP graph based on fitting.
	# See Fitting_Results.txt in Initial Point folder.
	# Reading the data file
	filename = './SiTFTTRANSDK.txt'
	print "Reading the Transfer Curve from ..."
	print filename
	# Reading data from file.
	volt, ids_orig = trans_readfile(filename) 
	
	# Setting up plot data
	print "Generating predicted plots."
	vgs_fit_min = 0
	vgs_fit_max = 30
	vgs_fit = np.array(range(vgs_fit_min,vgs_fit_max+1))
	ids_fit_orig = np.array(get_val_trans(vgs_fit, 0))
	Ids_FIT = dict()
	for ront, delta_Vth_val in DELTA_Vth.iteritems():
		Ids_FIT[ront] = get_val_trans(vgs_fit, delta_Vth_val)
	
	# Plotting Trasfer Curves
	Ids_norm = 1e6
	plt.figure()
	orig_plt, = plt.plot(volt, ids_orig*Ids_norm, 'rs', markersize=10)
	orig_fit, = plt.plot(vgs_fit, ids_fit_orig*Ids_norm, 'b-')

	PLT = dict()
	xls_array_vth = [[]]
	ront_modifier = 1e3
	i = 0
	PLT_LEGEND_pointer = [orig_plt, orig_fit]
	PLT_LEGEND_text = [r'Original', r'Least Square Fit']
	sorted_Ids_FIT = \
		sorted(Ids_FIT.iteritems(), key=operator.itemgetter(0))
	xls_array_vth[0] = ['Vgs (V)']+list(vgs_fit)
	for ront, ids_fit_list in sorted_Ids_FIT:
		xls_array_vth\
			.append([str(ront*ront_modifier)+' mR']+list(ids_fit_list))
		PLT[ront], = \
			plt.plot(vgs_fit, np.array(ids_fit_list)*Ids_norm, \
			marker_bank(i), markersize=10)
		PLT_LEGEND_pointer.append(PLT[ront])
		PLT_LEGEND_text.append(\
			r'Exposure ' + str(ront*ront_modifier) + ' mR')
		i += 1
	plt.legend(PLT_LEGEND_pointer, PLT_LEGEND_text, \
			loc='best', prop=plt_font.FontProperties(size=8))
	plt.xlim((0, 30))
	plt.ylim((0, 400))
	plt.title(\
		r'Transfer Curves for Various X-ray Exposures at Vds = 20 V', fontsize=16)
	plt.xlabel(r'Gate Voltage (V)', fontsize=16)
	plt.ylabel(r'Drain Current ($\mu$A/$\mu$m)', fontsize=16)
	ax = plt.gca()
	ax.set_yscale('log')
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	#plt.show()
	
	# Saving the graph data to xls format
	xlsname = r'Transfer_Curve_Shift.xls'
	save_excel(xls_array_vth, xlsname, 'Transfer_Curves_under_Exposure')
	
	# Sensitivity Evaluation
	print "Sensitivity Evaluation ...."
	DIFF_IDS = dict()
	for ront, ids_list in sorted_Ids_FIT:
		DIFF_IDS[ront] = np.array(ids_list) - np.array(ids_fit_orig)
	sorted_DIFF_IDS = sorted(DIFF_IDS.iteritems(), \
		key=operator.itemgetter(0))
		
	# Drawing Sensitivity Plot
	vgs_sensitivity_point = [1, 10, 20, 30]
	ronts_sensitivity = []
	for ront, temp in sorted_Ids_FIT:
		ronts_sensitivity.append(ront)
	ronts_sensitivity = np.array(ronts_sensitivity)
	
	SENSI_LIST = dict()
	for vgs_val in vgs_sensitivity_point:
		SENSI_LIST[vgs_val] = []
		for ront, diff_ids_list in sorted_DIFF_IDS:
			SENSI_LIST[vgs_val]\
				.append(eval_func(vgs_val,vgs_fit,diff_ids_list))
		SENSI_LIST[vgs_val] = np.array(SENSI_LIST[vgs_val])
	
	# Plotting Sensitivity
	print "Plotting Sensitivity."
	xls_array_sensitivity = [[]]
	xls_array_sensitivity[0] = ['Exposure (R)']+list(ronts_sensitivity)
	plt.figure()
	# Generating sorted list for dict SENSI_LIST
	sorted_SENSI_LIST = sorted(SENSI_LIST.iteritems(),\
			key=operator.itemgetter(0))
	i = 0
	for vgs_val, sensi_list_var in sorted_SENSI_LIST:
		print "Plotting at Vgs = ", str(vgs_val), " V"
		xls_array_sensitivity\
			.append([str(vgs_val)+" V"]+list(sensi_list_var))
		plt.plot(ronts_sensitivity*ront_modifier, \
			sensi_list_var*Ids_norm, marker_bank(i),\
			markersize=10, label='at Vgs = '+str(vgs_val)+' V')
		i += 1
	plt.title(r'Sensitivity Plot on Various Exposures at Vds = 20 V', fontsize=16)
	plt.xlabel(r'Exposures (mR)', fontsize=16)
	plt.ylabel(\
		r'EHP Induced Drain Current ($\mu$A/$\mu$m)', fontsize=16)
	plt.legend(loc='best', prop=plt_font.FontProperties(size=8))
	ax = plt.gca()
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	plt.show()
	
	# Saving the graph data to xls format
	xlsname = r'Sensitivity.xls'
	save_excel(xls_array_sensitivity, xlsname, 'Sensitivity_Plot')
	
# running the program   
if __name__ == "__main__":
	main()
