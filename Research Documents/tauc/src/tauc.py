#!/usr/bin/env python
"""
	tauc.py
	
	Purpose:
		To draw Tauc's plot using transmittance data from
		optical spectra.
	
	Dependency:
		tauc_func, matplotlib, Tkinter, tkFileDialog, os, sys

"""
from tauc_func import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from Tkinter import *
import tkFileDialog as tkFD
import tkSimpleDialog as tkSD
import os, sys

def main(): 
	# Reading file, filename contains full path to the selected spectra source.
	root = Tk()
	root.withdraw()
	filename = \
		tkFD.askopenfilename( \
			defaultextension="txt", initialdir="./data", \
			title="Select an Optical Spectra Data.")
	if os.path.isfile(filename) != True:
		raise ValueError, "Invalid file or cannot find designated file!"
	else:
		print "Reading the spectrum data from ..."
		print filename
		# Reading data from file.
		wl, tm = tauc_readfile(filename)
	
	# Thickness?
	root = Tk()
	root.withdraw()
	thick = tkSD.askstring( \
		"Thickness?", "Preferred Unit: nm", initialvalue="300")
	thickness, unit = thick_input(thick)
	print "The film thickness: "+str(thickness)+" nm"       

	# Setting outfile name
	path_divided = filename.split("/")
	if os.name == 'posix':
		appendix_dir = "./tauc_plot/"
	else:
		appendix_dir = ".\\tauc_plot\\"

	filename_divided = path_divided[len(path_divided)-1].split(".") # filename_divided[0] contains the 'name' part of the .txt file. 
	pdf_outf = filename_divided[0] + '.pdf'
	input_file_name = path_divided[len(path_divided)-1]
		
	# Plotting Transmittance data
	pp = PdfPages(appendix_dir+pdf_outf)
	plt.figure(0)
	plt.plot(wl,tm, linewidth=2)
	plt.title(r'Optical Transimission Spectra of '+str(input_file_name) \
		+r" ("+str(thickness)+r"-nm-thick)",fontsize=16)
	plt.xlabel(r'Wavelength (nm)', fontsize=16)
	plt.ylabel(r'Transmittance (%)', fontsize=16)
	plt.axis([min(wl)+1,max(wl),0,100])
	ax = plt.gca()
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	for tick in ax.yaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	plt.savefig(pp, format='pdf')
	
	# TODO: Fix refractive index part with threshold definition.
	# Setting up threshold.
	root = Tk()
	root.withdraw()
	threshold_input = tkSD.askstring( \
		"Threshold Setup Window", \
		"Preferred Unit: nm", initialvalue="0 300 600").split()
	for threshold_val in threshold_input:           
		threshold = [ float(threshold_input[0]), \
					  float(threshold_input[1]), \
					  float(threshold_input[2]) ]
	# re-arranging threshold if input is mistaken
	threshold.sort()
	
	print "Threshold setting for each regions."
	print "Strong absorption: ", str(wl[0]), "to", str(threshold[0])        
	print "Medium/Weak absorption: ", str(threshold[0]), \
		  "to", str(threshold[1])
	print "Transparent absorption: ", str(threshold[1]), \
		  "to", str(threshold[2])
	
	# Evaluating Absorption Coefficient
	alpha = absorption( \
		wl, tm, threshold=threshold, \
		thickness=(thickness*1e-9),s=1.52,fit_type='poly')
	
	# Drawing Tauc Plot
	energy, tauc_alpha = tauc_plot(wl,alpha)
	
	plt.figure(1)
	plt.plot(energy,tauc_alpha, linewidth=2.0)
	plt.title(r'Tauc Plot of '+str(input_file_name) \
		+r' ('+str(thickness)+r'-nm-thick)', fontsize=16)
	plt.xlabel(r'Energy (eV)', fontsize=16)
	plt.ylabel(r'Absorption Coeff.' \
		+ r'$\sqrt{\alpha h \nu}$'+r' (A.U.)', fontsize=16)
	ax = plt.gca()
	for tick in ax.xaxis.get_major_ticks():
		tick.label1.set_fontsize(14)
	for tick in ax.yaxis.get_major_ticks():
		tick.set_visible(False)
	plt.savefig(pp, format='pdf')
	pp.close()

	# writing tauc data to file
	outf_name = appendix_dir+filename_divided[0]+"_tauc.txt"
	fp = open(outf_name, "w")
	tauc_data = []
	for i in range(0,len(tauc_alpha)-1,1):
		tauc_data.append(str(energy[i])+"\t"+str(tauc_alpha[i])+"\r\n")
	fp.writelines(tauc_data)
	fp.close()
	
	print "...Tauc plot Data Saved at, ", outf_name

	# Opening PDF format results.
	print "Opening Plot File (PDF format.)"
	if os.name == 'posix':
		os.system("/usr/bin/xdg-open ./"+appendix_dir+pdf_outf)
	else:
		os.filestart(pdf_outf)

	print "Execution Finished!!!!"
	
	return 0

# Running main program
if __name__ == "__main__":
	main()
