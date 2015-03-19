'''
Created on Apr. 7, 2011

@author: Taylor Shin

Threshold voltage processing and plotting for SiDetector Simulation.

'''

import os, sys
import numpy as np
# MS Excel I/O modules, make sure xlrd, xlwt, xlutils are installed.
import xlrd
import xlwt
import xlutils

# Stores list of markers for matplotlib plotting
def marker_bank(n=15, col=0):
	mk_bank = ['1', '2', '3', '4', '+', ',', '_', '|', '*', '>', '<', \
			  'D', 'H', 'd', 'h', 'o', 'p', 's', 'v', 'x']
	col_bank = ['', 'k', 'r', 'g', 'b', 'c', 'y']
	return col_bank[col] + mk_bank[n]

# Polynomial fit model using least square method. Returns Ids value
def trans_orig_fit(v_gs, delta_vth=0):
	vth = 5.0886
	if v_gs > vth+delta_vth:
		return .41*1e-6*(v_gs-(vth+delta_vth))**2 
	else:
		return 0

# Generating array based data for transfer curve Ids
def get_val_trans(v_gs, delta_vth=0):
	Ids = []
	for volt_val in v_gs:
		Ids.append(float(trans_orig_fit(volt_val,delta_vth)))
	return Ids

# Save 2D array as MS Excel 97 compatible spreadsheet.
def save_excel(ary_in, xlsname, sheet_name='Sheet 1'):
	'''
	Saves 2D array as MS Excel 97 compatible spreadsheet.
	Ver. 0.1 at this moment. Only saves one sheet.
	
	Input:
		2D Array, xls file name to be saved.
		
	Output:
		Filename for saved xls sheet.
	'''
	if str.find(xlsname, '.xls') == -1:
		xlsname = xlsname+'.xls'

	xls_wb = xlwt.Workbook()
	xls_ws = xls_wb.add_sheet(sheet_name)
		
	if ary_in == []:
		print "Input array was empty. Generating an empty xls file."
	else:		
		for index_col in range(0, len(ary_in)):
			for index_row in range(0, len(ary_in[index_col])):
				xls_ws.write(index_row, index_col, ary_in[index_col][index_row])
		
	xls_wb.save(xlsname)	
	print str(xlsname), " has been saved."
	return xlsname
	
# returns whether the input string is numeric format or not.
def isnumeric(value):
	return str(value)\
	.replace("e","").replace("E","").replace(".", "")\
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

# Reads data from spectroscopy output to python workspace.
def ront_readfile(fn):
	"""
	Reads csv form of Rontgen data sheet.

	Input: 
		Filename (including path)

	Output: 
		depth, conc

	"""
	originalFile = open(fn, "r")
	# adjusting filename for further use
	fn = fn.replace(".txt", "")
	depth_from_surface = []
	concentration = [] 
	lineIndex = 0
	list_exposure = []

	for line in originalFile.readlines():
		depth = []
		conc = []
		if line.find("\r") != -1:
			line = line.replace("\r","")
		if line.find("\n") != -1:
			line = line.replace("\n","")			
		if lineIndex == 0:
			tmp_list = line.split('\t')
			for list_element in tmp_list[2:]:
				# Double Quote Removal
				list_element.replace('\'', '', 2)
				list_element.replace('\"', '', 2)
				'''
				# Converting to float
				if list_element.find('R') != -1:
					list_element.replace('R', '')
					print list_element
					if list_element.find('m') != -1:
						adjust = 1e-3
						list_element.replace('m','')
						list_element.replace(' ','')
						list_element = float(list_element)*adjust
					elif list_element.find('u') != -1:
						adjust = 1e-6
						list_element.replace('u','')
						list_element.replace(' ','')
						list_element = float(list_element)*adjust
					else:
						pass	
				'''
				list_exposure.append(list_element)
		else:
			line = line.replace('\"', '', 10000)
			line = line.replace('  ', ' ', 10000)
			line = line.replace(' ', '\t')
			line = line.replace(',', '\t')
			line = line.replace("\t\t", "\t", 10000)
			data = line.split()
			depth = float(data[0])
			for datum in data[2:len(data)]:
				conc.append(float(datum))
				
		#print line.split("\t")
		if isnumeric(depth) == False:
			print "Invalid point at ", depth, " (m), false data."
		else:
			depth_from_surface.append(float(depth))
			concentration.append(conc)
		
		lineIndex += 1

	print "Total data points: " + str(lineIndex-1) + "\n"

	return list_exposure, np.array(depth_from_surface), np.array(concentration)

# Evaluate Function with linear interpolation.
def eval_func(x, fnx, fn): 
	"""
	Input:
		x: X value input.
		fnx: X component of function.
		fn: function itself.
	
	Output:
		Y value at x in X of function.
		
	"""
	if np.array(fnx).ndim != 1 or np.array(fn).ndim != 1:
		raise ValueError, "Function should be defined as 1D array."
	if np.array(fnx).size != np.array(fn).size:
		raise ValueError, "X axis data and Y axis data should match in dimension"
	
	if x < fnx[0] or x > fnx[len(fnx)-1]:
		print "X = ", str(x)
		raise ValueError, "X is out of given range. Make sure x sits within fnx range."
	else:
		index = 0
		for lower in fnx:
			if index == len(fnx)-1:
				return fn[len(fn)-1]
			elif index == 0 and x == fnx[0]:
				return fn[0]
			else:
				higher = fnx[index + 1]
				lowfn = fn[index]
				highfn = fn[index + 1]
				if x > higher:
					index = index + 1
				else:
					return (highfn - lowfn) / (higher - lower) * (x - lower) + lowfn


