#!/usr/bin/env python
"""
	taucGUI.py

	Purpose: 
		To draw Tauc's plot using transmittance data from 
		optical spectra.
		** GUI version based wxWigets python interpretation **

	Dependency:
		tauc_func_GUI, tauc_lib_GUI, matplotlib, wx, os, sys

"""
import wx
import numpy as np
import datetime
import os.path
import matplotlib
matplotlib.use('PDF')
import tauc_func_GUI as tf
import tauc_lib_GUI as tl
import matplotlib.pyplot as plt
import tauc_setup as setupGUI

ID_SAVELOG = 10001
ID_RUN = 10010
ID_SETUP = 10011
ID_PLOTSPEC = 10100

class MainWindow(wx.Frame):
	
	def __init__(self, filename='noname.log'):
		super(MainWindow, self).__init__(None, size=(640,400))
		self.filename = filename
		self.tauc_data_filename = []
		self.versionNum = "0.0.1"
		self.dirname = '.'
		self.fittype = 'log'
		self.energy = np.array([])
		self.tauc_alpha = np.array([])
		self.wavelength = []
		self.transmittance = []
		self.threshold = [0, 600, 800]
		self.thickness = 250
		self.IsLoaded = False
		self.CreateInteriorWindowComponents()
		self.CreateExteriorWindowComponents()
		self.displayLogo()

	def CreateInteriorWindowComponents(self):
		''' Create "interior" window components '''
		self.control = wx.TextCtrl(self, style=wx.TE_MULTILINE)
		self.control.SetEditable(False)
		#self.control = wx.ListBox(self)

	def CreateExteriorWindowComponents(self):
		''' Create "Exterior" window components '''
		self.CreateMenu()
		self.CreateStatusBar()
		self.SetTitle()
		
	def SetTitle(self):
		super(MainWindow, self).SetTitle('Reading spectra of %s'%self.filename)

	def CreateMenu(self):
		fileMenu = wx.Menu()
		for id, label, helpText, handler in \
			[(wx.ID_ABOUT, '&About', 'Information about this program', self.OnAbout),
			 (wx.ID_OPEN, '&Open', 'Open a new spectra file', self.OnOpen),
			 (wx.ID_SAVE, '&Save', 'Save current tauc data', self.OnSave),
			 (wx.ID_EXIT, 'E&xit', 'Terminate the program', self.OnExit)]:
			if id == None:
				fileMenu.AppendSeparator()
			else:
				item = fileMenu.Append(id, label, helpText)
				self.Bind(wx.EVT_MENU, handler, item)

		executeMenu = wx.Menu()
		for id, label, helpText, handler in \
				[(ID_PLOTSPEC, '&Plot Spectra', 'Plot fresh spectra from file', self.OnPlotSpec),
				 (ID_RUN, '&Run', 'Run the execution procedure', self.OnRun),
				 (ID_SETUP, 'Se&tup', 'Setup the procedure', self.OnSetup)]:
			if id == None:
				executeMenu.ApendSeparator()
			else:
				item = executeMenu.Append(id, label, helpText)
				self.Bind(wx.EVT_MENU, handler, item)

		menuBar = wx.MenuBar()
		for menuItems, label in \
			[(fileMenu, '&File'), (executeMenu, 'Execu&te')]:
			menuBar.Append(menuItems, label)
		self.SetMenuBar(menuBar)

	def defaultFileDialogOptions(self):
		return dict(message='Choose a spectra data', 
							defaultDir=self.dirname, wildcard='*.*')

	def askUserForFilename(self, **dialogOptions):
		dialog = wx.FileDialog(self, **dialogOptions)
		if dialog.ShowModal() == wx.ID_OK:
			userProvideFilename = True
			self.filename = dialog.GetFilename()
			self.dirname = dialog.GetDirectory()
			self.SetTitle()
		else:
			userProvideFilename = False
		dialog.Destroy()
		return userProvideFilename

	def DataCheck(self):
		if self.wavelength == [] and self.transmittance == []:
			self.control.AppendText("\n")
			self.control.AppendText("ERROR: No spectra data found!!!\n")
			self.control.AppendText("Select File --> Open to locate a suitable spectra.\n")
			return False
		else:
			return True

	def InitVariablesOnLoad(self):
		self.control.AppendText('\nInitializing...\n')
		self.energy = np.array([])
		self.tauc_alpha = np.array([])
		self.wavelength = []
		self.transmittance = []
		self.tauc_data_filename = []
		self.fittype = 'log'
		self.IsLoaded = True

	def displayLogo(self):
		self.control.AppendText("Tauc Plot Extractor Ver. "+self.versionNum+"\n")
		self.control.AppendText("Written by Kyung-Wook (Tylor) Shin\n"
								+"(kwshin@venus.uwaterloo.ca)\n"
								+"\n")
		self.control.AppendText("Select File --> Open for a new spectra.\n")

	#
	# Event Handlers:
	#
	def OnAbout(self, even):
		dialog = wx.MessageDialog(self, 'Tauc Plot extractor'
				 ' '+self.versionNum+'\n'
				 'Kyung-Wook Shin\n'
				 'kwshin@venus.uwaterloo.ca',
				 'About UW Tauc Plot', wx.OK)
		dialog.ShowModal()
		dialog.Destroy()

	def OnSave(self, even):
		self.DataCheck()
		if self.energy == np.array([]) and self.tauc_alpha == np.array([]):
			self.control.AppendText("\n")
			self.control.AppendText("ERROR: Tauc Data have not calculated! Please consider Execute -> Run.\n")
		else:
			self.tauc_data_filename = 'tauc_data_'+str(self.filename)+'.txt'
			tauc_data_file = open(os.path.join(self.dirname, self.tauc_data_filename), 'w')
			
			# Save calculated Tauc values into a file
			tauc_data_output = []
			for eV, alpha in self.energy, self.tauc_alpha:
				tauc_data_ouptput.append(str(eV)+'\t'+str(alpha)+'\n')	
			tauc_data_file.writelines(tauc_data_output)
			
			self.control.AppendText("Tauc Data was saved at :\n"
				+str(os.path.join(self.dirname, self.tauc_data_filename))
				+"\n\n")
			tauc_data_file.close()

	def OnOpen(self, even):
		# Initialize all variable to prevent confusion
		self.InitVariablesOnLoad()
		if self.askUserForFilename(style=wx.OPEN, **self.defaultFileDialogOptions()):
			self.control.AppendText('\n')
			self.control.AppendText('Reading data from ' 
									+ str(os.path.join(self.dirname, self.filename)) + '\n')
			spectrafile = open(os.path.join(self.dirname, self.filename), 'r')
			
			# Reading file to class variable wavelength, transmittance
			lineIndex = 0
			for line in spectrafile.readlines():
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
					if tl.isnumeric(wl) == False or tl.isnumeric(tm) == False:
						pass
					else:
						self.wavelength.append(float(wl))
						self.transmittance.append(float(tm))
						lineIndex += 1

			self.control.AppendText('Total data point(s) found: ' + str(lineIndex-1) + '\n')
			self.control.AppendText('\n')
			
			spectrafile.close()
			
		else:
			pass

	def OnRun(self, even):
		if self.DataCheck() != True:
			self.control.AppendText("\nCannot extract Tauc Plot from a void spectra\n")
		else:
			# Setting up parameters
			wlen = np.array(self.wavelength)
			tran = np.array(self.transmittance)
			# Notifying execution
			self.control.AppendText('\n\n Running Tauc Plot extraction \n')
			now = datetime.datetime.now()
			self.control.AppendText(now.strftime("%Y-%m-%d %H:%M:%S"))
			self.control.AppendText('\n')
			# Setting up thickness
			self.control.AppendText('Thickness is, ' + str(self.thickness) + ' nm\n')
			self.control.AppendText('Wavelength region separation :' 
									+ str(self.threshold[0]) + ' '
									+ str(self.threshold[1]) + ' '
									+ str(self.threshold[2]) + '\n')
			# TODO: Change the crude 'print' algorithm to more suitable for GUI applications
			self.energy, self.tauc_alpha = tf.tauc_plot(
					wlen,  
					tf.absorption(wlen, tran, self.threshold, self.thickness, 1.52, self.fittype))

			# Drawing Tauc Plot
			plt.figure()
			plt.plot(self.energy, self.tauc_alpha)
			plt.title(r"Tauc Plot of SiNx"+r" ("+str(self.thickness)+r"-nm-thick)", fontsize=16)
			plt.xlabel(r"Energy (eV)", fontsize=16)
			plt.ylabel(r'Absorption Coeff. '+ r'$\sqrt{\alpha h \nu}$'+r' (A.U.)', fontsize=16)
			ax = plt.gca()
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(14)
			for tick in ax.yaxis.get_major_ticks():
				tick.set_visible(False)
			plt.savefig('tauc')
			plt.draw()
			self.control.AppendText('Tauc Plot saved to '+'tauc.pdf\n')
			
	def OnSetup(self, even):
		# TODO: We need setup subroutine.
		if self.DataCheck() != True:
			self.control.AppendText('Setup is irrelevant if no spectra defined.\n')
		else:
			setup_frame = setupGUI.TaucSetup(None, -1, 'Tauc Setup Dialog')
			if setup_frame.Show() == ID_OK:
				# Apply the setup variables here
				pass
				
			# Then destroy the setup_frame since it has been completed its use.
			setup_frame.Destroy()

	def OnPlotSpec(self, even):
		# Plotting Transmittance data
		if self.DataCheck() != True:
			self.control.AppendText("Cannot draw the spectra plot!\n")
		else:
			# Drawing the spectra
			self.control.AppendText('\n\nDrawing Spectra...\n')
			plt.figure()
			plt.plot(np.array(self.wavelength),np.array(self.transmittance))
			plt.title(r'Optical Transimission Spectra of SiNx'
					  +r" ("+str(self.thickness)
					  +r"-nm-thick)",fontsize=16)
			plt.xlabel(r'Wavelength (nm)', fontsize=16)
			plt.ylabel(r'Transmittance (%)', fontsize=16)
			plt.axis([min(self.wavelength)+1,max(self.wavelength),0,100])
			ax = plt.gca()
			for tick in ax.xaxis.get_major_ticks():
				tick.label1.set_fontsize(14)
			for tick in ax.yaxis.get_major_ticks():
				tick.label1.set_fontsize(14)
			plt.savefig('spectra')
			plt.draw()
			self.control.AppendText('The Spectra was saved to '+'tauc.pdf\n')
			
	def OnExit(self, even):
		self.Close() # Close the Main Window
		
def main():
	app = wx.App(0)
	frame = MainWindow()
	frame.Show()
	app.MainLoop()
	
if __name__ == '__main__':
	main()
