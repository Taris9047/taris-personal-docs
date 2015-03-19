"""

	tauc_setup.py
	
	Purpose:
		Provides a setup window for taucGUI.py application.
	
	Dependency:
		wx

"""

import wx
import numpy as np

class TaucSetup(wx.Frame):
	
	def __init__(self, parent, id, title):
		wx.Frame.__init__(self, parent, id, title, wx.DefaultPosition, wx.Size(325, 320))
		self.threshold_scope = np.array([0, 600, 800])
		self.thickness = 250
		self.refractive_index_substrate = 1.51
		self.fitting_method = 'log'
		
		