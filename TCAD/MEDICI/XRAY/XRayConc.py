"""
	XRayConc.py
	
	Purpose:
		To calculate the concentration of electron hole pair
		from PENELOPE calculation. The final purpose is to
		generate a MEDICI input deck code according to the
		PENELOPE calculation.
	
	Dependency:
		numpy, matplotlib, os, sys

"""
import os, sys
import numpy as np
import XrayLib as XL
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Debug Preprocessor
debug_mode = 1

# index mapping with depth
def index_map_depth(index):
	if index == 0:
		return '300'
	elif index == 1:
		return '400'
	elif index == 2:
		return '500'
	elif index == 3:
		return '1000'
	else:
		raise ValueError, "Index must be one of: 0,1,2,3!"
		
# index mapping with energy
def index_map_energy(index):
	if index == 0:
		return '5keV'
	elif index == 0.5:
		return '10keV'	
	elif index == 1:
		return '15keV'
	elif index == 1.5:
		return '20keV'
	elif index == 2:
		return '25keV'
	elif index == 3:
		return '35keV'
	else:
		raise ValueError, "Index must be one of: 0,1,2,3!"

def main():
	# Showing version No.
	print "X-ray Concentration Ver. 0.0.3"
	print "Initializing..."
			
	# Processing EHP concentration according to attenuation coeff.
	# TODO: Familarize with the mechanism.
	# for now, we'll just use bruge force method
	
	# Select data (kai, yuan, nick)
	use_what = 'kai'
	
	
	if use_what == 'yuan':
		"""
			Info. on data variables for yuan's data
	
			depth: 
				2D array contains coordinate data from four data sources.
		
				depth[i,:] where i = 0,1,2,3 corresponds to 300, 400, 500, 1000
		
			X_photon: 
				3D array contains number of X-ray photons absorbed according
				to the coordinate data, depth.
			
				X_photon[i,j,:]
				where i = 0,1,2,3 corresponds to 300, 400, 500, 1000
				and, j = 0,1,2,3 corresponds to 5k, 15k, 25k, 30keV 
			
		"""
	
	
		print "Reading X-ray photon flux data. (from 300 to 1000 um) ..."
		temp_x_data = np.array([\
			XL.readfile('./yuan/EHP300um.txt', 5),\
			XL.readfile('./yuan/EHP400um.txt', 5),\
			XL.readfile('./yuan/EHP500um.txt', 5),\
			XL.readfile('./yuan/EHP1000um.txt', 5) \
			])
		depth = np.array([temp_x_data[0,0,:],\
				temp_x_data[1,0,:],temp_x_data[2,0,:],temp_x_data[3,0,:]],'f')
		X_photon = np.array([ \
				[temp_x_data[0,1,:],temp_x_data[0,2,:],temp_x_data[0,3,:],temp_x_data[0,4,:]], \
			    [temp_x_data[1,1,:],temp_x_data[1,2,:],temp_x_data[1,3,:],temp_x_data[1,4,:]], \
				[temp_x_data[2,1,:],temp_x_data[2,2,:],temp_x_data[2,3,:],temp_x_data[2,4,:]], \
				[temp_x_data[3,1,:],temp_x_data[3,2,:],temp_x_data[3,3,:],temp_x_data[3,4,:]], \
				],'f')
	
	
		# Converting cm^-2 to cm^-3
		#X_photon = X_photon*100*1e20
	
	# Uses Nick's data (Not yet implemented)
	elif use_what == 'nick':
		pass
	
	# Uses Kai's data
	else:
		temp_x_data = np.array(XL.readfile('./kai/ehp.txt'))
		if temp_x_data.shape[1] >= 100:
			depth, xptemp1 = XL.shrinkdp(temp_x_data[0,:], temp_x_data[3,:])
			depth, xptemp2 = XL.shrinkdp(temp_x_data[0,:], temp_x_data[2,:])
			depth, xptemp3 = XL.shrinkdp(temp_x_data[0,:], temp_x_data[1,:])
			depth = depth * 1e6
			X_photon = np.array([xptemp1, xptemp2, xptemp3])
		else:
			depth = np.array(temp_x_data[0,:], 'f') * 1e6
			X_photon = np.array([ \
				temp_x_data[3,:], temp_x_data[2,:], temp_x_data[1,:] \
				], 'f')
		
	
	# Checking the input data (EHP pairs) with visual inspection
	if debug_mode == 1:
		print "Plotting X-ray Photon Distribution..."
		plt.figure(2)
		if use_what == 'yuan':
			plt.plot(depth[0,:], X_photon[0,0,:], 'r-', \
			   		 depth[0,:], X_photon[0,1,:], 'g-', \
					 depth[0,:], X_photon[0,2,:], 'b-', \
					 depth[0,:], X_photon[0,3,:], 'c-')
			plt.legend(('5 keV', '15 keV', '25 keV', '35 keV'), loc='upper right')
		elif use_what == 'nick':
			pass
		else:
			plt.plot(depth, X_photon[0,:], 'r-', \
					 depth, X_photon[1,:], 'g-', \
					 depth, X_photon[2,:], 'b-')
			plt.legend(('5 keV', '10 keV', '20 keV'), loc='upper right')
			
		plt.title("Depth vs. X-ray Photon")
		plt.xlabel("Depth ($\mu m$)")
		plt.ylabel("Number of X-ray photons ($cm^{-3}$)")
		ax = plt.gca()
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(14)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(14)
		#ax.yaxis.set_major_formatter(FuncFormatter(XL.tick_num_format(x,3)))

	if use_what == 'yuan':
		# Lets re-arrange dataset into textfiles. ../XRAY/sorted/
		print "Re-sorting the data to readable forms ...\n"
		for i in range(depth.shape[0]):
			for j in range(X_photon.shape[1]):
				if depth.shape[1] != X_photon.shape[2]:
					raise ValueError, "Wrong Matching between depth and X_photon"
				contents = [] 
			
				# for bottom illumination, we need to arrange the data inversely.
				for k in range(depth.shape[1]-1,-1,-1):
					contents.append(str(depth[i,k])+'\t'+str(X_photon[i,j,k])+'\n')
			
				filename = \
					'./sorted/'+\
					index_map_depth(i)+'nm'+index_map_energy(j)+'.txt'
				fp = open(filename,'w')
				fp.writelines(contents)
				fp.close()
		print "Check ./sorted if the data are stored correctly.\n\n"
	
		print "Preparing to write the Medici codes ..."
	
		print "Generating codes for 10000 um case ..."
		medici_code = []
		for j in range(X_photon.shape[1]):
			input_deck_name = './inp/PHOTOGEN_'+index_map_energy(j)+'.inp' 	
			print "Generating "+input_deck_name+" ... "
			fp = open(input_deck_name,'w')
			light_source_name = '../XRAY/sorted/'+ \
			   index_map_depth(3)+'nm'+index_map_energy(j)+'.txt'
			medici_code = XL.Medici_photogen_codegen(depth[3,:], X_photon[3,j,:], \
					      light_source_name, 0, 45, 1002, 1, 46, 250)
			fp.writelines(medici_code)
			fp.close()
	
	# Using Nick's data to generate the code.
	elif use_what == 'nick':
		pass
	
	# Using Kai's data
	else:
		print "Re-sorting the data to readable forms ...\n"
		print "Applying alternative indices for 5, 10, 20 keV data. \n"
		indices = np.array([0, 0.5, 1.5])
		
		for j in range(X_photon.shape[0]):
			if depth.shape[0] != X_photon.shape[1]:
				raise ValueError, "Wrong Matching between depth and X_photon"
			contents = [] 
						
			# for bottom illumination, we need to arrange the data inversely.
			for k in range(depth.shape[0]-1,-1,-1):
				contents.append(str(depth[k])+'\t'+str(X_photon[j,k])+'\n')
			
			filename = \
				'./sorted/'+\
				'10000'+'um'+index_map_energy(indices[j])+'.txt'
			fp = open(filename,'w')
			fp.writelines(contents)
			fp.close()
		
		print "Check ./sorted if the data are stored correctly.\n\n"

		print "Preparing to write the Medici codes ...\n"

		print "Generating codes for 10000 um case ...\n"
	
		medici_code = []
		for j in range(X_photon.shape[0]):
			input_deck_name = './inp/PHOTOGEN_'+index_map_energy(indices[j])+'.inp' 	
			print "Generating "+input_deck_name+" ... "
			fp = open(input_deck_name,'w')
			light_source_name = '../XRAY/sorted/'+ \
			   '10000'+'um'+index_map_energy(indices[j])+'.txt'
			medici_code = XL.Medici_photogen_codegen(depth[0], X_photon[j,:], \
					      light_source_name, 0, 45, 1000, 1, 46, 250)
			fp.writelines(medici_code)
			fp.close()

	
	
	
	
	print "The code generation finished!\n"
	print "Bye!\n\n"
	
	
	
	
	
	
	
	if debug_mode == 1:
		plt.show()
		
	
	return 0

# Running the main program
if __name__ == "__main__":
    main()


	