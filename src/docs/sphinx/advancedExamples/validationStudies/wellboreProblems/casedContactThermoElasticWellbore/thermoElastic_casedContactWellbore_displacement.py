import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt

def main():

	# File paths
	hdf5FilePath = "displacementHistory.hdf5"

	# Plot GEOSX results
	hf = h5py.File(hdf5FilePath, 'r')
	time = np.asarray( hf.get('totalDisplacement Time') )
	center = np.asarray( hf.get('totalDisplacement ReferencePosition') )
	displacement = np.asarray( hf.get('totalDisplacement') )
	
	nNodes = center.shape[1]
	xCoord = center[0, 0:nNodes, 0]
	yCoord = center[0, 0:nNodes, 1]
	
	ux_10000s = displacement[10, 0:nNodes, 0] # index 10 is for 1e4s
	uy_10000s = displacement[10, 0:nNodes, 1]
	
	ux_100000s = displacement[100, 0:nNodes, 0] # index 100 is for 1e5s
	uy_100000s = displacement[100, 0:nNodes, 1]

	rCoord, ur_10000s, ur_100000s = [], [], []

	# Extract radial data at theta = 0
	for idx in range(nNodes):
		if (yCoord[idx] < 1e-6):
			rCoord.append(xCoord[idx])
			ur_10000s.append(ux_10000s[idx])
			ur_100000s.append(ux_100000s[idx])
	
	# Reference results
	displacement_radial_analytic_1e4s = genfromtxt('displacement_radial_analytic_10000s.txt')
	displacement_radial_analytic_1e5s = genfromtxt('displacement_radial_analytic_100000s.txt')

	# Plot radial displacement at 1e4 (s)
	plt.plot( rCoord,
			  ur_10000s,        
			  'r+',
			  label='GEOS: t = 1e4 (s)')
	
	plt.plot( displacement_radial_analytic_1e4s[:,0],
		      displacement_radial_analytic_1e4s[:,1],
		      'r-',
		      label='Analytic: t = 1e4 (s)')
	
	# Plot radial displacement at 1e5 (s)
	plt.plot( rCoord,
			  ur_100000s,        
			  'b+',
			  label='GEOS: t = 1e5 (s)')
	
	plt.plot( displacement_radial_analytic_1e5s[:,0],
		      displacement_radial_analytic_1e5s[:,1],
		      'b-',
		      label='Analytic: t = 1e5 (s)')
	
	plt.grid()
	plt.ylabel(r'Displacement [m]')
	plt.xlabel('Radial coordinate [m]')
	plt.xlim(0.15,0.4)
	plt.ylim(-600.0e-6,100.0e-6)

	plt.legend(loc='lower right')

	plt.savefig('displacement.png')

if __name__ == "__main__":
	main()
