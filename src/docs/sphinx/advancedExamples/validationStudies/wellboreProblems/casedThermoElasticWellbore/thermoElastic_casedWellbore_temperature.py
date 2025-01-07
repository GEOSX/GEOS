import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt

def main():

	# File paths
	hdf5FilePath = [ "temperatureHistory_casing.hdf5",
		             "temperatureHistory_cement.hdf5",
		             "temperatureHistory_rock.hdf5" ]

	# Plot GEOS results
	hasLabel = True

	for filePath in hdf5FilePath:
		hf = h5py.File(filePath, 'r')
		time = np.asarray( hf.get('temperature Time') )
		center = np.asarray( hf.get('temperature elementCenter') )
		temperature = np.asarray( hf.get('temperature') )

		nElements = center.shape[1]
		xCoord = center[0, 0:nElements, 0]
		yCoord = center[0, 0:nElements, 1]
		rCoord = np.sqrt( xCoord*xCoord + yCoord*yCoord )

		if (hasLabel):
			plt.plot( rCoord,
					  temperature[9, 0:nElements],        
					  'r+',
					  label='GEOS: t = 1e4 (s)')

			plt.plot( rCoord,
					  temperature[99, 0:nElements],        
					  'b+',
					  label='GEOS: t = 1e5 (s)')

			hasLabel = False
		else:
			plt.plot( rCoord,
					  temperature[9, 0:nElements],        
					  'r+')

			plt.plot( rCoord,
					  temperature[99, 0:nElements],        
					  'b+')

	temperature_radial_analytic_1e4s = genfromtxt('temperature_radial_analytic_10000s.txt')
	temperature_radial_analytic_1e5s = genfromtxt('temperature_radial_analytic_100000s.txt')

	plt.plot( temperature_radial_analytic_1e4s[:,0],
		      temperature_radial_analytic_1e4s[:,1],        
		      'r-',
		      label='Analytic: t = 1e4 (s)')

	plt.plot( temperature_radial_analytic_1e5s[:,0],
		      temperature_radial_analytic_1e5s[:,1],        
		      'b-',
		      label='Analytic: t = 1e5 (s)')

	plt.grid()
	plt.ylabel('Temperature [°C]')
	plt.xlabel('Radial coordinate [m]')
	plt.xlim(0.15,0.4)

	plt.legend(loc='lower right')
	plt.savefig('temperature.png')

if __name__ == "__main__":
	main()
