import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
import analyticalResults

def main():

	# File paths
	hdf5FilePath = "displacementJumpHistory.hdf5"

	# Plot GEOSX results
	hf = h5py.File(hdf5FilePath, 'r')
	time = np.asarray( hf.get('displacementJump Time') )
	displacementJump = np.asarray( hf.get('displacementJump') )

	nTime = time.shape[0]

	displacementJump_normal = displacementJump[0:nTime, 0, 0]	
	
	# Reference results
	listTime = np.linspace(0, 0.5, 21, endpoint=False)
	listP0 = np.linspace(0, -10e6, 21, endpoint=False)
	listTime2 = np.linspace(0.5, 2.0, 21, endpoint=True)
	listP02 = np.linspace(-10e6, 10e6, 21, endpoint=True)
	for x in listTime2:
		listTime = np.append(listTime, x)
	for x in listP02:
		listP0 = np.append(listP0, x)

	list_displacementJump_anal = []
	for P0 in listP0:
		displacementJump_anal = analyticalResults.analyticDisplacementJump(P0) 
		list_displacementJump_anal.append(displacementJump_anal * 1e6) # converted to um

	# Plot displacement jump versus time
	plt.plot( time,
			  displacementJump_normal * 1e6,        # converted to um
			  'k+',
			  label='GEOS')

	plt.plot( listTime,
			  list_displacementJump_anal,
			  'r-',
			  label='Analytical')
	
	plt.grid()
	plt.ylabel(r'Normal displacement jump [$\mu$m]')
	plt.xlabel('Time [s]')
	plt.xlim(0, 2)

	plt.legend(loc='upper left')
	plt.show()

if __name__ == "__main__":
	main()
