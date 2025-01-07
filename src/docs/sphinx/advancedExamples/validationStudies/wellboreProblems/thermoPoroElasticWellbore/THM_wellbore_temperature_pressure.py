import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
from numpy import genfromtxt
import analyticalResults

def main():

	# Plot GEOSX results
	hf = h5py.File("temperatureHistory_rock.hdf5", 'r')
	time = np.asarray( hf.get('temperature Time') )
	center = np.asarray( hf.get('temperature elementCenter') )
	temperature = np.asarray( hf.get('temperature') )

	hf = h5py.File("pressureHistory_rock.hdf5", 'r')
	pressure = np.asarray( hf.get('pressure') )
	
	nElements = center.shape[1]
	xCoord = center[0, 0:nElements, 0]
	yCoord = center[0, 0:nElements, 1]
	rCoord = np.sqrt( xCoord*xCoord + yCoord*yCoord )

	plt.figure(figsize=(10,5))

	plt.subplot(1,2,1)

	plt.plot( rCoord,
			  temperature[10, 0:nElements], # index of the time step at 60s is 10       
			  'r+',
			  label='GEOS: t = 1 (min)')

	plt.plot( rCoord,
			  temperature[24, 0:nElements], # index of the time step at 3600s is 24       
			  'b+',
			  label='GEOS: t = 1 (hour)')
	
	plt.subplot(1,2,2)

	plt.plot( rCoord,
			  pressure[10, 0:nElements],        
			  'r+')

	plt.plot( rCoord,
			  pressure[24, 0:nElements],        
			  'b+')
	
	# Plot analytical results at one minute
	t = 60 #s
	r_anal, T_anal, p_anal, ur_anal, sig_rr_anal, sig_tt_anal, sig_zz_anal = analyticalResults.analyticalResults(t)

	plt.subplot(1,2,1)
	plt.plot( r_anal,
			  T_anal,        
			  'r-',
			  label='Analytical: t = 1 (min)')

	plt.subplot(1,2,2)
	plt.plot( r_anal,
			  p_anal,        
			  'r-')

	# Plot analytical results at one hour
	t = 3600 #s
	r_anal, T_anal, p_anal, ur_anal, sig_rr_anal, sig_tt_anal, sig_zz_anal = analyticalResults.analyticalResults(t)

	plt.subplot(1,2,1)
	plt.plot( r_anal,
			  T_anal,        
			  'b-',
			  label='Analytical: t = 1 (hour)')

	plt.subplot(1,2,2)
	plt.plot( r_anal,
			  p_anal,        
			  'b-')

	plt.subplot(1,2,1)
	plt.grid()
	plt.ylabel('Temperature [°C]')
	plt.xlabel('Radial coordinate [m]')
	plt.xlim(0.1,0.3)
	plt.ylim(0,100)

	plt.legend(loc='upper right')

	plt.subplot(1,2,2)
	plt.grid()
	plt.ylabel('Pore pressure [Pa]')
	plt.xlabel('Radial coordinate [m]')
	plt.xlim(0.1,0.3)
	#plt.ylim(0,70e6)
	plt.tight_layout()
	plt.savefig('temperature_pressure.png')
	plt.isoff()

if __name__ == "__main__":
	main()
	
