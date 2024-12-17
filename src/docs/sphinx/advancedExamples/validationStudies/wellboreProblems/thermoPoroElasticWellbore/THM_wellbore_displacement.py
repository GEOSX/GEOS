import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
from numpy import genfromtxt
import analyticalResults

def main():
	# Plot GEOSX results
	hf = h5py.File("displacementHistory.hdf5", 'r')
	time = np.asarray( hf.get('totalDisplacement Time') )
	center = np.asarray( hf.get('totalDisplacement ReferencePosition') )
	displacement = np.asarray( hf.get('totalDisplacement') )

	nNodes = center.shape[1]
	xCoord = center[0, 0:nNodes, 0]
	yCoord = center[0, 0:nNodes, 1]	
	ux_60s = displacement[10, 0:nNodes, 0] # the index for time step at 60s is 10
	ux_3600s = displacement[24, 0:nNodes, 0] # the index for time step at 3600s is 24

	rCoord, ur_60s, ur_3600s = [], [], []

	# Extract radial data at theta = 0
	for idx in range(nNodes):
		if (yCoord[idx] < 1e-6):
			rCoord.append(xCoord[idx])
			ur_60s.append(ux_60s[idx]*1e6) # converted to um
			ur_3600s.append(ux_3600s[idx]*1e6) # converted to um

	plt.figure(figsize=(7,5))

	plt.plot( rCoord,
			  ur_60s,        
			  'r+',
			  label='GEOS: t = 1 (min)')

	plt.plot( rCoord,
			  ur_3600s,        
			  'b+',
			  label='GEOS: t = 1 (hour)')
	
	# Plot analytical results at one minute
	t = 60 #s
	r_anal, T_anal, p_anal, ur_anal, sig_rr_anal, sig_tt_anal, sig_zz_anal = analyticalResults.analyticalResults(t)

	plt.plot( r_anal,
			  ur_anal * 1e6, # converted to um        
			  'r-',
			  label='Analytical: t = 1 (min)')

	# Plot analytical results at one hour
	t = 3600 #s
	r_anal, T_anal, p_anal, ur_anal, sig_rr_anal, sig_tt_anal, sig_zz_anal = analyticalResults.analyticalResults(t)

	plt.plot( r_anal,
			  ur_anal * 1e6, # converted to um        
			  'b-',
			  label='Analytical: t = 1 (hour)')

	plt.grid()
	plt.ylabel(r'Displacement [$\mu$m]')
	plt.xlabel('Radial coordinate [m]')
	plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
	plt.xlim(0.1,0.3)
	plt.ylim(0,500)
	plt.legend(loc='upper left')
	plt.savefig('displacement.png')
	plt.isoff()
	
if __name__ == "__main__":
	main()
	
