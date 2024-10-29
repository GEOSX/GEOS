import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt

def stressRotation(stress, phi_x):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])
    return np.dot(np.dot(np.transpose(rotx), stress), rotx)

def main():

	bulkModuli = [ 159.4202899e9, 2.298850575e9, 5.535714286e9 ]
	thermalExpansionCoefficients = [ 1.2e-5, 2.0e-5, 2.0e-5 ]

	regions = [ 'casing', 'cement', 'rock' ]

	# Plot GEOSX results
	hasLabel = True

	plt.figure(figsize=(10,5))

	for idx, region in enumerate(regions):
		filePath_stress = 'stressHistory_' + region + '.hdf5'
		filePath_temperature = 'temperatureHistory_' + region + '.hdf5'

		stress_field_name = region + 'Solid_stress'
		
		# Get stress, time and element center
		hf_stress = h5py.File(filePath_stress, 'r')
		time = np.array( hf_stress.get(stress_field_name + ' Time') )
		center = np.array( hf_stress.get(stress_field_name + ' elementCenter') )
		stress = np.array( hf_stress.get(stress_field_name) )

		# Get temperature
		hf_temperature = h5py.File(filePath_temperature, 'r')
		temperature = np.array( hf_temperature.get('temperature') )

		# Compute total stress
		stress_xx_total = stress[:,:,0] - 3 * bulkModuli[idx] * thermalExpansionCoefficients[idx] * temperature
		stress_yy_total = stress[:,:,1] - 3 * bulkModuli[idx] * thermalExpansionCoefficients[idx] * temperature
		stress_zz_total = stress[:,:,2] - 3 * bulkModuli[idx] * thermalExpansionCoefficients[idx] * temperature
		stress_yz_total = stress[:,:,3]
		stress_xz_total = stress[:,:,4]
		stress_xy_total = stress[:,:,5]
		
		# Coordinate of elemnt center
		nElements = center.shape[1]
		xCoord = center[0, :, 0]
		yCoord = center[0, :, 1]
		rCoord = np.sqrt( xCoord*xCoord + yCoord*yCoord )

		# Compute stress components in cylindrical coordinate system
		stress_rr_total = np.zeros(stress_xx_total.shape)
		stress_tt_total = np.zeros(stress_xx_total.shape)

		for idx_time in range(stress.shape[0]):
			for idx_elem in range(stress.shape[1]):
				stressMatrix_cartesian = np.array([[stress_xx_total[idx_time][idx_elem],stress_xy_total[idx_time][idx_elem],stress_xz_total[idx_time][idx_elem]],\
                                 				   [stress_xy_total[idx_time][idx_elem],stress_yy_total[idx_time][idx_elem],stress_yz_total[idx_time][idx_elem]],\
                                 				   [stress_xz_total[idx_time][idx_elem],stress_yz_total[idx_time][idx_elem],stress_zz_total[idx_time][idx_elem]]])

				if(yCoord[idx_elem] != 0):
					phi_x = np.arctan( xCoord[idx_elem]/yCoord[idx_elem] )
				else:
					phi_x = 0

				stressMatrix_cylindirical = stressRotation(stressMatrix_cartesian, phi_x)
				stress_rr_total[idx_time][idx_elem] = stressMatrix_cylindirical[1][1]
				stress_tt_total[idx_time][idx_elem] = stressMatrix_cylindirical[0][0]

				
		# Plot radial stress
		plt.subplot(1,2,1)
		plt.plot( rCoord,
				  stress_rr_total[10, :],        
				  'r+')

		plt.plot( rCoord,
				  stress_rr_total[100, :],        
				  'b+')

		# Plot tangent stress
		plt.subplot(1,2,2)
		if (hasLabel):
			plt.plot( rCoord,
					  stress_tt_total[10, :],        
					  'r+',
					  label='GEOS: t = 1e4 (s)')

			plt.plot( rCoord,
					  stress_tt_total[100, :],        
					  'b+',
					  label='GEOS: t = 1e5 (s)')

			hasLabel = False
		else:
			plt.plot( rCoord,
					  stress_tt_total[10, :],        
					  'r+')

			plt.plot( rCoord,
					  stress_tt_total[100, :],        
					  'b+')
		
	# Reference results
	stress_radial_analytic_1e4s = genfromtxt('stress_radial_analytic_10000s.txt')
	stress_radial_analytic_1e5s = genfromtxt('stress_radial_analytic_100000s.txt')
	stress_tangent_analytic_1e4s = genfromtxt('stress_tangent_analytic_10000s.txt')
	stress_tangent_analytic_1e5s = genfromtxt('stress_tangent_analytic_100000s.txt')

	plt.subplot(1,2,1)
	plt.plot( stress_radial_analytic_1e4s[:,0],
		      stress_radial_analytic_1e4s[:,1],        
		      'r-',
		      label='Analytic: t = 1e4 (s)')

	plt.plot( stress_radial_analytic_1e5s[:,0],
		      stress_radial_analytic_1e5s[:,1],        
		      'b-',
		      label='Analytic: t = 1e5 (s)')

	plt.grid()
	plt.ylabel('Radial stress [Pa]')
	plt.xlabel('Radial coordinate [m]')
	plt.xlim(0.15,0.4)
	plt.ylim(-2.0e6,5.5e6)

	plt.subplot(1,2,2)
	plt.plot( stress_tangent_analytic_1e4s[:,0],
		      stress_tangent_analytic_1e4s[:,1],        
		      'r-',
		      label='Analytic: t = 1e4 (s)')

	plt.plot( stress_tangent_analytic_1e5s[:,0],
		      stress_tangent_analytic_1e5s[:,1],        
		      'b-',
		      label='Analytic: t = 1e5 (s)')

	plt.grid()
	plt.ylabel('Tangent stress [Pa]')
	plt.xlabel('Radial coordinate [m]')
	plt.xlim(0.15,0.4)
	plt.ylim(-10e6,30e6)
	plt.legend(loc='upper right')

	plt.savefig('stress.png')
	
if __name__ == "__main__":
	main()
