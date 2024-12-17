import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
from numpy import genfromtxt
from xml.etree import ElementTree
import os
import argparse

import analyticalResults

def stressRotation(stress, phi_x):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])
    return np.dot(np.dot(np.transpose(rotx), stress), rotx)

def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()
    outputDir = args.outputDir
    geosDir = args.geosDir

    # Material properties from input XML file
    xmlFilePathPrefix = geosDir + "/inputFiles/wellbore/ThermoPoroElasticWellbore"

    xmlData = analyticalResults.getDataFromXML(xmlFilePathPrefix)

    bulkModulus = xmlData[2] # drained bulk modulus of the porous medium
    thermalExpansionCoefficients = xmlData[4] # drained linear thermal expansion coefficient of the porous medium
    grainBulkModulus = xmlData[6]
    BiotCoefficient = 1.0 - bulkModulus / grainBulkModulus

    # Get solid stress, time and element center from GEOSX results
    stress_field_name = 'rockSolid_stress'
    hf_stress = h5py.File(outputDir + '/stressHistory_rock.hdf5', 'r')
    time = np.asarray( hf_stress.get(stress_field_name + ' Time') )
    center = np.asarray( hf_stress.get(stress_field_name + ' elementCenter') )
    stress = np.asarray( hf_stress.get(stress_field_name) )

    # Get temperature and pore pressure for computing the total stress
    hf_temperature = h5py.File(outputDir + '/temperatureHistory_rock.hdf5', 'r')
    temperature = np.asarray( hf_temperature.get('temperature') )

    hf_pressure = h5py.File(outputDir + '/pressureHistory_rock.hdf5', 'r')
    pressure = np.asarray( hf_pressure.get('pressure') )

    # Compute total stress
    stress_xx_total = stress[:,:,0] - BiotCoefficient * pressure - 3 * bulkModulus * thermalExpansionCoefficients * temperature
    stress_yy_total = stress[:,:,1] - BiotCoefficient * pressure - 3 * bulkModulus * thermalExpansionCoefficients * temperature
    stress_zz_total = stress[:,:,2] - BiotCoefficient * pressure - 3 * bulkModulus * thermalExpansionCoefficients * temperature
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

    # Plot GEOSX results    
    plt.figure(figsize=(10,5))

    plt.subplot(1,2,1)

    plt.plot( rCoord,
              stress_rr_total[10, :],        
              'r+',
              label='GEOS: t = 1 (min)')

    plt.plot( rCoord,
              stress_rr_total[24, :],        
              'b+',
              label='GEOS: t = 1 (hour)')
    
    plt.subplot(1,2,2)

    plt.plot( rCoord,
              stress_tt_total[10, :],        
              'r+',
              label='GEOS: t = 1 (min)')

    plt.plot( rCoord,
              stress_tt_total[24, :],        
              'b+',
              label='GEOS: t = 1 (hour)')
    
    # Plot analytical results at one minute
    t = 60 #s
    r_anal, T_anal, p_anal, ur_anal, sig_rr_anal, sig_tt_anal, sig_zz_anal = analyticalResults.analyticalResults(t)

    plt.subplot(1,2,1)
    plt.plot( r_anal,
              sig_rr_anal,        
              'r-',
              label='Analytical: t = 1 (min)')

    plt.subplot(1,2,2)
    plt.plot( r_anal,
              sig_tt_anal,        
              'r-',
              label='Analytical: t = 1 (min)')

    # Plot analytical results at one hour
    t = 3600 #s
    r_anal, T_anal, p_anal, ur_anal, sig_rr_anal, sig_tt_anal, sig_zz_anal = analyticalResults.analyticalResults(t)

    plt.subplot(1,2,1)
    plt.plot( r_anal,
              sig_rr_anal,        
              'b-',
              label='Analytical: t = 1 (hour)')

    plt.subplot(1,2,2)
    plt.plot( r_anal,
              sig_tt_anal,        
              'b-',
              label='Analytical: t = 1 (hour)')

    plt.subplot(1,2,1)
    plt.grid()
    plt.ylabel('Radial stress [Pa]')
    plt.xlabel('Radial coordinate [m]')
    plt.xlim(0.1,0.3)
    #plt.ylim(0,100)

    plt.subplot(1,2,2)
    plt.grid()
    plt.ylabel('Tangential stress [Pa]')
    plt.xlabel('Radial coordinate [m]')
    plt.xlim(0.1,0.3)
    #plt.ylim(0,70e6)
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig('stress.png')
    plt.isoff()

if __name__ == "__main__":
    main()
    
