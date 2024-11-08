import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt
import os
import argparse

# This function allows rotating the stress matrix (with Voigt notation) around an angle on the plane x-y
def stressRotation(stress, phi_x):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])
    return np.dot(np.dot(np.transpose(rotx), stress), rotx)

def getHydromechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)
    param1 = tree.find('Constitutive/ElasticIsotropic')
   
    hydromechanicalParameters = dict.fromkeys([
        "bulkModulus", "shearModulus", "defaultDrainedTEC", "dDrainedTEC_dT", "referenceTemperature"
    ])

    hydromechanicalParameters["bulkModulus"] = float(param1.get("defaultBulkModulus"))
    hydromechanicalParameters["shearModulus"] = float(param1.get("defaultShearModulus"))   
    hydromechanicalParameters["defaultDrainedTEC"] = float(param1.get("defaultDrainedTEC"))
    hydromechanicalParameters["dDrainedTEC_dT"] = float(param1.get("dDrainedTEC_dT"))
    hydromechanicalParameters["referenceTemperature"] = float(param1.get("referenceTemperature"))

    return hydromechanicalParameters

def main():
   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # File path
    outputDir = args.outputDir
    geosDir = args.geosDir
    hdf5File1Path = outputDir + "/stressHistory_rock.hdf5"
    hdf5File2Path = outputDir + "/temperatureHistory_rock.hdf5"
    xmlFilePath = geosDir + "/inputFiles/thermoPoromechanics/ThermoElasticOedometric_linearTEC_base.xml"

    # Extract info from XML
    hydromechanicalParameters = getHydromechanicalParametersFromXML(xmlFilePath)

    bulkModulus = hydromechanicalParameters["bulkModulus"]
    shearModulus = hydromechanicalParameters["shearModulus"]
    poissonRatio = (3.0*bulkModulus - 2.0*shearModulus) / (6.0*bulkModulus + 2.0*shearModulus)
    youngModulus = 3.0 * bulkModulus * (1.0 - 2.0*poissonRatio)

    defaultDrainedTEC = hydromechanicalParameters["defaultDrainedTEC"]
    dDrainedTEC_dT = hydromechanicalParameters["dDrainedTEC_dT"]
    referenceTemperature = hydromechanicalParameters["referenceTemperature"]

    # Extract stress
    hf_stress = h5py.File(hdf5File1Path, 'r')
    time = np.array( hf_stress.get('rockSolid_linearTEC_stress Time') )
    center = np.array( hf_stress.get('rockSolid_linearTEC_stress elementCenter') )
    stress = np.array( hf_stress.get('rockSolid_linearTEC_stress') )

    hf_temperature = h5py.File(hdf5File2Path, 'r')
    temperature = np.array( hf_temperature.get('temperature') )

    nElements = center.shape[1]
    nTimes = time.shape[0]   

    element_idx = 0
    temp = temperature[:,element_idx]

    # Get stress: the contribution of temperature to stress is already computed within this GEOS version.
    stress_xx_total = stress[:,element_idx,0] 
    stress_yy_total = stress[:,element_idx,1]
    stress_zz_total = stress[:,element_idx,2] 
    stress_yz_total = stress[:,element_idx,3]
    stress_xz_total = stress[:,element_idx,4]
    stress_xy_total = stress[:,element_idx,5]

    # Analytical results
    drainedTEC = defaultDrainedTEC +  dDrainedTEC_dT*(temp - referenceTemperature)
    pconf_anal = youngModulus/(1.0-poissonRatio) * (defaultDrainedTEC + drainedTEC)/2.0 * (temp - referenceTemperature)

    # Plot 
    plt.plot( temp,
              -stress_xx_total * 1e-5, # convert to bar    
              'k+',
              label=f'GEOS with linear TEC')

    plt.plot( temp,
              pconf_anal * 1e-5, # convert to bar  
              'b-',
              label='Analytical results')
	
    plt.grid()
    plt.ylabel('Pconf [Bar]')
    plt.xlabel('Temperature [C]')
    plt.xlim(20.0, 90.0)
    plt.ylim(0.0, 200.0)   

    plt.legend(loc='lower right')
    plt.show()
	
if __name__ == "__main__":
    main()
