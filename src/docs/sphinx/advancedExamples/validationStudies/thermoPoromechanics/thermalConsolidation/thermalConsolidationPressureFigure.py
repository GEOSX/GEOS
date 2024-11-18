import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt

def main():

    # File paths
    hdf5FilePathPressure = "pressureHistory.hdf5"

    
    # Read simulation output from HDF5 file
    hf = h5py.File(hdf5FilePathPressure, 'r')
    timePressure = hf.get('pressure Time')
    timePressure = np.asarray(timePressure)
    centerPressure = hf.get('pressure elementCenter')
    centerPressure = np.asarray(centerPressure)
    pressure = hf.get('pressure')
    pressure = np.asarray(pressure)

    time = 1
    posElement1 = -1
    posElement2 = -1 
    posElement3 = -1
    last = -1
    for j in range(0, centerPressure.shape[1]):
        if centerPressure[0,j,1] >= 0 and posElement1 == -1:
            posElement1 = j
        if centerPressure[0,j,1] >= 4.2 and posElement2 == -1:
            posElement2 = j
        if centerPressure[0,j,1] >= 5.6 and posElement3 == -1:
            posElement3 = j

    for j in range(0, timePressure.shape[0]):            
        if j > 0 and timePressure[j] < 1e-12 and last == -1:     
            last = j

    pressure_0 = genfromtxt('thermoConsolidationPressure_0m.csv', delimiter=',')
    pressure_4p2 = genfromtxt('thermoConsolidationPressure_4p2m.csv', delimiter=',')
    pressure_5p6 = genfromtxt('thermoConsolidationPressure_5p6m.csv', delimiter=',')    

    plt.plot( timePressure[0:last],
              pressure[0:last,posElement1],        
              'r--',
              label='GEOSX: z = 0.0 m')
    plt.plot( timePressure[0:last],
              pressure[0:last,posElement2],        
              'b--',
              label='GEOSX: z = 4.2 m')
    plt.plot( timePressure[0:last],
              pressure[0:last,posElement3],        
              'k--',
              label='GEOSX: z = 5.6 m')

    plt.plot( pressure_0[:,0],
              pressure_0[:,1],        
              'r-',
              label='Analytical: z = 0.0 m')
    plt.plot( pressure_4p2[:,0],
              pressure_4p2[:,1],        
              'b-',
              label='Analytical: z = 4.2 m')
    plt.plot( pressure_5p6[:,0],
              pressure_5p6[:,1],
              'k-',
              label='Analytical: z = 5.6 m')

    plt.xscale("log")
    
    plt.grid()
    plt.ylabel('Pressure [Pa]')
    plt.xlabel('Time [s]')
    plt.xlim(0.01,100000)
    
    plt.legend(loc='lower left')
    plt.show()
    
if __name__ == "__main__":
    main()
