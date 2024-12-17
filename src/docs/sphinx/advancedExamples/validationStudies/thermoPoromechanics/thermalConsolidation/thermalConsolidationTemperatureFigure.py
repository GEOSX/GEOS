import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt

def main():

    # File paths
    hdf5FilePathTemperature = "temperatureHistory.hdf5"

    
    # Read simulation output from HDF5 file
    hf = h5py.File(hdf5FilePathTemperature, 'r')
    timeTemperature = hf.get('temperature Time')
    timeTemperature = np.asarray(timeTemperature)
    centerTemperature = hf.get('temperature elementCenter')
    centerTemperature = np.asarray(centerTemperature)
    temperature = hf.get('temperature')
    temperature = np.asarray(temperature)
    
    time = 1
    posElement1 = -1
    posElement2 = -1 
    posElement3 = -1
    last = -1
    
    for j in range(0, centerTemperature.shape[1]):
        if centerTemperature[0,j,1] >= 0 and posElement1 == -1:
            posElement1 = j
        if centerTemperature[0,j,1] >= 4.2 and posElement2 == -1:
            posElement2 = j
        if centerTemperature[0,j,1] >= 5.6 and posElement3 == -1:
            posElement3 = j

    for j in range(0, timeTemperature.shape[0]):            
        if j > 0 and timeTemperature[j] < 1e-12 and last == -1:     
            last = j
            
    temperature_0 = genfromtxt('thermoConsolidationTemp_0m.csv', delimiter=',')
    temperature_4p2 = genfromtxt('thermoConsolidationTemp_4p2m.csv', delimiter=',')
    temperature_5p6 = genfromtxt('thermoConsolidationTemp_5p6m.csv', delimiter=',')    

    plt.plot( timeTemperature[0:last],
              temperature[0:last,posElement1],        
              'r--',
              label='GEOSX: z = 0.0 m')
    plt.plot( timeTemperature[0:last],
              temperature[0:last,posElement2],        
              'b--',
              label='GEOSX: z = 4.2 m')
    plt.plot( timeTemperature[0:last],
              temperature[0:last,posElement3],        
              'k--',
              label='GEOSX: z = 5.6 m')

    plt.plot( temperature_0[:,0],
              273+temperature_0[:,1],        
              'r-',
              label='Analytical: z = 0.0 m')
    plt.plot( temperature_4p2[:,0],
              273+temperature_4p2[:,1],        
              'b-',
              label='Analytical: z = 4.2 m')
    plt.plot( temperature_5p6[:,0],
              273+temperature_5p6[:,1],
              'k-',
              label='Analytical: z = 5.6 m')

    plt.xscale("log")
    
    plt.grid()
    plt.ylabel('Temperature [K]')
    plt.xlabel('Time [s]')
    plt.ylim(273,323)
    plt.xlim(0.01,100000)    
    
    plt.legend(loc='upper left')
    plt.show()
    #plt.savefig(' temperature.png')
    
if __name__ == "__main__":
    main()
