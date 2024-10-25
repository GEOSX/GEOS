import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
from numpy import genfromtxt

def main():

    # File paths
    hdf5FilePathDisplacement = "displacementHistory.hdf5"

    
    # Read simulation output from HDF5 file
    hf = h5py.File(hdf5FilePathDisplacement, 'r')
    timeDisplacement = hf.get('totalDisplacement Time')
    timeDisplacement = np.asarray(timeDisplacement)
    centerDisplacement = hf.get('totalDisplacement ReferencePosition')
    centerDisplacement = np.asarray(centerDisplacement)
    displacement = hf.get('totalDisplacement')
    displacement = np.asarray(displacement)

    time = 1
    posVertex1 = -1
    posVertex2 = -1 
    posVertex3 = -1
    last = -1
    
    for j in range(0, 284):            
        if centerDisplacement[0,j,1] >= 1.4 and posVertex1 == -1:
            posVertex1 = j
        if centerDisplacement[0,j,1] >= 4.2 and posVertex2 == -1:
            posVertex2 = j
        if centerDisplacement[0,j,1] >= 7 and posVertex3 == -1:
            posVertex3 = j

    for j in range(0, timeDisplacement.shape[0]):            
        if j > 0 and timeDisplacement[j] < 1e-12 and last == -1:     
            last = j

    displacement_1p4 = genfromtxt('thermoConsolidationDisp_1p4m.csv', delimiter=',')
    displacement_4p2 = genfromtxt('thermoConsolidationDisp_4p2m.csv', delimiter=',')
    displacement_7 = genfromtxt('thermoConsolidationDisp_7m.csv', delimiter=',')    

    plt.plot( timeDisplacement[0:last],
              -displacement[0:last,posVertex1,1],        
              'r--',
              label='GEOSX: z = 1.4 m')
    plt.plot( timeDisplacement[0:last],
              -displacement[0:last,posVertex2,1],
              'b--',
              label='GEOSX: z = 4.2 m')
    plt.plot( timeDisplacement[0:last],
              -displacement[0:last,posVertex3,1],        
              'k--',
              label='GEOSX: z = 5.6 m')

    plt.plot( displacement_1p4[:,0],
              displacement_1p4[:,1],        
              'r-',
              label='Analytical: z = 1.4 m')
    plt.plot( displacement_4p2[:,0],
              displacement_4p2[:,1],        
              'b-',
              label='Analytical: z = 4.2 m')
    plt.plot( displacement_7[:,0],
              displacement_7[:,1],
              'k-',
              label='Analytical: z = 7 m')
                           
    plt.xscale("log")
    
    plt.grid()
    plt.ylabel('Displacement [m]')
    plt.xlabel('Time [s]')
    plt.xlim(0.01,100000)
    
    plt.legend(loc='upper left')
    plt.show()
    # plt.savefig('displacement.png')


    
if __name__ == "__main__":
    main()
