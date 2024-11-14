import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
from mpmath import *
import math
import argparse
import os


def main(directory):
    # File paths
    filename = "temperature_history.hdf5"
    hdf5FilePathTemperature = os.path.join( directory, filename )    
    
    # Read simulation output from HDF5 file
    hf = h5py.File(hdf5FilePathTemperature, 'r')
    timeTemperature = hf.get('temperature Time')
    timeTemperature = np.asarray(timeTemperature)
    centerTemperature = hf.get('temperature elementCenter')
    centerTemperature = np.asarray(centerTemperature)
    xcord_elm = centerTemperature[0,:,0]
    ycord_elm = centerTemperature[0,:,1]
    zcord_elm = centerTemperature[0,:,2]
    temperature = hf.get('temperature')
    temperature = np.asarray(temperature)
    print(ycord_elm)
    
    filename = "pressure_history.hdf5"
    hdf5FilePathPressure = os.path.join( directory, filename ) 

    hf = h5py.File(hdf5FilePathPressure, 'r')
    pressure = np.asarray( hf.get('pressure') )
    
    t_index = [0, 1, 2, 9]
    tstar = [0, 100, 200, 900]     
    numFiles=len(tstar)
    # Extract Curve
    midY = 10
    xlist = []
    pplist = []
    temperaturelist = []
    for tt in range(len(t_index)):
        xtemp = []
        ptemp = []
        ttemp = []
        for i in range(0,len(zcord_elm)):
            if abs(ycord_elm[i] - midY) < 0.01:
               print(i)
               print(xcord_elm[i], pressure[tt,i]/1.0e6, temperature[tt,i])
               xtemp.append(xcord_elm[i])
               ptemp.append(pressure[tt,i]/1.0e6)
               ttemp.append(temperature[tt,i])

        xlist.append(xtemp)
        pplist.append(ptemp)
        temperaturelist.append(ttemp)

    xlist = np.asarray(xlist)
    pplist = np.asarray(pplist)
    temperaturelist = np.asarray(temperaturelist)


    #Visualization
    N1 = 1
    fsize = 32
    msize = 10
    lw = 4
    mew = 2
    malpha = 0.6
    lalpha = 0.8

    fig, ax = plt.subplots(2, 2, figsize=(32, 18))
    colorlist = ['lime', 'orangered', 'black', 'deepskyblue']
    markerlist = ['o', 's', 'P', 'D']

    for i in range(0,numFiles):
        ax[0,0].plot(xlist[i][0::N1], pplist[i][0::N1], linestyle='None', marker=markerlist[i], alpha=0.4, markerfacecolor=colorlist[i], fillstyle='full', markersize=msize, label='GEOS_t='+str(tstar[i])+'s')
    #ax[0,0].set_xlim([1, 10])
    #ax[0,0].set_ylim([-1, 11])
    ax[0,0].set_xlabel(r'r_d', size=fsize, weight="bold")
    ax[0,0].set_ylabel(r'Pore Pressure (MPa)', size=fsize, weight="bold")
    ax[0,0].legend(loc='upper right',fontsize=fsize*0.8)
    ax[0,0].grid(True, which="both", ls="-")
    ax[0,0].xaxis.set_tick_params(labelsize=fsize)
    ax[0,0].yaxis.set_tick_params(labelsize=fsize)


    for i in range(0,numFiles):
        ax[1,0].plot(xlist[i][0::N1], temperaturelist[i][0::N1], linestyle='None', marker=markerlist[i], alpha=0.4, markerfacecolor=colorlist[i], fillstyle='full', markersize=msize, label='GEOS_t='+str(tstar[i])+'s')
    #ax[1,0].set_xlim([1, 10])
    #ax[1,0].set_ylim([-1, 11])
    ax[1,0].set_xlabel(r'r_d', size=fsize, weight="bold")
    ax[1,0].set_ylabel(r'Temperature (K)', size=fsize, weight="bold")
    ax[1,0].legend(loc='upper right',fontsize=fsize*0.8)
    ax[1,0].grid(True, which="both", ls="-")
    ax[1,0].xaxis.set_tick_params(labelsize=fsize)
    ax[1,0].yaxis.set_tick_params(labelsize=fsize)
    
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

    fig.savefig('Verification_tutorial.png')

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot solution of fractureMatrixThermalFlow")
    parser.add_argument("-d", "--workingDir", type=str, help="Directory containing the solution")

    args, unknown_args = parser.parse_known_args()
    if unknown_args:
        print("unknown arguments %s" % unknown_args)
    
    directory = args.directory

    main( directory )
