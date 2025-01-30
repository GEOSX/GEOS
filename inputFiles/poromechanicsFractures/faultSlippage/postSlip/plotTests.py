import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
import math
from math import sin, cos, tan, exp, atan, asin


def main():
    # File path
    hdf5File1Path = "traction_history.hdf5"
    hdf5File2Path = "displacementJump_history.hdf5"    

    # Read HDF5
    # Global Coordinate of Fracture Element Center
    hf = h5py.File(hdf5File1Path, 'r')
    xl = hf.get('traction elementCenter')
    xl = np.array(xl)
    xcord = xl[0, :, 0]
    ycord = xl[0, :, 1]
    zcord = xl[0, :, 2]

    # Local Normal Traction
    hf = h5py.File(hdf5File1Path, 'r')
    trac = hf.get('traction')
    trac = np.array(trac)
    normalTraction = trac[-1, :, 0]
    shearTraction = trac[-1, :, 2]

    # Local Shear Displacement
    hf = h5py.File(hdf5File2Path, 'r')
    jump = hf.get('displacementJump')
    jump = np.array(jump)
    normalJump = jump[-1, :, 0]
    shearJump = jump[-1, :, 2]

    depth = []
    xloc = []
    tnlist = []
    tslist = []
    dnlist = []
    dslist = []
    for i in range(0, len(zcord)):
        if abs(zcord[i] / 12.5 - 1.) < 0.01:            
            depth.append(ycord[i])
            tnlist.append(normalTraction[i] / 1.0e6)
            tslist.append(shearTraction[i] / 1.0e6)
            dnlist.append(normalJump[i])
            dslist.append(shearJump[i])

    # Load reference data
    depth_ref, ts_ref = np.loadtxt('Sxy_Analytical.txt', skiprows=1, unpack=True)
    depth_ref2, ds_ref = np.loadtxt('slip_Analytical.txt', skiprows=1, unpack=True)    


    fsize = 30
    msize = 12
    lw = 8
    malpha = 0.6
    lalpha = 0.6
    fig, ax = plt.subplots(1, 2, figsize=(24, 16))
    cmap = plt.get_cmap("tab10")

    #ax[0].plot(tnlist, depth, color=cmap(0), label='Normal', lw=lw)
    ax[0].plot(tslist, depth, 'o', color=cmap(2), markersize=msize, alpha=malpha, label='GEOS', mec='k', mew=2)
    ax[0].plot(ts_ref/1.0e6, depth_ref, color=cmap(2), alpha=lalpha, lw=lw, label='Analytical')
    ax[0].grid()
    ax[0].set_xlim(-20, 20)
    ax[0].set_ylim(-250, 250)
    ax[0].set_xlabel('Shear Stress [MPa]', size=fsize, weight="bold")
    ax[0].set_ylabel('Depth [m]', size=fsize, weight="bold")
    ax[0].legend(loc='lower right', fontsize=fsize)
    ax[0].xaxis.set_tick_params(labelsize=fsize)
    ax[0].yaxis.set_tick_params(labelsize=fsize)

    #ax[1].plot(dnlist, depth, color=cmap(0), label='Normal', lw=lw)    
    ax[1].plot(dslist, depth, 'o', color=cmap(2), markersize=msize, alpha=malpha, label='GEOS', mec='k', mew=2)
    ax[1].plot(ds_ref, depth_ref2, color=cmap(2), alpha=lalpha, lw=lw, label='Analytical')    
    ax[1].grid()
    #ax[1].set_xlim(-1, 1)
    ax[1].set_ylim(-250, 250)
    ax[1].set_xlabel('DisplacementJump [m]', size=fsize, weight="bold")
    ax[1].set_ylabel('Depth [m]', size=fsize, weight="bold")
    ax[1].legend(loc='lower right', fontsize=fsize)
    ax[1].xaxis.set_tick_params(labelsize=fsize)
    ax[1].yaxis.set_tick_params(labelsize=fsize)

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

    fig.savefig('Simulation_Results.png')

    plt.show()


if __name__ == "__main__":
    main()
