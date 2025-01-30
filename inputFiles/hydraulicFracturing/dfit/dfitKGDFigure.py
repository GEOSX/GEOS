import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin,cos,tan,exp
import h5py
import xml.etree.ElementTree as ElementTree


def getMeshSettings(xmlFilePath):
    
    tree = ElementTree.parse(xmlFilePath)
    meshElement = tree.find('Mesh/InternalMesh')

    nodeXCoords = meshElement.get("xCoords")
    nodeYCoords = meshElement.get("yCoords")
    nodeZCoords = meshElement.get("zCoords")
    
    nodeXCoords = [float(i) for i in nodeXCoords[1:-1].split(",")]
    nodeYCoords = [float(i) for i in nodeYCoords[1:-1].split(",")]
    nodeZCoords = [float(i) for i in nodeZCoords[1:-1].split(",")]    

    xMin = nodeXCoords[0]
    xMax = nodeXCoords[-1]
    yMin = nodeYCoords[0]
    yMax = nodeYCoords[-1]
    zMin = nodeZCoords[0]
    zMax = nodeZCoords[-1]

    nXElem = meshElement.get("nx")
    nYElem = meshElement.get("ny")
    nZElem = meshElement.get("nz")
    nXElem = [float(i) for i in nXElem[1:-1].split(",")]
    nYElem = [float(i) for i in nYElem[1:-1].split(",")]
    nZElem = [float(i) for i in nZElem[1:-1].split(",")]  
    nx = nXElem[0]
    ny = nYElem[0]
    nz = nZElem[0]
    
    return xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz

def Gfunction(delt_D):
    gfun = 4./3.*( (1.+ delt_D)**1.5 - delt_D**1.5 )
    g0 = 4./3.
    Gfun = 4./math.pi*(gfun - g0 )
    return Gfun

def main():
    # Pumping schedule
    q = 3.25 # kg/s
    t_shutDown = 10.0
    epsilon = 0.0001
    Rate = [0.0, q, q, 0.0, 0.0 ]
    Rate_time = [0.0, 0.0+ epsilon, t_shutDown - epsilon,  t_shutDown + epsilon, 300 ]

    # Load and process GEOSX results
    # File path
    hdf5File1Path = "/usr/workspace/cusini1/geosx/geosx_dev/GEOSX/Output/dfit_kgd_output.hdf5"
    xmlFile1Path = "/usr/workspace/cusini1/geosx/geosx_dev/GEOSX/inputFiles/hydraulicFracturing/dfit/dfit_kgd_base.xml"
    xmlFile2Path = "/usr/workspace/cusini1/geosx/geosx_dev/GEOSX/inputFiles/hydraulicFracturing/dfit/dfit_kgd_benchmark.xml"

    # Read simulation parameters from XML file
    xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz = getMeshSettings(xmlFile2Path)
    fracHeight = abs(zMax-zMin)
    dy = abs(yMax-yMin)/ny
    dz = abs(zMax-zMin)/nz 

    # Read simulation output from HDF5 file
    # Global Coordinate of Element Center
    hf = h5py.File(hdf5File1Path, 'r')
    xl = hf.get('pressure elementCenter')
    xl = np.array(xl)
    xcord = xl[0,:,0]
    ycord = xl[0,:,1]
    zcord = xl[0,:,2]
    tl = hf.get('pressure Time')
    tl = np.array(tl)
    # Load pressure
    fpre = hf.get('pressure')  
    fpre = np.array(fpre)
    # Load elementAperture
    aper = hf.get('elementAperture')  
    aper = np.array(aper)
    # Load elementArea
    area = hf.get('elementArea')  
    area = np.array(area)

    # Find injection location
    for i in range(0,len(zcord)):
        if abs(xcord[i]) < 0.01 and abs(ycord[i]/(dy/2.) - 1) < 0.01 and abs(zcord[i]/(dz/2.) - 1) < 0.01:
           well_ind = i

    # Find shut-in time
    for i in range(0,len(tl)):
        if abs(tl[i]/t_shutDown - 1.0) < 0.001:
           t_ind = i

    # Compute frac length
    length = np.zeros([len(tl)])     
    for i in range(0,len(tl)):
        temp = 0
        for j in range(0,len(aper[0,:])):
            if aper[i,j]>1.0e-5:
               temp += area[i,j]

        length[i] = temp/fracHeight

    # Compute G-function
    dtD = []
    G_tD = []
    P_tD = []
    for i in range(t_ind,len(tl)):
        tD = (tl[i] - tl[t_ind])/tl[t_ind]
        dtD.append(tD)
        G_tD.append(Gfunction(tD))
        P_tD.append(fpre[i,well_ind])

    dtD = np.array(dtD)
    G_tD = np.array(G_tD)
    P_tD = np.array(P_tD)

    dPdG = np.zeros(len(G_tD)-1)
    GdPdG = np.zeros(len(G_tD)-1)
    for i in range(1,len(dPdG)):
        dPdG[i] = -(P_tD[i+1] - P_tD[i-1])/(G_tD[i+1] - G_tD[i-1])
        GdPdG[i] = G_tD[i] * dPdG[i]

    # Search the location of the tangent point 
    slope = np.zeros(len(GdPdG))
    mslop = 0
    for i in range(1,len(slope)):
        slope[i] = (GdPdG[i]-GdPdG[0])/(G_tD[i] - G_tD[0])
        if G_tD[i]<4.5 and G_tD[i]>0.0 and slope[i]>mslop:
           mslop = slope[i]
           tangent_ind = i

    tangentx = np.linspace(G_tD[0], G_tD[tangent_ind]*1.2, 10)
    tangenty = mslop*tangentx + GdPdG[0]

    vlinex = [G_tD[tangent_ind], G_tD[tangent_ind]]
    vliney = [0.0, GdPdG[tangent_ind]*1.2]


    # Visualization
    N1 = 2
    fsize = 30
    msize = 16
    lw=6
    malpha = 1.0
    fig, ax = plt.subplots(2,2,figsize=(32, 18))

    ax[0,0].plot(tl, fpre[:,well_ind]/1.0e6, lw=lw, color='deepskyblue', alpha=0.8)
    ax[0,0].plot(tl[t_ind], fpre[t_ind,well_ind]/1.0e6, linestyle='None', marker='o', alpha=0.5, markerfacecolor='r', fillstyle='full', markersize=msize)
    #ax[0,0].set_xlim([1, 10])
    ax[0,0].set_ylim([30, 60])
    ax[0,0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0,0].set_ylabel(r'Pressure (MPa)', size=fsize, weight="bold", color='deepskyblue')
    ax[0,0].grid(True, which="both", ls="-")
    ax[0,0].xaxis.set_tick_params(labelsize=fsize)
    ax[0,0].yaxis.set_tick_params(labelsize=fsize, colors='deepskyblue')
    ax2 = ax[0,0].twinx()
    ax2.plot(Rate_time, Rate, lw=lw, linestyle='--', color='r', alpha=0.8 )
    ax2.set_ylim([0, 4.0])
    ax2.set_ylabel(r'Pumping Rate (kg/s)', size=fsize, weight="bold", color='r')
    ax2.yaxis.set_tick_params(labelsize=fsize, colors='r')


    ax[0,1].plot(tl, aper[:,well_ind]*1.0e3, lw=lw, color='cyan', alpha=0.8)
    ax[0,1].plot(tl[t_ind], aper[t_ind,well_ind]*1.0e3, linestyle='None', marker='o', alpha=0.5, markerfacecolor='r', fillstyle='full', markersize=msize)
    #ax[0,1].set_xlim([1, 10])
    #ax[0,1].set_ylim([-6.4, -5.8])
    ax[0,1].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0,1].set_ylabel(r'Frac Mouth Opeing (mm)', size=fsize, weight="bold")
    ax[0,1].grid(True, which="both", ls="-")
    ax[0,1].xaxis.set_tick_params(labelsize=fsize)
    ax[0,1].yaxis.set_tick_params(labelsize=fsize)


    ax[1,0].plot(tl, length, lw=lw, color='orangered', alpha=0.8)
    ax[1,0].plot(tl[t_ind], length[t_ind], linestyle='None', marker='o', alpha=0.5, markerfacecolor='r', fillstyle='full', markersize=msize)
    #ax[0,1].set_xlim([1, 10])
    #ax[1,0].set_ylim([-6.4, -5.8])
    ax[1,0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[1,0].set_ylabel(r'Frac Length (m)', size=fsize, weight="bold")
    ax[1,0].grid(True, which="both", ls="-")
    ax[1,0].xaxis.set_tick_params(labelsize=fsize)
    ax[1,0].yaxis.set_tick_params(labelsize=fsize)


    ax[1,1].plot(G_tD, P_tD/1.0e6, lw=lw, color='deepskyblue', alpha=0.8)
    ax[1,1].plot(G_tD[tangent_ind], P_tD[tangent_ind]/1.0e6, linestyle='None', marker='o', alpha=0.5, markerfacecolor='g', fillstyle='full', markersize=msize*1.5)
    ax[1,1].set_xlim([0, 4.5])
    ax[1,1].set_ylim([55, 58])
    ax[1,1].set_xlabel(r'G (Time)', size=fsize, weight="bold")
    ax[1,1].set_ylabel(r'Pressure (MPa)', size=fsize, weight="bold", color='deepskyblue')
    #ax[0,1].legend(loc='lower right',fontsize=fsize*0.8)
    ax[1,1].grid(True, which="both", ls="-")
    ax[1,1].xaxis.set_tick_params(labelsize=fsize)
    ax[1,1].yaxis.set_tick_params(labelsize=fsize, colors='deepskyblue')
    ax2 = ax[1,1].twinx()
    ax2.plot(G_tD[0:len(G_tD)-1], GdPdG, lw=lw, color='lime', alpha=0.8 )
    ax2.plot(tangentx, tangenty, lw=lw*0.5, linestyle='--', color='black', alpha=0.8 )
    ax2.plot(vlinex, vliney, lw=lw*0.5, linestyle='--', color='g', alpha=0.8 )
    ax2.set_ylim([0, 3.0e6])
    ax2.set_ylabel(r'GdP/dG', size=fsize, weight="bold", color='lime')
    ax2.yaxis.set_tick_params(labelsize=fsize, colors='lime')

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.3, hspace=0.3)

    fig.savefig('DFIT_KGD.png')

    plt.show()

if __name__ == "__main__":
    main()

