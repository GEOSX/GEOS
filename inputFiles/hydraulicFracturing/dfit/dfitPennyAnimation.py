import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin,cos,tan,exp
import h5py
import xml.etree.ElementTree as ElementTree
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter


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
    xMax = nodeXCoords[1]
    yMin = nodeYCoords[0]
    yMax = nodeYCoords[1]
    zMin = nodeZCoords[0]
    zMax = nodeZCoords[1]

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
    t_shutDown = 120.0
    epsilon = 0.0001
    Rate = [0.0, q, q, 0.0, 0.0 ]
    Rate_time = [0.0, 0.0+ epsilon, t_shutDown - epsilon,  t_shutDown + epsilon, 400 ]

    # Load and process GEOSX results
    # File path
    hdf5File1Path = "dfit_penny_output.hdf5"
    xmlFile1Path = "dfit_pennyShaped_base.xml"
    xmlFile2Path = "dfit_pennyShaped_benchmark.xml"

    # Read simulation parameters from XML file
    xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz = getMeshSettings(xmlFile2Path)    
    dx = abs(xMax-xMin)/nx
    dy = abs(yMax-yMin)/ny

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
        if abs(xcord[i]/(dx/2.) - 1) < 0.01 and abs(ycord[i]/(dy/2.) - 1) < 0.01 and abs(zcord[i]) < 0.01:
           well_ind = i

    # Find shut-in time
    for i in range(0,len(tl)):
        if abs(tl[i]/t_shutDown - 1.0) < 0.01:
           t_ind = i

    # Compute frac length
    radius = np.zeros([len(tl)])     
    for i in range(0,len(tl)):
        temp = 0
        for j in range(0,len(aper[0,:])):
            if aper[i,j]>1.0e-5:
               temp += area[i,j]

        radius[i] = (temp*4.0/math.pi)**0.5

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
        if G_tD[i]<5.0 and G_tD[i]>0.0 and slope[i]>mslop:
           mslop = slope[i]
           tangent_ind = i

    tangentx = np.linspace(G_tD[0], G_tD[tangent_ind]*1.2, 10)
    tangenty = mslop*tangentx + GdPdG[0]

    vlinex = [G_tD[tangent_ind], G_tD[tangent_ind]]
    vliney = [0.0, GdPdG[tangent_ind]*1.2]

    diag_ind = 21
    vx1 = [tl[diag_ind], tl[diag_ind]]
    vy1 = [0, 80]

    # Visulization
    N1 = 2
    fsize = 30
    msize = 150
    lw=6
    malpha = 1.0
    fig, ax = plt.subplots(2,2,figsize=(32, 18), dpi=36)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.3, hspace=0.3)

    def animate_func(num):
        ax[0,0].clear()
        ax[0,0].plot(tl[ :num+1], fpre[ :num+1, well_ind]/1.0e6, lw=lw, color='deepskyblue', alpha=0.8)
        ax[0,0].scatter(tl[num], fpre[num, well_ind]/1.0e6, c='k', marker='o', s=msize)
        ax[0,0].set_xlim([0, tl[num+1]])
        ax[0,0].set_ylim([30, 60]) 
        ax[0,0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
        ax[0,0].set_ylabel(r'Pressure (MPa)', size=fsize, weight="bold", color='deepskyblue')
        ax[0,0].grid(True, which="both", ls="-")
        ax[0,0].xaxis.set_tick_params(labelsize=fsize)
        ax[0,0].yaxis.set_tick_params(labelsize=fsize, colors='deepskyblue')
        ax2 = ax[0,0].twinx()
        ax2.plot(Rate_time, Rate, lw=lw, linestyle='--', color='r', alpha=0.8 )
        ax2.set_xlim([0, tl[num+1]])
        ax2.set_ylim([0, 4.0])
        ax2.set_ylabel(r'Pumping Rate (kg/s)', size=fsize, weight="bold", color='r')
        ax2.yaxis.set_tick_params(labelsize=fsize, colors='r') 
     
        ax[0,1].clear()
        ax[0,1].plot(tl[ :num+1], aper[ :num+1, well_ind]*1.0e3, lw=lw, color='cyan', alpha=0.8)
        ax[0,1].scatter(tl[num], aper[num, well_ind]*1.0e3, c='k', marker='o', s=msize)
        ax[0,1].set_xlim([0, tl[num+1]])
        ax[0,1].set_ylim([0, 2.5])
        ax[0,1].set_xlabel(r'Time (s)', size=fsize, weight="bold")
        ax[0,1].set_ylabel(r'Frac Mouth Opeing (mm)', size=fsize, weight="bold")
        ax[0,1].grid(True, which="both", ls="-")
        ax[0,1].xaxis.set_tick_params(labelsize=fsize)
        ax[0,1].yaxis.set_tick_params(labelsize=fsize)  
       
        ax[1,0].clear()
        ax[1,0].plot(tl[ :num+1], radius[ :num+1], lw=lw, color='orangered', alpha=0.8)
        ax[1,0].scatter(tl[num], radius[num], c='k', marker='o', s=msize)        
        ax[1,0].set_xlim([0, tl[num+1]])
        ax[1,0].set_ylim([0, 21])
        ax[1,0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
        ax[1,0].set_ylabel(r'Frac Radius (m)', size=fsize, weight="bold")
        ax[1,0].grid(True, which="both", ls="-")
        ax[1,0].xaxis.set_tick_params(labelsize=fsize)
        ax[1,0].yaxis.set_tick_params(labelsize=fsize) 


    ax[1,1].plot(G_tD, P_tD/1.0e6, lw=lw, color='deepskyblue', alpha=0.8)
    ax[1,1].plot(G_tD[tangent_ind], P_tD[tangent_ind]/1.0e6, linestyle='None', marker='o', alpha=0.5, markerfacecolor='g', fillstyle='full', markersize=15)
    ax[1,1].set_xlim([0, 5])
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
 

    ani = animation.FuncAnimation(fig, animate_func, interval=20, frames=t_ind+20, repeat=False)

    # Save the animation as an animated GIF
    ani.save("dfit_Penny.gif", dpi=20, writer=PillowWriter(fps=1))
    
    plt.show()

if __name__ == "__main__":
    main()

