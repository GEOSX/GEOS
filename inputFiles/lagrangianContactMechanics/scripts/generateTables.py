import numpy as np
import os
import sys
import xml.etree.ElementTree as ElementTree
import matplotlib
import matplotlib.pyplot as plt
    
class SingularCrackSlip:

    def __init__(self, mechanicalParameters, length ):
        K = mechanicalParameters["bulkModulus"]
        G = mechanicalParameters["shearModulus"]
        poisson_ratio= (3 * K - 2 * G) / (2 * (3 * K + G))
        
        mu_star = G /( 1 - poisson_ratio)
        self.tau_0 = 0.0
        self.tau_r = -1.0

        self.scaling = 2*(self.tau_0 - self.tau_r)/mu_star
        self.halfLength = length

    def computeSlip(self, x):
        return self.scaling * np.sqrt(self.halfLength**2 - x**2)

    def computeTraction(self, x):
        if x < -self.halfLength or x > self.halfLength:
            return self.tau_0 + (self.tau_0-self.tau_r) * ( np.abs(x)/np.sqrt(x**2 - self.halfLength**2) - 1 )
        else:
            return self.tau_r
class GaussianSlip:

    def __init__(self, peakStrength, length ):
        self.scaling =  peakStrength
        self.halfLength = length

    def computeSlip(self, x):
        denom = 1 / (self.halfLength/2)
        return self.scaling*np.exp(-0.5*((x)/denom)**2)

def getMechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.find('Constitutive/ElasticIsotropic')

    mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus"])
    mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
    mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
    return mechanicalParameters

def getFractureLengthFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    rectangle = tree.find('Geometry/Box')
    xmin = rectangle.get("xMin")
    xmax = rectangle.get("xMax")
    xmin = [float(i) for i in xmin[1:-1].split(",")]
    xmax = [float(i) for i in xmax[1:-1].split(",")]
    length = ( xmax[0] - xmin[0] ) / 2
    origin = 0.0
   
    return length, origin

def plot_traction_solution():
    # Read HDF5
    import hdf5_wrapper
    hdf5File1Path = "Output/traction.hdf5"

    # Read HDF5
    data = hdf5_wrapper.hdf5_wrapper(hdf5File1Path).get_copy()
    traction = data['traction']
    traction = np.asarray(traction)
    traction_geos = traction[0, :, 1]
    x = data['traction elementCenter']
    x_geos = x[0, :, 0]

     #-------- Extract info from XML
    xmlFilePath = "/Users/cusini1/geosx/GEOSX/inputFiles/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_base.xml"

    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)

    # Get length of the fracture
    xmlFilePath = "/Users/cusini1/geosx/GEOSX/inputFiles/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_smoke.xml"
    totalHalfLength, originShift = getFractureLengthFromXML(xmlFilePath)
    halfLength = 2.0

    singularCrackSlipSolution = SingularCrackSlip(mechanicalParameters, halfLength)
    x = np.linspace(-totalHalfLength, totalHalfLength, 10000, endpoint=True)
    traction_analytical = np.zeros(len(x))
    i = 0
    for xCell in x:
        traction_analytical[i] = singularCrackSlipSolution.computeTraction(xCell)
        i += 1

    fsize = 30
    msize = 15
    lw = 2
    fig, ax = plt.subplots(1, figsize=(16, 12))
    cmap = plt.get_cmap("tab10")
    
    # Plot analytical (continuous line) and numerical (markers) aperture solution
    ax.plot(x, traction_analytical, color='r', label='Traction analytical', lw=lw)
    ax.plot(x_geos, traction_geos, color='k', label='geos', marker="o", lw=lw)

    ax.set_xlabel('Fault coordinate [m]', size=fsize, weight="bold")
    ax.set_ylabel('Shear traction', size=fsize, weight="bold")
    ax.legend(bbox_to_anchor=(0.75, 0.9), loc='center', borderaxespad=0., fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)
    plt.savefig("traction.png")

def output_tables(x, slip, name):
    # Save x to x.csv with one value per row
    np.savetxt('x.csv', x, fmt='%f')

    # Save aperture_analytical to jump.csv with one value per row
    np.savetxt(f'{name}.csv', slip, fmt='%f')      


def debug():
    #-------- Extract info from XML
    xmlFilePath = "/Users/cusini1/geosx/GEOSX/inputFiles/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_base.xml"

    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
    appliedPressure = 1.0

    # Get length of the fracture
    xmlFilePath = "/Users/cusini1/geosx/GEOSX/inputFiles/lagrangianContactMechanics/LagrangeContactBubbleStab_FixedSlip_smoke.xml"
    totalHalfLength, originShift = getFractureLengthFromXML(xmlFilePath)
    halfLength = 2.0

    # Initialize Sneddon's analytical solution
    singularCrackSlipSolution = SingularCrackSlip(mechanicalParameters, halfLength )
    peakStrength = 3.0
    gaussianSlipSolution = GaussianSlip( peakStrength, halfLength)

    # Plot analytical (continuous line) and numerical (markers) aperture solution
    x = np.linspace(-totalHalfLength, totalHalfLength, 100000, endpoint=True)
    singularCrackSlip = np.zeros(len(x))
    gaussianSlip = np.zeros(len(x))
    i = 0
    for xCell in x:
        if xCell > -halfLength and xCell < halfLength:
            singularCrackSlip[i] = singularCrackSlipSolution.computeSlip(xCell)
        gaussianSlip[i] = gaussianSlipSolution.computeSlip(xCell)
        i += 1

    fsize = 24
    msize = 15
    lw = 6
    fig, ax = plt.subplots(1, figsize=(16, 12))
    cmap = plt.get_cmap("tab10")

    ax.plot(x, singularCrackSlip , color='k', label='Analytical Solution', lw=lw)
    ax.grid()
    ax.set_xlabel('Fault coordinate [m]', size=fsize, weight="bold")
    ax.set_ylabel('slip [m]', size=fsize, weight="bold")
    ax.legend(bbox_to_anchor=(0.7, 1), loc='center', borderaxespad=0., fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)
    plt.savefig("singularCrackSlip.png")

    fig, ax = plt.subplots(1, figsize=(16, 12))
    cmap = plt.get_cmap("tab10")

    ax.plot(x, gaussianSlip , color='k', label='Analytical Solution', lw=lw)
    ax.grid()
    ax.set_xlabel('Fault coordinate [m]', size=fsize, weight="bold")
    ax.set_ylabel('slip [m]', size=fsize, weight="bold")
    ax.legend(bbox_to_anchor=(0.75, 0.9), loc='center', borderaxespad=0., fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)
    plt.savefig("gaussianSlip.png")

    output_tables(x, singularCrackSlip, "singularCrackSlip")
    output_tables(x, gaussianSlip, "gaussianSlip")

if __name__ == "__main__":
    debug()
    plot_traction_solution()
