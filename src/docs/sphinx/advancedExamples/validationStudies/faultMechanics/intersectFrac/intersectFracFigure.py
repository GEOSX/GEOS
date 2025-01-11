import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
from mpmath import *
import math
from math import sin, cos, tan, exp, atan, asin
from mpl_toolkits.mplot3d import axes3d
import os
import argparse


class Sneddon:

    def __init__(self, mechanicalParameters, length, pressure):
        K = mechanicalParameters["bulkModulus"]
        G = mechanicalParameters["shearModulus"]
        E = (9 * K * G) / (3 * K + G)
        nu = E / (2 * G) - 1

        self.scaling = (4 * (1 - nu**2)) * pressure / E
        self.halfLength = length

    def computeAperture(self, x):
        return self.scaling * (self.halfLength**2 - x**2)**0.5


def getMechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.find('Constitutive/ElasticIsotropic')

    mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus", "frictionCoefficient"])
    mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
    mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))

    param = tree.find('Constitutive/Coulomb')
    mechanicalParameters["frictionCoefficient"] = float(param.get("frictionCoefficient"))
    return mechanicalParameters


def getCompressiveStressFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.findall('FieldSpecifications/FieldSpecification')

    found_stress = False
    for elem in param:
        if elem.get("fieldName") == "rock_stress" and elem.get("component") == "1":
            stress = float(elem.get("scale")) * (-1)
            found_stress = True
        if found_stress: break

    return stress


def getFracturePressureFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.findall('FieldSpecifications/Traction')

    found_pressure = False
    for elem in param:
        if elem.get("name") == "NormalTraction" and elem.get("tractionType") == "normal":
            pressure = float(elem.get("scale")) * (-1)
            found_pressure = True
        if found_pressure: break

    return pressure


def getFractureGeometryFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.findall('Geometry/Rectangle')

    for elem in param:
        if elem.get("name") == "fracture1":
            dimensions = elem.get("dimensions")
            dimensions = [float(i) for i in dimensions[1:-1].split(",")]
            length1 = dimensions[0] / 2
        elif elem.get("name") == "fracture2":
            dimensions = elem.get("dimensions")
            dimensions = [float(i) for i in dimensions[1:-1].split(",")]
            length2 = dimensions[0] / 2

    return length1, length2


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
    hdf5File1Path = outputDir + "/traction_history.hdf5"
    hdf5File2Path = outputDir + "/displacementJump_history.hdf5"
    xmlFilePath = geosDir + "/inputFiles/lagrangianContactMechanics/TFrac_base.xml"

    # Read HDF5
    # Global Coordinate of Fracture Element Center
    hf = h5py.File(hdf5File1Path, 'r')
    xl = hf.get('traction elementCenter')
    xl = np.asarray(xl)
    xcord = xl[0, :, 0]
    ycord = xl[0, :, 1]
    zcord = xl[0, :, 2]

    # Local Normal Traction
    trac = hf.get('traction')
    trac = np.asarray(trac)
    normalTraction = trac[-1, :, 0]

    # Local Shear Displacement
    hf = h5py.File(hdf5File2Path, 'r')
    jump = hf.get('displacementJump')
    jump = np.asarray(jump)
    displacementJump = jump[-1, :, 1]
    aperture = jump[-1, :, 0]

    # Extract Local Inform for The Horizontal Fracture
    xlist = []
    tnlist = []
    gtlist = []
    for i in range(0, len(zcord)):
        if abs(ycord[i] / 50.0 - 1.) < 0.01:
            xlist.append(xcord[i])
            tnlist.append(-normalTraction[i] * 1.0e-6)
            gtlist.append(displacementJump[i] * 1.e3)

    # Extract Local Inform for The Vertical Fracture
    ylist = []
    apertlist = []
    for i in range(0, len(zcord)):
        if abs(xcord[i] - 0.) < 0.01:
            ylist.append(ycord[i])
            apertlist.append(aperture[i] * 1.e3)

    # Load numerical solutions from literature work
    r1, tn_literature = np.loadtxt('NormalTraction.txt', skiprows=0, unpack=True)
    r2, gt_literature = np.loadtxt('Slip.txt', skiprows=0, unpack=True)
    r3, ap_literature = np.loadtxt('Aperture.txt', skiprows=0, unpack=True)

    # Extract Mechanical Properties and Fracture Geometry from XML
    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
    compressiveStress = getCompressiveStressFromXML(xmlFilePath)
    length1, length2 = getFractureGeometryFromXML(xmlFilePath)
    appliedPressure = getFracturePressureFromXML(xmlFilePath)

    # Initialize Sneddon's analytical solution
    sneddonAnalyticalSolution = Sneddon(mechanicalParameters, length1, appliedPressure)

    # Plot analytical (continuous line) and numerical (markers) aperture solution
    x_analytical = np.linspace(-length1, length1, 101, endpoint=True)
    aperture_analytical = np.empty(len(x_analytical))
    i = 0
    for xCell in x_analytical:
        aperture_analytical[i] = sneddonAnalyticalSolution.computeAperture(xCell)
        i += 1

    fsize = 30
    msize = 12
    lw = 6
    fig, ax = plt.subplots(2, 2, figsize=(32, 18))
    cmap = plt.get_cmap("tab10")

    ax[0, 0].plot(r1 - length2, tn_literature, color=cmap(-1), label='Phan et al.(2003)', lw=lw)
    ax[0, 0].plot(xlist, tnlist, 'o', alpha=0.6, color=cmap(2), mec='k', label='GEOSX Results', markersize=msize)
    ax[0, 0].grid()
    ax[0, 0].set_xlabel('Horizontal Frac Length [m]', size=fsize, weight="bold")
    ax[0, 0].set_ylabel('Normal Traction [MPa]', size=fsize, weight="bold")
    ax[0, 0].legend(loc='lower left', fontsize=fsize * 0.8)
    ax[0, 0].xaxis.set_tick_params(labelsize=fsize)
    ax[0, 0].yaxis.set_tick_params(labelsize=fsize)

    ax[0, 1].plot(r2 - length2, gt_literature, color=cmap(-1), label='Phan et al.(2003)', lw=lw)
    ax[0, 1].plot(xlist, gtlist, 'o', alpha=0.6, color=cmap(2), mec='k', label='GEOSX Results', markersize=msize)
    ax[0, 1].grid()
    ax[0, 1].set_xlabel('Horizontal Frac Length [m]', size=fsize, weight="bold")
    ax[0, 1].set_ylabel('Slip [mm]', size=fsize, weight="bold")
    ax[0, 1].legend(loc='lower left', fontsize=fsize * 0.8)
    ax[0, 1].xaxis.set_tick_params(labelsize=fsize)
    ax[0, 1].yaxis.set_tick_params(labelsize=fsize)

    N1 = 2
    ax[1, 0].plot(-r3 + length1, ap_literature, color=cmap(-1), label='Phan et al.(2003)', lw=lw)
    ax[1, 0].plot(ylist[0::N1],
                  apertlist[0::N1],
                  'o',
                  alpha=0.6,
                  color=cmap(2),
                  mec='k',
                  label='GEOSX Results',
                  markersize=msize)
    ax[1, 0].plot(x_analytical, aperture_analytical * 1.0e3, '--', color=cmap(1), label='Sneddon Solution', lw=lw)
    ax[1, 0].grid()
    ax[1, 0].set_xlabel('Vertical Frac Length [m]', size=fsize, weight="bold")
    ax[1, 0].set_ylabel('Aperture [mm]', size=fsize, weight="bold")
    ax[1, 0].legend(loc='lower right', fontsize=fsize * 0.8)
    ax[1, 0].xaxis.set_tick_params(labelsize=fsize)
    ax[1, 0].yaxis.set_tick_params(labelsize=fsize)

    ax[1, 1].axis('off')

    plt.show()


if __name__ == "__main__":
    main()
