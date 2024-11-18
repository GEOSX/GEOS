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


class Analytical:

    def __init__(self, mechanicalParameters, length, inclination, stress):
        K = mechanicalParameters["bulkModulus"]
        G = mechanicalParameters["shearModulus"]
        E = (9 * K * G) / (3 * K + G)
        nu = E / (2 * G) - 1

        self.halfLength = length
        self.inc = inclination
        self.stress = stress
        self.scaling = (4 * (1 - nu**2)) / E
        self.frictionCoefficient = mechanicalParameters["frictionCoefficient"]

    def computeNormalTraction(self, x):
        return -self.stress * pow(sin(self.inc), 2)

    def computeShearDisplacement(self, x):
        return self.scaling * (self.stress * sin(self.inc) *
                               (cos(self.inc) - sin(self.inc) * self.frictionCoefficient)) * pow(
                                   (self.halfLength**2 - (self.halfLength - x - 1.)**2), 0.5)


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
        if elem.get("fieldName") == "rock_stress" and elem.get("component") == "0":
            stress = float(elem.get("scale")) * (-1)
            found_stress = True
        if found_stress: break

    return stress


def getFractureGeometryFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    rectangle = tree.find('Geometry/Rectangle')
    dimensions = rectangle.get("dimensions")
    dimensions = [float(i) for i in dimensions[1:-1].split(",")]
    length = dimensions[0] / 2
    origin = rectangle.get("origin")
    origin = [float(i) for i in origin[1:-1].split(",")]
    direction = rectangle.get("lengthVector")
    direction = [float(i) for i in direction[1:-1].split(",")]
    inclination = atan(direction[1] / direction[0])

    return length, origin[0], inclination


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
    xmlFile1Path = geosDir + "/inputFiles/lagrangianContactMechanics/SingleFracCompression_base.xml"
    xmlFile2Path = geosDir + "/inputFiles/lagrangianContactMechanics/SingleFracCompression_benchmark.xml"

    # Read HDF5
    # Global Coordinate of Fracture Element Center
    hf = h5py.File(hdf5File1Path, 'r')
    xl = hf.get('traction elementCenter')
    xl = np.asarray(xl)
    xcord = xl[0, :, 0]
    ycord = xl[0, :, 1]
    zcord = xl[0, :, 2]

    # Local Normal Traction
    hf = h5py.File(hdf5File1Path, 'r')
    trac = hf.get('traction')
    trac = np.asarray(trac)
    normalTraction = trac[0, :, 0]

    # Local Shear Displacement
    hf = h5py.File(hdf5File2Path, 'r')
    jump = hf.get('displacementJump')
    jump = np.asarray(jump)
    displacementJump = jump[0, :, 1]

    # Extract Local Inform for The Middle Layer
    xlist = []
    ylist = []
    xloc = []
    tnlist = []
    gtlist = []
    for i in range(0, len(zcord)):
        if abs(zcord[i] / 0.025 - 1.) < 0.01:
            xlist.append(xcord[i])
            ylist.append(ycord[i])
            xloc.append(pow(xcord[i]**2 + ycord[i]**2, 0.5) * xcord[i] / abs(xcord[i]))
            tnlist.append(normalTraction[i] / 1.0e6)
            gtlist.append(displacementJump[i] * 1.e3)

    # Extract Mechanical Properties and Fracture Geometry from XML
    mechanicalParameters = getMechanicalParametersFromXML(xmlFile1Path)
    compressiveStress = getCompressiveStressFromXML(xmlFile1Path)
    length, origin, inclination = getFractureGeometryFromXML(xmlFile2Path)

    # Initialize analytical solution
    AnalyticalSolution = Analytical(mechanicalParameters, length, inclination, compressiveStress)

    # Plot Analytical (continuous line) and Numerical (markers) Solution
    x_analytical = np.linspace(-length, length, 101, endpoint=True)
    tn_analytical = np.empty(len(x_analytical))
    gt_analytical = np.empty(len(x_analytical))
    i = 0
    for xCell in x_analytical:
        tn_analytical[i] = AnalyticalSolution.computeNormalTraction(xCell) / 1.0e6
        gt_analytical[i] = AnalyticalSolution.computeShearDisplacement(xCell) * 1.e3
        i += 1

    fsize = 30
    msize = 10
    lw = 8
    fig, ax = plt.subplots(1, 2, figsize=(32, 12))
    cmap = plt.get_cmap("tab10")

    ax[0].plot(x_analytical, tn_analytical, color=cmap(-1), label='Analytical Solution', lw=lw)
    ax[0].plot(xloc, tnlist, 'o', alpha=0.6, color=cmap(2), mec='k', label='Numerical Solution', markersize=msize)
    ax[0].grid()
    ax[0].set_xlim(-1, 1)
    ax[0].set_ylim(-18, 2)
    ax[0].set_xlabel('Length [m]', size=fsize, weight="bold")
    ax[0].set_ylabel('Normal Traction [MPa]', size=fsize, weight="bold")
    ax[0].legend(bbox_to_anchor=(0.5, 0.2), loc='center', borderaxespad=0., fontsize=fsize)
    ax[0].xaxis.set_tick_params(labelsize=fsize)
    ax[0].yaxis.set_tick_params(labelsize=fsize)

    ax[1].plot(x_analytical, gt_analytical, color=cmap(-1), label='Analytical Solution', lw=lw)
    ax[1].plot(xloc, gtlist, 'o', alpha=0.6, color=cmap(2), mec='k', label='Numerical Solution', markersize=msize)
    ax[1].grid()
    ax[1].set_xlim(-1, 1)
    ax[1].set_ylim(0, 4)
    ax[1].set_xlabel('Length [m]', size=fsize, weight="bold")
    ax[1].set_ylabel('Relative Shear Displacement [mm]', size=fsize, weight="bold")
    ax[1].legend(bbox_to_anchor=(0.5, 0.2), loc='center', borderaxespad=0., fontsize=fsize)
    ax[1].xaxis.set_tick_params(labelsize=fsize)
    ax[1].yaxis.set_tick_params(labelsize=fsize)

    plt.show()


if __name__ == "__main__":
    main()
