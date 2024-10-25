import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
from mpmath import *
import math
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

    mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus"])
    mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
    mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
    return mechanicalParameters


def getFracturePressureFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.findall('FieldSpecifications/FieldSpecification')
    pressure = 0.0
    found_traction = False
    for elem in param:
        if elem.get("fieldName") == "traction" and elem.get("component") == "0":
            pressure = float(elem.get("scale")) * (-1)
            found_traction = True
        if found_traction: break

    return pressure


def getFractureLengthFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    rectangle = tree.find('Geometry/Rectangle')
    dimensions = rectangle.get("dimensions")
    dimensions = [float(i) for i in dimensions[1:-1].split(",")]
    length = dimensions[0] / 2
    origin = rectangle.get("origin")
    origin = [float(i) for i in origin[1:-1].split(",")]

    return length, origin[0]


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()

    #-------- EmbeddeFrac File path
    outputDir = args.outputDir
    geosDir = args.geosDir
    hdf5File1Path = outputDir + "/displacementJump_embeddedFrac.hdf5"

    # Read HDF5
    hf = h5py.File(hdf5File1Path, 'r')
    jump = hf.get('displacementJump')
    jump = np.asarray(jump)
    aperture_EmbeddeFrac = jump[0, :, 0]
    x = hf.get('displacementJump elementCenter')
    loc_EmbeddeFrac = x[0, :, 1]

    #-------- Lagrange Contact File path
    hdf5File2Path = "displacementJump_contactMechanics.hdf5"

    # Read HDF5
    hf = h5py.File(hdf5File2Path, 'r')
    jump = hf.get('displacementJump')
    jump = np.asarray(jump)
    aperture_Contact = jump[0, :, 0]
    x = hf.get('displacementJump elementCenter')
    loc_Contact = x[0, :, 1]

    #-------- HydroFrac File path
    hdf5File3Path = "displacementJump_hydroFrac.hdf5"

    # Read HDF5
    hf = h5py.File(hdf5File3Path, 'r')
    jump = hf.get('elementAperture')
    jump = np.asarray(jump)
    aperture_HydroFrac = jump[0, :]
    x = hf.get('elementAperture elementCenter')
    loc_HydroFrac = x[0, :, 1]

    #-------- Extract info from XML
    xmlFilePath = geosDir + "/inputFiles/efemFractureMechanics/Sneddon_embeddedFrac"

    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath+"_base.xml")
    appliedPressure = getFracturePressureFromXML(xmlFilePath+"_base.xml")

    # Get length of the fracture
    length, origin = getFractureLengthFromXML(xmlFilePath+"_verification.xml")

    # Initialize Sneddon's analytical solution
    sneddonAnalyticalSolution = Sneddon(mechanicalParameters, length, appliedPressure)

    # Plot analytical (continuous line) and numerical (markers) aperture solution
    x_analytical = np.linspace(-length, length, 501, endpoint=True)
    aperture_analytical = np.empty(len(x_analytical))
    i = 0
    for xCell in x_analytical:
        aperture_analytical[i] = sneddonAnalyticalSolution.computeAperture(xCell)
        i += 1

    fsize = 30
    msize = 15
    lw = 6
    fig, ax = plt.subplots(1, figsize=(16, 12))
    cmap = plt.get_cmap("tab10")

    N1 = 1
    ax.plot(x_analytical, aperture_analytical * 1.0e3, color='k', label='Analytical Solution', lw=lw)
    ax.plot(loc_EmbeddeFrac[0::N1],
            aperture_EmbeddeFrac[0::N1] * 1.0e3,
            'o',
            color=cmap(0),
            label='Embedded Fracture',
            markersize=msize * 0.6,
            alpha=0.8)
    ax.plot(loc_Contact[0::N1],
            aperture_Contact[0::N1] * 1.0e3,
            's',
            color=cmap(4),
            label='Lagrange Contact',
            markersize=msize * 0.6,
            alpha=0.6)
    ax.plot(loc_HydroFrac[0::N1],
            aperture_HydroFrac[0::N1] * 1.0e3,
            'X',
            color=cmap(1),
            label='HydroFrac Solver',
            markersize=msize * 0.6,
            alpha=0.4)

    ax.grid()
    ax.set_xlabel('Fracture Length [m]', size=fsize, weight="bold")
    ax.set_ylabel('Fracture Aperture [mm]', size=fsize, weight="bold")
    ax.legend(bbox_to_anchor=(0.5, 0.2), loc='center', borderaxespad=0., fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)

    plt.show()


if __name__ == "__main__":
    main()
