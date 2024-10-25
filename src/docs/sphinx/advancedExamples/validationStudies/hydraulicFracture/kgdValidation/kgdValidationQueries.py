import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin, cos, tan, exp
import hdf5_wrapper
import xml.etree.ElementTree as ElementTree


# Get mesh settings from xml file
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


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Load and process GEOSX results
    # File path
    prefix = geosDir + "/inputFiles/hydraulicFracturing/"
    hdf5File = prefix + "KGD_validation_output.hdf5"
    xmlFile1Path = prefix + "kgdValidation_base.xml"
    xmlFile2Path = prefix + "kgdValidation_benchmark.xml"

    # Read simulation parameters from XML file
    xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz = getMeshSettings(xmlFile2Path)
    fracHeight = abs(zMax - zMin)
    dx = abs(xMax - xMin) / nx
    dy = abs(yMax - yMin) / ny
    dz = abs(zMax - zMin) / nz
    zmean = (zMax + zMin) / 2.

    # Read simulation output from HDF5 file
    # Global Coordinate of Element Center
    hf = hdf5_wrapper.hdf5_wrapper(hdf5File)
    xl = hf['pressure elementCenter']
    xl = np.asarray(xl)
    xcord = xl[-1, :, 0]
    ycord = xl[-1, :, 1]
    zcord = xl[-1, :, 2]
    tl = hf['pressure Time']
    tl = np.asarray(tl)
    # Load pressure
    fpre = hf['pressure']
    fpre = np.asarray(fpre)
    # Load elementAperture
    aper = hf['elementAperture']
    aper = np.asarray(aper)
    # Load elementArea
    area = hf['elementArea']
    area = np.asarray(area)

    # Query simulation results
    xloc_58 = 0.015
    xloc_57 = 0.041
    xloc_lvdt = 0.0285
    wellPressure = np.zeros([len(tl)])
    G58Pressure = np.zeros([len(tl)])
    G57Pressure = np.zeros([len(tl)])
    LVDTAperture = np.zeros([len(tl)])
    fracArea = np.zeros([len(tl)])
    for j in range(0, len(tl)):
        xcord = xl[j, :, 0]
        ycord = xl[j, :, 1]
        zcord = xl[j, :, 2]
        temp = 0
        for i in range(0, len(aper[0, :])):
            if abs(xcord[i] / (dx / 2.) - 1) < 0.01 and abs(ycord[i]) < 0.01 and abs(zcord[i] /
                                                                                     (zmean - dz / 2.) - 1) < 0.01:
                wellPressure[j] = fpre[j, i]

            if abs(xcord[i] - xloc_58) < 0.001 and abs(ycord[i]) < 0.01 and abs(zcord[i] /
                                                                                (zmean - dz / 2.) - 1) < 0.01:
                G58Pressure[j] = fpre[j, i]

            if abs(xcord[i] - xloc_57) < 0.001 and abs(ycord[i]) < 0.01 and abs(zcord[i] /
                                                                                (zmean - dz / 2.) - 1) < 0.01:
                G57Pressure[j] = fpre[j, i]

            if abs(xcord[i] - xloc_lvdt) < 0.001 and abs(ycord[i]) < 0.01 and abs(zcord[i] /
                                                                                  (zmean - dz / 2.) - 1) < 0.01:
                LVDTAperture[j] = aper[j, i]

            if aper[j, i] > 1.0e-5:
                temp += area[j, i]

        fracArea[j] = temp

    # Ouput fracture characteristics
    header = '      time  wpressure 58pressure 57pressure  Laperture       area'
    timehist = []
    for i in range(0, len(tl)):
        time = tl[i]
        pressure1 = wellPressure[i]
        pressure2 = G58Pressure[i]
        pressure3 = G57Pressure[i]
        aperture = LVDTAperture[i]
        farea = fracArea[i]
        timehist.append(
            [float(time),
             float(pressure1),
             float(pressure2),
             float(pressure3),
             float(aperture),
             float(farea)])

    np.savetxt('model_results.txt', timehist, header=header, fmt='%10.4g', comments='')


if __name__ == "__main__":
    main()
