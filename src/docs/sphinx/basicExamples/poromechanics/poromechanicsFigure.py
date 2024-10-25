import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
from mpmath import *
import math
import os
import argparse


class terzaghi:

    def __init__(self, hydromechanicalParameters, xMin, xMax, appliedTraction):
        E = hydromechanicalParameters["youngModulus"]
        nu = hydromechanicalParameters["poissonRation"]
        b = hydromechanicalParameters["biotCoefficient"]
        mu = hydromechanicalParameters["fluidViscosity"]
        cf = hydromechanicalParameters["fluidCompressibility"]
        phi = hydromechanicalParameters["porosity"]
        k = hydromechanicalParameters["permeability"]

        K = E / 3.0 / (1.0 - 2.0 * nu)    # bulk modulus
        Kv = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))    # uniaxial bulk modulus
        Se = (b - phi) * (1.0 - b) / K + phi * cf    # constrained specific storage

        self.characteristicLength = xMax - xMin
        self.appliedTraction = abs(appliedTraction)
        self.loadingEfficiency = b / (Kv * Se + b**2)
        self.consolidationCoefficient = (k / mu) * Kv / (Se * Kv + b**2)
        self.consolidationTime = self.characteristicLength**2 / self.consolidationCoefficient
        self.initialPressure = self.loadingEfficiency * self.appliedTraction

    def computePressure(self, x, t):
        if t == 0.0:
            return self.initialPressure
        else:
            cc = self.consolidationCoefficient
            L = self.characteristicLength
            p = nsum(
                lambda m: 1 / (2 * m + 1) * exp(-((2 * m + 1)**2) * (math.pi**2) * cc * t / 4 / L / L) * sin(
                    (2 * m + 1) * math.pi * x / 2 / L), [0, inf])
            return 4 * self.initialPressure / math.pi * p


def getHydromechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param1 = tree.find('Constitutive/ElasticIsotropic')
    param2 = tree.find('Constitutive/BiotPorosity')
    param3 = tree.find('Constitutive/CompressibleSinglePhaseFluid')
    param4 = tree.find('Constitutive/ConstantPermeability')

    hydromechanicalParameters = dict.fromkeys([
        "youngModulus", "poissonRation", "biotCoefficient", "fluidViscosity", "fluidCompressibility", "porosity",
        "permeability"
    ])

    hydromechanicalParameters["youngModulus"] = float(param1.get("defaultYoungModulus"))
    hydromechanicalParameters["poissonRation"] = float(param1.get("defaultPoissonRatio"))

    E = hydromechanicalParameters["youngModulus"]
    nu = hydromechanicalParameters["poissonRation"]
    K = E / 3.0 / (1.0 - 2.0 * nu)
    Kg = float(param2.get("defaultGrainBulkModulus"))

    hydromechanicalParameters["biotCoefficient"] = 1.0 - K / Kg
    hydromechanicalParameters["porosity"] = float(param2.get("defaultReferencePorosity"))
    hydromechanicalParameters["fluidViscosity"] = float(param3.get("defaultViscosity"))
    hydromechanicalParameters["fluidCompressibility"] = float(param3.get("compressibility"))

    perm = param4.get("permeabilityComponents")
    perm = np.array(perm[1:-1].split(','), float)
    hydromechanicalParameters["permeability"] = perm[0]

    return hydromechanicalParameters


def getAppliedTractionFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)
    param = tree.find('FieldSpecifications/Traction')
    return float(param.get("scale"))


def getDomainMaxMinXCoordFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)
    meshElement = tree.find('Mesh/InternalMesh')
    nodeXCoords = meshElement.get("xCoords")
    nodeXCoords = [float(i) for i in nodeXCoords[1:-1].split(",")]
    xMin = nodeXCoords[0]
    xMax = nodeXCoords[-1]
    return xMin, xMax


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()
    outputDir = args.outputDir
    geosDir = args.geosDir

    # File path
    hdf5FilePath = outputDir + "/pressure_history.hdf5"
    xmlBaseFilePath = geosDir + "/inputFiles/poromechanics/PoroElastic_Terzaghi_base_direct.xml"
    xmlSmokeFilePath = geosDir + "/inputFiles/poromechanics/PoroElastic_Terzaghi_smoke.xml"

    # Read HDF5
    hf = h5py.File(hdf5FilePath, 'r')
    time = hf.get('pressure Time')
    pressure = hf.get('pressure')
    x = hf.get('pressure elementCenter')

    # Extract info from XML
    hydromechanicalParameters = getHydromechanicalParametersFromXML(xmlBaseFilePath)
    appliedTraction = getAppliedTractionFromXML(xmlBaseFilePath)

    # Get domain min/max coordinate in the x-direction
    xMin, xMax = getDomainMaxMinXCoordFromXML(xmlSmokeFilePath)

    # Initialize Terzaghi's analytical solution
    terzaghiAnalyticalSolution = terzaghi(hydromechanicalParameters, xMin, xMax, appliedTraction)

    # Plot analytical (continuous line) and numerical (markers) pressure solution
    x_analytical = np.linspace(xMin, xMax, 51, endpoint=True)
    pressure_analytical = np.empty(len(x_analytical))

    cmap = plt.get_cmap("tab10")
    iplt = -1
    for k in range(0, len(time), 2):
        iplt += 1
        t = time[k, 0]
        i = 0
        for xCell in x_analytical:
            xScaled = terzaghiAnalyticalSolution.characteristicLength * (xCell - xMin) / (xMax - xMin)
            pressure_analytical[i] = terzaghiAnalyticalSolution.computePressure(xScaled, t)
            i += 1
        plt.plot(x_analytical, pressure_analytical, color=cmap(iplt), label='t = ' + str(t) + ' s')
        plt.plot(x[k, :, 0], pressure[k, :], 'o', color=cmap(iplt))

    plt.grid()
    plt.xlabel('$x$ [m]')
    plt.ylabel('pressure [Pa]')
    plt.legend(bbox_to_anchor=(0.1, 0.55), loc='lower left', borderaxespad=0.)
    plt.show()


if __name__ == "__main__":
    main()
