import matplotlib
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import numpy as np
import h5py
import os
import argparse


def getMaxTime(xmlFilePath):

    tree = ElementTree.parse(xmlFilePath)
    eventElement = tree.find('Events')

    return float(eventElement.get('maxTime'))


def getDomainMaxMinCoords(xmlFilePath):

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

    return xMin, xMax, yMin, yMax, zMin, zMax


def getPorosity(xmlFilePath):

    tree = ElementTree.parse(xmlFilePath)
    refPorosity = 0
    porosityElement = tree.find('Constitutive/PressurePorosity')
    refPorosity = porosityElement.get('defaultReferencePorosity')
    return float(refPorosity)


def getInjectionRate(xmlFilePath):

    tree = ElementTree.parse(xmlFilePath)
    sourceFluxElement = tree.find('FieldSpecifications/SourceFlux')

    return abs(float(sourceFluxElement.get("scale")))


def getBrooksCoreyRelperm(xmlFilePath):

    tree = ElementTree.parse(xmlFilePath)
    brooksCoreyElement = tree.find('Constitutive/BrooksCoreyRelativePermeability')
    phaseMinVolFrac = brooksCoreyElement.get("phaseMinVolumeFraction")
    phaseRelPermExp = brooksCoreyElement.get("phaseRelPermExponent")
    phaseRelPermMaxValue = brooksCoreyElement.get("phaseRelPermMaxValue")

    phaseMinVolFrac = [float(i) for i in phaseMinVolFrac[1:-1].split(",")]
    phaseRelPermExp = [float(i) for i in phaseRelPermExp[1:-1].split(",")]
    phaseRelPermMaxValue = [float(i) for i in phaseRelPermMaxValue[1:-1].split(",")]

    return phaseMinVolFrac, phaseRelPermExp, phaseRelPermMaxValue

def getViscosity( pvdgFilePath, \
                  pvtwFilePath ):

    pvdgFile = open(pvdgFilePath, 'r')
    pvtwFile = open(pvtwFilePath, 'r')

    for i, line in enumerate(pvdgFile):
        if i == 1:
            visco = float(line.split(" ")[2])
    for i, line in enumerate(pvtwFile):
        if i == 1:
            viscw = float(line.split(" ")[3])

    return [visco, viscw]


def getDensity(xmlFilePath):

    tree = ElementTree.parse(xmlFilePath)
    brooksCoreyElement = tree.find('Constitutive/DeadOilFluid')
    surfaceDensities = brooksCoreyElement.get("surfaceDensities")

    surfaceDensities = [float(i) for i in surfaceDensities[1:-1].split(",")]

    return surfaceDensities

def computeFractionalFlow( phaseVolFrac, \
                           phaseMinVolFrac, \
                           phaseRelPermExp, \
                           phaseRelPermMaxValue, \
                           phaseVisc ):

    scale = 1 - phaseMinVolFrac[0] - phaseMinVolFrac[1]
    phaseMob = [0, 0]
    dPhaseMob_dS = [0, 0]

    for i in range(0, 2):
        scaledPhaseVolFrac = (phaseVolFrac[i] - phaseMinVolFrac[i]) / scale
        phaseMob[i] = phaseRelPermMaxValue[i] / phaseVisc[i] * scaledPhaseVolFrac**phaseRelPermExp[i]
        dPhaseMob_dS[i] = phaseRelPermExp[i] * phaseRelPermMaxValue[i] / phaseVisc[i] * (scaledPhaseVolFrac**(
            phaseRelPermExp[i] - 1)) / scale

    fractionalFlow = phaseMob[0] / (phaseMob[0] + phaseMob[1])
    dFractionalFlow_dS = (dPhaseMob_dS[0] * phaseMob[1] + phaseMob[0] * dPhaseMob_dS[1]) / (
        (phaseMob[0] + phaseMob[1])**2)

    return fractionalFlow, dFractionalFlow_dS


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()


    # File paths
    outputDir = args.outputDir
    geosDir = args.geosDir
    hdf5FilePath = outputDir + "/saturationHistory.hdf5"
    pvdgFilePath = geosDir + "/inputFiles/compositionalMultiphaseFlow/benchmarks/buckleyLeverettProblem/buckleyLeverett_table/pvdg.txt"
    pvtwFilePath = geosDir + "/inputFiles/compositionalMultiphaseFlow/benchmarks/buckleyLeverettProblem/buckleyLeverett_table/pvtw.txt"
    xmlFile1Path = geosDir + "/inputFiles/compositionalMultiphaseFlow/benchmarks/buckleyLeverettProblem/buckleyLeverett_base.xml"
    xmlFile2Path = geosDir + "/inputFiles/compositionalMultiphaseFlow/benchmarks/buckleyLeverettProblem/buckleyLeverett_benchmark.xml"

    # Read simulation parameters from XML file
    xMin, xMax, yMin, yMax, zMin, zMax = getDomainMaxMinCoords(xmlFile2Path)
    maxTime = getMaxTime(xmlFile2Path)
    porosity = getPorosity(xmlFile1Path)
    phaseMinVolFrac, phaseRelPermExp, phaseRelPermMaxValue = getBrooksCoreyRelperm(xmlFile1Path)
    phaseVisc = getViscosity(pvdgFilePath, pvtwFilePath)
    phaseDens = getDensity(xmlFile1Path)
    area = (yMax - yMin) * (zMax - zMin)
    q = getInjectionRate(xmlFile1Path) / phaseDens[0]

    # Read simulation output from HDF5 file
    hf = h5py.File(hdf5FilePath, 'r')
    time = hf.get('phaseVolumeFraction Time')
    time = np.asarray(time)
    center = hf.get('phaseVolumeFraction elementCenter')
    center = np.asarray(center)
    phaseVolFracFromGEOSX = hf.get('phaseVolumeFraction')
    phaseVolFracFromGEOSX = np.asarray(phaseVolFracFromGEOSX)

    n = 100000    # sampling value, will decide the accuracy of front detection in the BL analytical solution
    minDiff = 1e99
    iShock = n - 1
    volFrac = np.linspace(0, 1, num=n)
    dFractionalFlow_dS = np.zeros((n, ))

    for i in range(0, len(volFrac), 1):

        # step 1: compute the fractional flow and its derivative
        fractionalFlow, \
        dFractionalFlow_dS[i] = computeFractionalFlow( [ volFrac[i], (1-volFrac[i]) ], \
                                                       phaseMinVolFrac, \
                                                       phaseRelPermExp, \
                                                       phaseRelPermMaxValue, \
                                                       phaseVisc )

        # step 2: compute the location of the shock by matching the slope of the tangent to dfw_dS
        if (volFrac[i] > 0 and dFractionalFlow_dS[i]):
            slopeFromOrigin = fractionalFlow / volFrac[i]
            if (abs((slopeFromOrigin - dFractionalFlow_dS[i]) / dFractionalFlow_dS[i]) < minDiff):
                minDiff = abs((slopeFromOrigin - dFractionalFlow_dS[i]) / dFractionalFlow_dS[i])
                iShock = i

    # step 3: beyond the shock, constant dfw_dS as a function of saturation
    for i in range(0, iShock, 1):
        dFractionalFlow_dS[i] = dFractionalFlow_dS[iShock]

    fsize = 30
    lw = 6
    fig, ax = plt.subplots(1, 2, figsize=(24, 10))
    cmap = plt.get_cmap("tab10")
    iplt = -1
    # loop over the report times
    for iRptTime in range(1, time.shape[0], 1):

        iplt += 1

        # plot analytical solution
        td = q / area * (time[iRptTime, 0]) / porosity / xMax
        ax[0].plot( (time[iRptTime,0])*q*dFractionalFlow_dS/(area*porosity)/xMax, \
                  volFrac, \
                  lw=lw, linestyle='-', alpha=0.5, \
                  color=cmap(iplt), label='Analytical_t*='+str(round(td, 3)) )

        # plot numerical solution
        ax[0].plot( center[0,:,0]/xMax, \
                  phaseVolFracFromGEOSX[iRptTime,:,0], \
                  lw=lw, linestyle=':', \
                  color=cmap(iplt), label='Numerical_t*='+str(round(td, 3)) )

        # plot analytical solution
        ax[1].plot( (time[iRptTime,0])*q*dFractionalFlow_dS/(area*porosity)/xMax, \
                  1.0-volFrac, \
                  lw=lw, linestyle='-', alpha=0.5, \
                  color=cmap(iplt), label='Analytical_t*='+str(round(td, 3)) )

        # plot numerical solution
        ax[1].plot( center[0,:,0]/xMax, \
                  phaseVolFracFromGEOSX[iRptTime,:,1], \
                  lw=lw, linestyle=':', \
                  color=cmap(iplt), label='Numerical_t*='+str(round(td, 3)) )

    ax[0].set_xlim((0, 1))
    ax[0].set_ylim((0, 1))
    ax[0].set_xlabel('$x_d$', size=fsize, weight="bold")
    ax[0].set_ylabel('$S_g$', size=fsize, weight="bold")
    ax[0].legend(loc='upper right', fontsize=fsize * 0.5)
    ax[0].grid(True, which="both", ls="-")
    ax[0].xaxis.set_tick_params(labelsize=fsize)
    ax[0].yaxis.set_tick_params(labelsize=fsize)

    ax[1].set_xlim((0, 1))
    ax[1].set_ylim((0, 1))
    ax[1].set_xlabel('$x_d$', size=fsize, weight="bold")
    ax[1].set_ylabel('$S_w$', size=fsize, weight="bold")
    ax[1].legend(loc='lower right', fontsize=fsize * 0.5)
    ax[1].grid(True, which="both", ls="-")
    ax[1].xaxis.set_tick_params(labelsize=fsize)
    ax[1].yaxis.set_tick_params(labelsize=fsize)

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.3, hspace=0.3)

    plt.show()


if __name__ == "__main__":
    main()
