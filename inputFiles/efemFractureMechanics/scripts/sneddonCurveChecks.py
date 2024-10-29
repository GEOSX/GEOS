import numpy as np
import os
import sys
import xml.etree.ElementTree as ElementTree
import matplotlib
import matplotlib.pyplot as plt


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

    return length, origin[1]


def sneddon_curve_check_solution(**kwargs):
    # Read HDF5
    localX = np.squeeze(kwargs['displacementJump elementCenter'])[:, 0]

    #-------- Extract info from XML
    xmlFilePath = "./Sneddon_embeddedFrac_base.xml"

    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
    appliedPressure = getFracturePressureFromXML(xmlFilePath)

    # Get length of the fracture
    xmlFilePath = "./Sneddon_embeddedFrac_benchmark.xml"
    length, originShift = getFractureLengthFromXML(xmlFilePath)

    localX = localX - originShift

    analyticalSolution = Sneddon(mechanicalParameters, length, appliedPressure)
    aperture_analytical = np.empty(len(localX))
    i = 0
    for xCell in localX:
        aperture_analytical[i] = analyticalSolution.computeAperture(xCell)
        i += 1

    dispJumpAnalytical = np.zeros(np.shape(kwargs['displacementJump']))
    dispJump = kwargs['displacementJump']

    dispJumpAnalytical[:, :, 0] = aperture_analytical

    return dispJumpAnalytical


def debug():
    #-------- EmbeddeFrac File path
    import hdf5_wrapper
    hdf5File1Path = "Output/displacementJump_embeddedFrac.hdf5"

    # Read HDF5
    data = hdf5_wrapper.hdf5_wrapper(hdf5File1Path).get_copy()
    jump = data['displacementJump']
    jump = np.array(jump)
    aperture_EmbeddeFrac = jump[0, :, 0]
    x = data['displacementJump elementCenter']
    loc_EmbeddeFrac = x[0, :, 0]

    #-------- Extract info from XML
    xmlFilePath = "../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_base.xml"

    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)
    appliedPressure = getFracturePressureFromXML(xmlFilePath)

    # Get length of the fracture
    xmlFilePath = "../inputFiles/efemFractureMechanics/Sneddon_embeddedFrac_benchmark.xml"
    length, originShift = getFractureLengthFromXML(xmlFilePath)

    print(length)

    loc_EmbeddeFrac = loc_EmbeddeFrac - originShift

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

    ax.grid()
    ax.set_xlabel('Fracture Length [m]', size=fsize, weight="bold")
    ax.set_ylabel('Fracture Aperture [mm]', size=fsize, weight="bold")
    ax.legend(bbox_to_anchor=(0.5, 0.2), loc='center', borderaxespad=0., fontsize=fsize)
    ax.xaxis.set_tick_params(labelsize=fsize)
    ax.yaxis.set_tick_params(labelsize=fsize)

    plt.savefig("sneddon.png")


if __name__ == "__main__":
    debug()
