import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
from mpmath import *
import math
from math import sin, cos, tan, exp, atan, asin
from scipy.optimize import newton
import os
import argparse

class Mandel:

    def __init__(self, hydromechanicalParameters, alen, blen, appliedTraction):
        F = appliedTraction
        B = hydromechanicalParameters["skemptonCoefficient"]
        nu = hydromechanicalParameters["poissonRatio"]
        nuu = hydromechanicalParameters["undrainedPoissonRatio"]
        p0 = 1. / 3. / alen * B * (1. + nuu) * F

        alpha_n = []
        eps = np.finfo(np.float64).eps
        coef = (1 - nu) / (nuu - nu)

        n = 1
        while True:
            root = newton(func=lambda xi: tan(xi) - coef * xi,
                          x0=-math.pi / 2 + n * math.pi - 100 * eps,
                          fprime=lambda xi: 1 + pow(tan(xi), 2) - coef,
                          tol=eps)
            if root < 50:
                alpha_n.append(root)
                n += 1
            else:
                break
        self.alpha_n = alpha_n

        self.consolidationCoefficient = hydromechanicalParameters["consolidationCoefficient"]
        self.xlength = alen
        self.p0 = p0

        G = hydromechanicalParameters["shearModulus"]
        self.scaling1 = -F * (1.0 - nu) / 2.0 / G / alen
        self.scaling2 = F * (1.0 - nuu) / G / alen

    def computePressure(self, x, t):
        cc = self.consolidationCoefficient
        Lx = self.xlength
        alpha_n = self.alpha_n
        p0 = self.p0
        solution = 0
        for k in range(len(alpha_n)):
            an = alpha_n[k]
            solution = solution + 2. * sin(an) / (an - sin(an) * cos(an)) * (cos(an * x / Lx) - cos(an)) * exp(
                -an**2 * cc * t / Lx**2)

        return p0 * solution

    def computeVerticalDisplacement(self, z, t):
        cc = self.consolidationCoefficient
        Lx = self.xlength
        alpha_n = self.alpha_n
        scaling1 = self.scaling1
        scaling2 = self.scaling2
        solution = 0
        for k in range(len(alpha_n)):
            an = alpha_n[k]
            solution = solution + sin(an) * cos(an) / (an - sin(an) * cos(an)) * exp(-an**2 * cc * t / Lx**2)

        return (scaling1 + scaling2 * solution) * z


def getHydromechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param1 = tree.find('Constitutive/ElasticIsotropic')
    param2 = tree.find('Constitutive/BiotPorosity')
    param3 = tree.find('Constitutive/CompressibleSinglePhaseFluid')
    param4 = tree.find('Constitutive/ConstantPermeability')

    hydromechanicalParameters = dict.fromkeys([
        "bulkModulus", "shearModulus", "biotCoefficient", "fluidViscosity", "fluidCompressibility", "porosity",
        "permeability", "skemptonCoefficient", "poissonRatio", "undrainedPoissonRatio", "consolidationCoefficient"
    ])

    hydromechanicalParameters["bulkModulus"] = float(param1.get("defaultBulkModulus"))
    hydromechanicalParameters["shearModulus"] = float(param1.get("defaultShearModulus"))

    K = hydromechanicalParameters["bulkModulus"]
    G = hydromechanicalParameters["shearModulus"]
    E = (9.0 * K * G) / (3.0 * K + G)
    nu = E / (2.0 * G) - 1.0
    Ks = float(param2.get("defaultGrainBulkModulus"))

    hydromechanicalParameters["biotCoefficient"] = 1.0 - K / Ks
    hydromechanicalParameters["porosity"] = float(param2.get("defaultReferencePorosity"))
    hydromechanicalParameters["fluidViscosity"] = float(param3.get("defaultViscosity"))
    hydromechanicalParameters["fluidCompressibility"] = float(param3.get("compressibility"))

    perm = param4.get("permeabilityComponents")
    perm = np.array(perm[1:-1].split(','), float)
    hydromechanicalParameters["permeability"] = perm[0]

    phi = hydromechanicalParameters["porosity"]
    cf = hydromechanicalParameters["fluidCompressibility"]
    bBiot = hydromechanicalParameters["biotCoefficient"]
    kp = hydromechanicalParameters["permeability"]
    mu = hydromechanicalParameters["fluidViscosity"]
    M = 1. / (phi * cf + (bBiot - phi) / Ks)
    Ku = K + bBiot**2 * M
    B = bBiot * M / Ku
    nuu = (3. * nu + bBiot * B * (1 - 2. * nu)) / (3. - bBiot * B * (1 - 2. * nu))
    cc = 2. * kp / mu * B**2 * G * (1. - nu) * (1. + nuu)**2 / 9. / (1. - nuu) / (nuu - nu)
    hydromechanicalParameters["skemptonCoefficient"] = B
    hydromechanicalParameters["poissonRatio"] = nu
    hydromechanicalParameters["undrainedPoissonRatio"] = nuu
    hydromechanicalParameters["consolidationCoefficient"] = cc

    return hydromechanicalParameters


def getGeometryFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    meshElement = tree.find('Mesh/InternalMesh')
    dimensions = meshElement.get("xCoords")
    dimensions = [float(i) for i in dimensions[1:-1].split(",")]
    alen = dimensions[1]
    dimensions = meshElement.get("zCoords")
    dimensions = [float(i) for i in dimensions[1:-1].split(",")]
    blen = dimensions[1]
    dimensions = meshElement.get("nx")
    dimensions = [float(i) for i in dimensions[1:-1].split(",")]
    nx = dimensions[0]
    dimensions = meshElement.get("nz")
    dimensions = [float(i) for i in dimensions[1:-1].split(",")]
    nz = dimensions[0]

    return alen, blen, nx, nz


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
    hdf5File1Path = outputDir + "/pressure_history.hdf5"
    hdf5File2Path = outputDir + "/displacement_history.hdf5"
    xmlFile1Path = geosDir + "/inputFiles/poromechanics/PoroElastic_Mandel_base.xml"
    xmlFile2Path = geosDir + "/inputFiles/poromechanics/PoroElastic_Mandel_benchmark_fim.xml"

    # Read HDF5
    # Global Coordinate of Element Center
    hf = h5py.File(hdf5File1Path, 'r')
    xl = hf.get('pressure elementCenter')
    xcord = xl[0, :, 0]
    ycord = xl[0, :, 1]
    zcord = xl[0, :, 2]
    pl = hf.get('pressure')
    tl = hf.get('pressure Time')

    # Global Coordinate of Nodal Point
    hf = h5py.File(hdf5File2Path, 'r')
    xl_node = hf.get('totalDisplacement ReferencePosition')
    xcord_node = xl_node[0, :, 0]
    ycord_node = xl_node[0, :, 1]
    zcord_node = xl_node[0, :, 2]
    # Load Displacement Components
    disp = hf.get('totalDisplacement')

    # Extract Mechanical Properties and Fracture Geometry from XML
    hydromechanicalParameters = getHydromechanicalParametersFromXML(xmlFile1Path)
    F = 1.0e4
    La, Lb, na, nb = getGeometryFromXML(xmlFile2Path)
    B = hydromechanicalParameters["skemptonCoefficient"]
    nuu = hydromechanicalParameters["undrainedPoissonRatio"]
    G = hydromechanicalParameters["shearModulus"]
    p0 = 1. / 3. / La * B * (1. + nuu) * F
    u0 = -F * Lb * (1.0 - nuu) / 2. / G / La

    xd_numerical = xcord / La
    zd_numerical = zcord_node / Lb
    t = [0.05, 0.5, 5.0, 10.0]
    pressure_numerical = np.zeros([len(t), len(xd_numerical)])
    displacement_numerical = np.zeros([len(t), len(zd_numerical)])
    for i in range(len(t)):
        for j in range(1, len(tl)):
            if tl[j] <= t[i]:
                pressure_numerical[i, :] = pl[j - 1, :] / p0
                displacement_numerical[i, :] = disp[j - 1, :, 2] / u0

    # Initialize Mandel's analytical solution
    mandelAnalyticalSolution = Mandel(hydromechanicalParameters, La, Lb, F)

    x = np.linspace(0, La, 101, endpoint=True)
    xd_analytical = x / La
    pressure_analytical = np.zeros([len(t), len(x)])
    displacement_analytical = np.zeros([len(t), len(x)])
    for i in range(len(t)):
        for j in range(len(x)):
            pressure_analytical[i][j] = mandelAnalyticalSolution.computePressure(x[j], t[i]) / p0
            displacement_analytical[i][j] = mandelAnalyticalSolution.computeVerticalDisplacement(x[j], t[i]) / u0

    fsize = 30
    msize = 15
    lw = 8
    fig, ax = plt.subplots(1, 2, figsize=(32, 12))
    cmap = plt.get_cmap("tab10")

    for i in range(len(t)):
        ax[0].plot(xd_analytical,
                   pressure_analytical[i][:],
                   color=cmap(i),
                   alpha=0.6,
                   label='t = ' + str(t[i]) + 's - Analytical Solution',
                   lw=lw)
        ax[0].plot(xd_numerical,
                   pressure_numerical[i][:],
                   'o',
                   alpha=0.8,
                   color=cmap(i),
                   mec='k',
                   label='t = ' + str(t[i]) + 's - Numerical Solution',
                   markersize=msize)
    ax[0].set_xlim(0, 1)
    ax[0].set_ylim(0, 1.2)
    ax[0].set_xlabel('Normalized Distance, '
                     r'$x$/$a$', size=fsize, weight="bold")
    ax[0].set_ylabel('Normalized Pressure, '
                     r'$p$/$p_{0}$', size=fsize, weight="bold")
    ax[0].legend(bbox_to_anchor=(0.02, 0.6), loc='center left', borderaxespad=0., fontsize=fsize * 0.7)
    ax[0].grid(True)
    ax[0].xaxis.set_tick_params(labelsize=fsize)
    ax[0].yaxis.set_tick_params(labelsize=fsize)

    for i in range(len(t)):
        ax[1].plot(xd_analytical,
                   displacement_analytical[i][:],
                   color=cmap(i),
                   alpha=0.6,
                   label='t = ' + str(t[i]) + 's - Analytical Solution',
                   lw=lw)
        ax[1].plot(zd_numerical,
                   displacement_numerical[i][:],
                   'o',
                   alpha=0.8,
                   color=cmap(i),
                   mec='k',
                   label='t = ' + str(t[i]) + 's - Numerical Solution',
                   markersize=msize)
    ax[1].set_xlim(0, 1)
    ax[1].set_ylim(0, 1.5)
    ax[1].set_xlabel('Normalized Distance, '
                     r'$z$/$b$', size=fsize, weight="bold")
    ax[1].set_ylabel('Normalized Displacement, '
                     r'$u_{z}$/$u_{z=b,0}$', size=fsize, weight="bold")
    ax[1].legend(loc='upper left', fontsize=fsize * 0.7)
    ax[1].grid(True)
    ax[1].xaxis.set_tick_params(labelsize=fsize)
    ax[1].yaxis.set_tick_params(labelsize=fsize)

    plt.show()

if __name__ == "__main__":
    main()
