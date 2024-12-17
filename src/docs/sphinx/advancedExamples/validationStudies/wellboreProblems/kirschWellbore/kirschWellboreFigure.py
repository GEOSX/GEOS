import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
import math
from math import sin, cos, tan, exp, atan, asin
import os
import argparse


class Analytical:

    def __init__(self, mechanicalParameters, rw, Stress, Pw, theta):
        K = mechanicalParameters["bulkModulus"]
        G = mechanicalParameters["shearModulus"]
        E = (9 * K * G) / (3 * K + G)
        nu = E / (2 * G) - 1

        self.rw = rw
        self.Stress = Stress
        self.Pw = Pw
        self.theta = theta
        self.G = G
        self.nu = nu

    def computeRadialStress(self, x):
        return (self.Stress[0] + self.Stress[1]) / 2. * (1. - (self.rw / x)**2) + (
            self.Stress[0] - self.Stress[1]) / 2. * (1. - 4. * (self.rw / x)**2 + 3. * (self.rw / x)**4) * cos(
                2. * self.theta) + self.Pw * (self.rw / x)**2

    def computeHoopStress(self, x):
        return (self.Stress[0] + self.Stress[1]) / 2. * (1. + (self.rw / x)**2) - (
            self.Stress[0] - self.Stress[1]) / 2. * (1. + 3. * (self.rw / x)**4) * cos(
                2. * self.theta) - self.Pw * (self.rw / x)**2

    def computeShearStress(self, x):
        return -(self.Stress[0] - self.Stress[1]) / 2. * (1. + 2. * (self.rw / x)**2 - 3. *
                                                          (self.rw / x)**4) * sin(2. * self.theta)

    def computeRadialDisp(self, x):
        return -(self.rw)**2 / x / 2. / self.G * ((self.Stress[0] + self.Stress[1]) / 2. +
                                                  (self.Stress[0] - self.Stress[1]) / 2. *
                                                  (4. * (1. - self.nu) -
                                                   (self.rw / x)**2) * cos(2. * self.theta) - self.Pw)

    def computeShearDisp(self, x):
        return (self.rw)**2 / x / 2. / self.G * (self.Stress[0] - self.Stress[1]) / 2. * (2. * (1. - 2. * self.nu) +
                                                                                          (self.rw / x)**2) * sin(
                                                                                              2. * self.theta)


def getMechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)
    param = tree.find('Constitutive/ElasticIsotropic')
    mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus"])
    mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
    mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
    return mechanicalParameters


def getCompressiveStressFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.findall('FieldSpecifications/FieldSpecification')
    Stress = np.empty(3)
    for elem in param:
        if elem.get("name") == "Sxx" and elem.get("component") == "0":
            Stress[0] = float(elem.get("scale")) * (-1)
        elif elem.get("name") == "Syy" and elem.get("component") == "1":
            Stress[1] = float(elem.get("scale")) * (-1)
        elif elem.get("name") == "Szz" and elem.get("component") == "2":
            Stress[2] = float(elem.get("scale")) * (-1)

    param = tree.findall('FieldSpecifications/Traction')
    found_stress = False
    for elem in param:
        if elem.get("name") == "WellLoad" and elem.get("tractionType") == "normal":
            Pw = float(elem.get("scale")) * (-1)
            found_stress = True
        if found_stress: break

    return Stress, Pw


def getGeometryFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    wellRadius = tree.find('Mesh/InternalWellbore')
    radius = wellRadius.get("radius")
    radius = [float(i) for i in radius[1:-1].split(",")]
    rw = radius[0]
    rout = radius[1]

    nt = wellRadius.get("nt")
    nt = [float(i) for i in nt[1:-1].split(",")]
    elenum_t = nt[0]

    return rw, rout, elenum_t


# Rotate stress tensor, counter-clockwise rotation direction defined as positive
def stressRotation(stress, phi_x):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])

    return np.dot(np.dot(rotx, stress), np.transpose(rotx))


# Rotate displacement vector, counter-clockwise rotation direction defined as positive
def dispRotation(disp, phi_x):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])

    return np.dot(rotx, disp)


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
    outputDir = args.outputDir
    geosDir = args.geosDir
    hdf5File1Path = outputDir + "/stress_history.hdf5"
    hdf5File2Path = outputDir + "/displacement_history.hdf5"
    xmlFile1Path = geosDir + "/inputFiles/solidMechanics/KirschProblem_base.xml"
    xmlFile2Path = geosDir + "/inputFiles/solidMechanics/KirschProblem_benchmark.xml"

    # Read HDF5
    # Global Coordinate of Element Center
    hf = h5py.File(hdf5File1Path, 'r')
    xl_elm = hf.get('rock_stress elementCenter')
    xl_elm = np.asarray(xl_elm)
    xcord_elm = xl_elm[0, :, 0]
    ycord_elm = xl_elm[0, :, 1]
    zcord_elm = xl_elm[0, :, 2]
    # Load Stress Components
    sigma = hf.get('rock_stress')
    sigma = np.asarray(sigma)
    sigma_Cart = np.zeros([len(sigma[0, :, 0]), 6])
    for i in range(0, len(sigma[0, :, 0])):
        for j in range(0, 6):
            for k in range(0, 8):
                sigma_Cart[i, j] += sigma[0, i, j + 6 * k] / 8.

    # Global Coordinate of Nodal Point
    hf = h5py.File(hdf5File2Path, 'r')
    xl_node = hf.get('totalDisplacement ReferencePosition')
    xl_node = np.asarray(xl_node)
    xcord_node = xl_node[0, :, 0]
    ycord_node = xl_node[0, :, 1]
    zcord_node = xl_node[0, :, 2]
    # Load Displacement Components
    disp_load = hf.get('totalDisplacement')
    disp_load = np.asarray(disp_load)
    disp_Cart = disp_load[0, :, :]

    # Extract Mechanical Properties and Fracture Geometry from XML
    mechanicalParameters = getMechanicalParametersFromXML(xmlFile1Path)
    Stress, Pw = getCompressiveStressFromXML(xmlFile1Path)
    rw, rout, nt = getGeometryFromXML(xmlFile2Path)

    # Extract Curve
    theta = (1. / 4. - 1. / 2. / nt) * np.pi
    rlist_elm = []
    sigmalist = []
    for i in range(0, len(zcord_elm)):
        if abs(zcord_elm[i] / 0.5 - 1.) < 0.01 and abs(xcord_elm[i] * tan(theta) / ycord_elm[i] - 1.) < 0.01:
            rlist_elm.append((xcord_elm[i]**2 + ycord_elm[i]**2)**0.5)
            sigmalist.append(sigma_Cart[i, :] / 1.0e6)

    sig_rr, sig_tt, sig_rt = [], [], []
    for i in range(0, len(rlist_elm)):
        stress = np.array([[sigmalist[i][0],sigmalist[i][5],sigmalist[i][4]],\
                           [sigmalist[i][5],sigmalist[i][1],sigmalist[i][3]],\
                           [sigmalist[i][4],sigmalist[i][3],sigmalist[i][2]]])

        stressLocal = stressRotation(stress, theta)
        sig_rr.append(stressLocal[0][0])
        sig_tt.append(stressLocal[1][1])
        sig_rt.append(stressLocal[0][1])

    theta = 1. / 4. * np.pi
    rlist_node = []
    displist = []
    for i in range(0, len(zcord_node)):
        if abs(zcord_node[i] / 1.0 - 1.) < 0.01 and abs(xcord_node[i] * tan(theta) /
                                                        (ycord_node[i] + 1e-6) - 1.) < 0.01:
            rlist_node.append((xcord_node[i]**2 + ycord_node[i]**2)**0.5)
            displist.append(disp_Cart[i, :] * 1.0e3)

    u_r, u_t = [], []
    for i in range(0, len(rlist_node)):
        disp = np.array([displist[i][0], displist[i][1], displist[i][2]])

        dispLocal = dispRotation(disp, theta)
        u_r.append(dispLocal[0])
        u_t.append(dispLocal[1])

    # Initialize analytical solution
    AnalyticalSolution = Analytical(mechanicalParameters, rw, Stress, Pw, theta)

    # Plot Analytical (continuous line) and Numerical (markers) Solution
    x_analytical = np.linspace(rw, rout, 1000, endpoint=True)
    srr_analytical = np.empty(len(x_analytical))
    s00_analytical = np.empty(len(x_analytical))
    sr0_analytical = np.empty(len(x_analytical))
    ur_analytical = np.empty(len(x_analytical))
    u0_analytical = np.empty(len(x_analytical))
    i = 0
    for xCell in x_analytical:
        srr_analytical[i] = AnalyticalSolution.computeRadialStress(xCell) / 1.0e6
        s00_analytical[i] = AnalyticalSolution.computeHoopStress(xCell) / 1.0e6
        sr0_analytical[i] = AnalyticalSolution.computeShearStress(xCell) / 1.0e6
        ur_analytical[i] = AnalyticalSolution.computeRadialDisp(xCell) * 1.0e3
        u0_analytical[i] = AnalyticalSolution.computeShearDisp(xCell) * 1.0e3
        i += 1

    #Visulization
    N1 = 1
    fsize = 32
    msize = 15
    lw = 8
    malpha = 1.0

    fig, ax = plt.subplots(1, 2, figsize=(32, 12))
    cmap = plt.get_cmap("tab10")

    ax[0].semilogx(x_analytical,
                   -srr_analytical,
                   lw=lw,
                   alpha=0.5,
                   color=cmap(0),
                   label=u"\u03C3"
                   r'$_{rr}$ - Analytical')
    ax[0].semilogx(rlist_elm,
                   sig_rr,
                   'o',
                   color=cmap(0),
                   markersize=msize,
                   alpha=malpha,
                   label=u"\u03C3"
                   r'$_{rr}$ - GEOSX')
    ax[0].semilogx(x_analytical,
                   -s00_analytical,
                   lw=lw,
                   alpha=0.5,
                   color=cmap(1),
                   label=u"\u03C3"
                   r'$_{\theta\theta}$ - Analytical')
    ax[0].semilogx(rlist_elm,
                   sig_tt,
                   'o',
                   color=cmap(1),
                   markersize=msize,
                   alpha=malpha,
                   label=u"\u03C3"
                   r'$_{\theta\theta}$ - GEOSX')
    ax[0].semilogx(x_analytical,
                   -sr0_analytical,
                   lw=lw,
                   alpha=0.5,
                   color=cmap(2),
                   label=u"\u03C3"
                   r'$_{r\theta}$ - Analytical')
    ax[0].semilogx(rlist_elm,
                   sig_rt,
                   'o',
                   color=cmap(2),
                   markersize=msize,
                   alpha=malpha,
                   label=u"\u03C3"
                   r'$_{r\theta}$ - GEOSX')
    ax[0].set_xlim(rw, rout)
    ax[0].set_xlabel(r'r (m)', size=fsize, weight="bold")
    ax[0].set_ylabel(u"\u03C3"
                     r' (MPa)', size=fsize, weight="bold")
    ax[0].legend(loc='lower right', fontsize=fsize * 0.8)
    ax[0].grid(True)
    ax[0].xaxis.set_tick_params(labelsize=fsize)
    ax[0].yaxis.set_tick_params(labelsize=fsize)

    ax[1].semilogx(x_analytical, ur_analytical, lw=lw, alpha=0.5, color=cmap(0), label=r'u$_{r}$ - Analytical')
    ax[1].semilogx(rlist_node, u_r, 'o', color=cmap(0), markersize=msize, alpha=malpha, label=r'u$_{r}$ - GEOSX')
    ax[1].semilogx(x_analytical, u0_analytical, lw=lw, alpha=0.5, color=cmap(1), label=r'u$_{\theta}$ - Analytical')
    ax[1].semilogx(rlist_node, u_t, 'o', color=cmap(1), markersize=msize, alpha=malpha, label=r'u$_{\theta}$ - GEOSX')
    ax[1].set_xlim(rw, rout)
    ax[1].set_xlabel(r'r (m)', size=fsize, weight="bold")
    ax[1].set_ylabel(r'Displacement (mm)', size=fsize, weight="bold")
    ax[1].legend(loc='lower right', fontsize=fsize * 0.8)
    ax[1].grid(True)
    ax[1].xaxis.set_tick_params(labelsize=fsize)
    ax[1].yaxis.set_tick_params(labelsize=fsize)

    plt.show()


if __name__ == "__main__":
    main()
