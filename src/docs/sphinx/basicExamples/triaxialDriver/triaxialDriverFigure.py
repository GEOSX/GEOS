import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin, cos, tan, exp, atan
import xml.etree.ElementTree as ElementTree
import os
import argparse


def yieldSurface(xmlFilePath, mechanicalParameters):
    tree = ElementTree.parse(xmlFilePath)

    model = tree.find('Tasks/TriaxialDriver')

    if model.get("material") == "DruckerPrager":

        param = tree.find('Constitutive/DruckerPrager')
        yieldParameters = dict.fromkeys(["p_Yield", "q_iniYield"])

        phi_i = mechanicalParameters["frictionAngle"]
        c_i = mechanicalParameters["cohesion"] / 1.0e6
        f_i = atan(6.0 * sin(phi_i / 180 * np.pi) / (3.0 - sin(phi_i / 180 * np.pi))) * 180 / np.pi
        d_i = 6.0 * c_i * cos(phi_i / 180 * np.pi) / (3.0 - sin(phi_i / 180 * np.pi))
        k_i = tan(f_i / 180 * np.pi)
        p_Yield = np.linspace(0, 50, 100)
        q_iniYield = k_i * p_Yield + d_i
        yieldParameters["p_Yield"] = p_Yield
        yieldParameters["q_iniYield"] = q_iniYield

    elif model.get("material") == "ExtendedDruckerPrager":

        param = tree.find('Constitutive/ExtendedDruckerPrager')
        yieldParameters = dict.fromkeys(["p_Yield", "q_iniYield", "q_resYield"])

        phi_i = mechanicalParameters["initialFrictionAngle"]
        phi_r = mechanicalParameters["residualFrictionAngle"]
        c_i = mechanicalParameters["cohesion"] / 1.0e6
        f_i = atan(6.0 * sin(phi_i / 180 * np.pi) / (3.0 - sin(phi_i / 180 * np.pi))) * 180 / np.pi
        f_r = atan(6.0 * sin(phi_r / 180 * np.pi) / (3.0 - sin(phi_r / 180 * np.pi))) * 180 / np.pi
        d_i = 6.0 * c_i * cos(phi_i / 180 * np.pi) / (3.0 - sin(phi_i / 180 * np.pi))
        po = d_i / tan(f_i / 180 * np.pi)
        d_r = po * tan(f_r / 180 * np.pi)

        k_i = tan(f_i / 180 * np.pi)
        p_Yield = np.linspace(0, 50, 100)
        q_iniYield = k_i * p_Yield + d_i
        k_r = tan(f_r / 180 * np.pi)
        q_resYield = k_r * p_Yield + d_r

        yieldParameters["p_Yield"] = p_Yield
        yieldParameters["q_iniYield"] = q_iniYield
        yieldParameters["q_resYield"] = q_resYield

    elif model.get("material") == "DelftEgg":

        param = tree.find('Constitutive/DelftEgg')
        yieldParameters = dict.fromkeys(["p_Yield", "q_iniYield", "p_CSL", "q_CSL"])

        pc0 = -mechanicalParameters["preConsolidationPressure"] / 1.0e6
        alpha = mechanicalParameters["shapeParameter"]
        M = mechanicalParameters["cslSlope"]

        p_CSL = np.linspace(0, pc0 * 2, 500)
        q_CSL = M * p_CSL
        qlist2 = np.zeros(len(p_CSL))
        for i in range(0, len(p_CSL)):
            if alpha**2 * p_CSL[i] * (2.0 * alpha * pc0 / (alpha + 1.0) -
                                      p_CSL[i]) - alpha**2 * (alpha - 1.0) / (alpha + 1.0) * pc0**2 < 0.0:
                qlist2[i] = 0.0
            else:
                qlist2[i] = M * pow(
                    alpha**2 * p_CSL[i] * (2.0 * alpha * pc0 / (alpha + 1.0) - p_CSL[i]) - alpha**2 * (alpha - 1.0) /
                    (alpha + 1.0) * pc0**2, 0.5)
        idx = np.argwhere(np.diff(np.sign(q_CSL - qlist2))).flatten()
        pc1 = p_CSL[idx[-1]] * 2
        plist_MCC = np.linspace(0, pc1 / 2.0, 500)
        qlist_MCC = M * pow(plist_MCC * (pc1 - plist_MCC), 0.5)
        plist_DE = np.linspace(pc1 / 2.0, pc0, 5000)
        qlist_DE = M * pow(
            alpha**2 * plist_DE * (2.0 * alpha * pc0 / (alpha + 1.0) - plist_DE) - alpha**2 * (alpha - 1.0) /
            (alpha + 1.0) * pc0**2, 0.5)
        p_Yield = np.concatenate((plist_MCC, plist_DE))
        q_iniYield = np.concatenate((qlist_MCC, qlist_DE))

        yieldParameters["p_Yield"] = p_Yield
        yieldParameters["q_iniYield"] = q_iniYield
        yieldParameters["p_CSL"] = p_CSL
        yieldParameters["q_CSL"] = q_CSL

    elif model.get("material") == "ModifiedCamClay":

        param = tree.find('Constitutive/DelftEgg')
        yieldParameters = dict.fromkeys(["p_Yield", "q_iniYield", "p_CSL", "q_CSL"])

        pc0 = -mechanicalParameters["preConsolidationPressure"] / 1.0e6
        alpha = 1.0
        M = mechanicalParameters["cslSlope"]

        p_CSL = np.linspace(0, pc0 * 2, 500)
        q_CSL = M * p_CSL
        p_Yield = np.linspace(0, pc0, 500)
        q_iniYield = M * pow(p_Yield * (pc0 - p_Yield), 0.5)

        yieldParameters["p_Yield"] = p_Yield
        yieldParameters["q_iniYield"] = q_iniYield
        yieldParameters["p_CSL"] = p_CSL
        yieldParameters["q_CSL"] = q_CSL

    return yieldParameters


def getMechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    model = tree.find('Tasks/TriaxialDriver')

    if model.get("material") == "DruckerPrager":
        param = tree.find('Constitutive/DruckerPrager')

        mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus", "cohesion", "frictionangle"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["cohesion"] = float(param.get("defaultCohesion"))
        mechanicalParameters["frictionAngle"] = float(param.get("defaultFrictionAngle"))

    elif model.get("material") == "ExtendedDruckerPrager":
        param = tree.find('Constitutive/ExtendedDruckerPrager')

        mechanicalParameters = dict.fromkeys(
            ["bulkModulus", "shearModulus", "cohesion", "initialFrictionAngle", "residualFrictionAngle"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["cohesion"] = float(param.get("defaultCohesion"))
        mechanicalParameters["initialFrictionAngle"] = float(param.get("defaultInitialFrictionAngle"))
        mechanicalParameters["residualFrictionAngle"] = float(param.get("defaultResidualFrictionAngle"))

    elif model.get("material") == "Elastic":
        param = tree.find('Constitutive/ElasticIsotropic')

        mechanicalParameters = dict.fromkeys(["bulkModulus", "shearModulus"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))

    elif model.get("material") == "DelftEgg":
        param = tree.find('Constitutive/DelftEgg')

        mechanicalParameters = dict.fromkeys(
            ["bulkModulus", "shearModulus", "preConsolidationPressure", "shapeParameter", "cslSlope"])
        mechanicalParameters["bulkModulus"] = float(param.get("defaultBulkModulus"))
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["preConsolidationPressure"] = float(param.get("defaultPreConsolidationPressure"))
        mechanicalParameters["shapeParameter"] = float(param.get("defaultShapeParameter"))
        mechanicalParameters["cslSlope"] = float(param.get("defaultCslSlope"))

    elif model.get("material") == "ModifiedCamClay":
        param = tree.find('Constitutive/ModifiedCamClay')

        mechanicalParameters = dict.fromkeys(
            ["shearModulus", "preConsolidationPressure", "cslSlope", "recompressionIndex"])
        mechanicalParameters["shearModulus"] = float(param.get("defaultShearModulus"))
        mechanicalParameters["preConsolidationPressure"] = float(param.get("defaultPreConsolidationPressure"))
        mechanicalParameters["cslSlope"] = float(param.get("defaultCslSlope"))
        mechanicalParameters["recompressionIndex"] = float(param.get("defaultRecompressionIndex"))

    return mechanicalParameters


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()

    # File path
    outputDir = args.outputDir
    geosDir = args.geosDir
    xmlFilePath = geosDir + "/inputFiles/triaxialDriver/triaxialDriver_ExtendedDruckerPrager_basicExample.xml"

    # Extract info from XML
    mechanicalParameters = getMechanicalParametersFromXML(xmlFilePath)

    # Load GEOSX results
    path = outputDir + "/simulationResults.txt"
    time, ax_strain, ra_strain1, ra_strain2, ax_stress, ra_stress1, ra_stress2, newton_iter, residual_norm = np.loadtxt(
        path, skiprows=5, unpack=True)
    p_num = -(ax_stress + 2.0 * ra_stress1) / 3.0 / 1.0e6
    q_num = -(ax_stress - ra_stress1) / 1.0e6
    strain_vol = ax_strain + 2.0 * ra_strain1

    #Visulization
    N1 = 1
    fsize = 30
    msize = 12
    lw = 6
    malpha = 0.5
    fig, ax = plt.subplots(3, 1, figsize=(15, 27))
    cmap = plt.get_cmap("tab10")

    ax[0].plot(-ax_strain * 100,
               q_num,
               'o',
               color=cmap(0),
               mec='b',
               markersize=msize,
               alpha=malpha,
               label='Triaxial Driver')
    ax[0].plot(-ra_strain1 * 100, q_num, 'o', color=cmap(0), mec='b', markersize=msize, alpha=malpha)
    ax[0].set_xlabel(r'Strain (%)', size=fsize, weight="bold")
    ax[0].set_ylabel(r'Deviatoric Stress (MPa)', size=fsize, weight="bold")
    ax[0].legend(loc='lower right', fontsize=fsize)
    ax[0].xaxis.set_tick_params(labelsize=fsize)
    ax[0].yaxis.set_tick_params(labelsize=fsize)

    ax[1].plot(p_num, q_num, 'o', color=cmap(0), mec='b', markersize=msize, alpha=malpha, label='Triaxial Driver')
    # Yield surface for plastic models
    tree = ElementTree.parse(xmlFilePath)
    model = tree.find('Tasks/TriaxialDriver')
    if model.get("material") == "DruckerPrager":
        yieldParameters = yieldSurface(xmlFilePath, mechanicalParameters)
        p_Yield = yieldParameters["p_Yield"]
        q_iniYield = yieldParameters["q_iniYield"]
        ax[1].plot(p_Yield, q_iniYield, lw=lw, alpha=0.8, color='k', linestyle='--', label='Initial Yield Surface')
    elif model.get("material") == "ExtendedDruckerPrager":
        yieldParameters = yieldSurface(xmlFilePath, mechanicalParameters)
        p_Yield = yieldParameters["p_Yield"]
        q_iniYield = yieldParameters["q_iniYield"]
        q_resYield = yieldParameters["q_resYield"]
        ax[1].plot(p_Yield, q_iniYield, lw=lw, alpha=0.8, color='k', linestyle='--', label='Initial Yield Surface')
        ax[1].plot(p_Yield,
                   q_resYield,
                   lw=lw,
                   alpha=0.8,
                   color='orange',
                   linestyle='--',
                   label='Residual Yield Surface')
    elif model.get("material") == "DelftEgg":
        yieldParameters = yieldSurface(xmlFilePath, mechanicalParameters)
        p_Yield = yieldParameters["p_Yield"]
        q_iniYield = yieldParameters["q_iniYield"]
        p_CSL = yieldParameters["p_CSL"]
        q_CSL = yieldParameters["q_CSL"]
        ax[1].plot(p_Yield, q_iniYield, lw=lw, alpha=0.8, color='k', linestyle='--', label='Initial Yield Surface')
        ax[1].plot(p_CSL, q_CSL, lw=lw, alpha=0.8, color='orange', linestyle='--', label='Critical State Line')
    elif model.get("material") == "ModifiedCamClay":
        yieldParameters = yieldSurface(xmlFilePath, mechanicalParameters)
        p_Yield = yieldParameters["p_Yield"]
        q_iniYield = yieldParameters["q_iniYield"]
        p_CSL = yieldParameters["p_CSL"]
        q_CSL = yieldParameters["q_CSL"]
        ax[1].plot(p_Yield, q_iniYield, lw=lw, alpha=0.8, color='k', linestyle='--', label='Initial Yield Surface')
        ax[1].plot(p_CSL, q_CSL, lw=lw, alpha=0.8, color='orange', linestyle='--', label='Critical State Line')

    ax[1].set_xlabel(r'p (MPa)', size=fsize, weight="bold")
    ax[1].set_ylabel(r'q (MPa)', size=fsize, weight="bold")
    ax[1].legend(loc='lower right', fontsize=fsize)
    ax[1].xaxis.set_tick_params(labelsize=fsize)
    ax[1].yaxis.set_tick_params(labelsize=fsize)

    ax[2].plot(-ax_strain * 100,
               -strain_vol * 100,
               'o',
               color=cmap(0),
               mec='b',
               markersize=msize,
               alpha=malpha,
               label='Triaxial Driver')
    ax[2].set_xlabel(r'Axial Strain (%)', size=fsize, weight="bold")
    ax[2].set_ylabel(r'Volumetric Strain (%)', size=fsize, weight="bold")
    ax[2].legend(loc='lower right', fontsize=fsize)
    ax[2].xaxis.set_tick_params(labelsize=fsize)
    ax[2].yaxis.set_tick_params(labelsize=fsize)

    plt.subplots_adjust(left=0.2, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

    plt.show()


if __name__ == "__main__":
    main()
