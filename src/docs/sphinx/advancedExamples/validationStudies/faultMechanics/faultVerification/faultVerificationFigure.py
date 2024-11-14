import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import xml.etree.ElementTree as ElementTree
import csv
import os
import argparse


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

    E = float(param1.get("defaultYoungModulus"))
    nu = float(param1.get("defaultPoissonRatio"))
    K = E / 3.0 / (1.0 - 2.0 * nu)
    G = E / 2.0 / (1.0 + nu)
    hydromechanicalParameters["bulkModulus"] = K
    hydromechanicalParameters["shearModulus"] = G
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


def getCompressiveStressFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param = tree.findall('FieldSpecifications/FieldSpecification')
    Stress = np.empty(3)
    for elem in param:
        if elem.get("name") == "stressXX" and elem.get("component") == "0":
            Stress[0] = float(elem.get("scale"))
        elif elem.get("name") == "stressYY" and elem.get("component") == "1":
            Stress[1] = float(elem.get("scale"))
        elif elem.get("name") == "stressZZ" and elem.get("component") == "2":
            Stress[2] = float(elem.get("scale"))
        elif elem.get("name") == "initialPressure" and elem.get("initialCondition") == "1":
            Pr_i = float(elem.get("scale"))

    return Stress, Pr_i


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
    xmlFilePath = geosDir + "/inputFiles/poromechanics/faultPoroelastic_base.xml"

    # Extract info from XML
    hydromechanicalParameters = getHydromechanicalParametersFromXML(xmlFilePath)
    bBiot = hydromechanicalParameters["biotCoefficient"]
    Stress, Pr_i = getCompressiveStressFromXML(xmlFilePath)

    # Load simulation result for case with impermeable fault
    file = open(outputDir + "/result_imp.csv")
    csvreader = csv.reader(file)
    header = next(csvreader)
    rows = []
    for row in csvreader:
        rows.append(row)
    file.close()

    rows = np.array(rows)
    sxx_imp = np.empty(len(rows[:, 23]))
    syy_imp = np.empty(len(rows[:, 23]))
    sxy_imp = np.empty(len(rows[:, 23]))
    pp_imp = np.empty(len(rows[:, 23]))
    y_imp = np.empty(len(rows[:, 23]))
    for i in range(0, len(rows[:, 23])):
        pp_imp[i] = float(rows[i, 22])
        sxx_imp[i] = ((float(rows[i, 14]) - bBiot * pp_imp[i]) - (Stress[0] - bBiot * Pr_i)) / 1.0e6
        syy_imp[i] = ((float(rows[i, 15]) - bBiot * pp_imp[i]) - (Stress[1] - bBiot * Pr_i)) / 1.0e6
        sxy_imp[i] = float(rows[i, 19]) / 1.0e6
        y_imp[i] = float(rows[i, 26])

    # Load simulation result for case with permeable fault
    file = open( outputDir + "/result_per.csv")
    csvreader = csv.reader(file)
    header = next(csvreader)
    rows = []
    for row in csvreader:
        rows.append(row)
    file.close()

    rows = np.array(rows)
    sxx_per = np.empty(len(rows[:, 23]))
    syy_per = np.empty(len(rows[:, 23]))
    sxy_per = np.empty(len(rows[:, 23]))
    pp_per = np.empty(len(rows[:, 23]))
    y_per = np.empty(len(rows[:, 23]))
    for i in range(0, len(rows[:, 23])):
        pp_per[i] = float(rows[i, 22])
        sxx_per[i] = ((float(rows[i, 14]) - bBiot * pp_per[i]) - (Stress[0] - bBiot * Pr_i)) / 1.0e6
        syy_per[i] = ((float(rows[i, 15]) - bBiot * pp_per[i]) - (Stress[1] - bBiot * Pr_i)) / 1.0e6
        sxy_per[i] = float(rows[i, 19]) / 1.0e6
        y_per[i] = float(rows[i, 26])

    # Load analytical solution
    y_ana, sxx_per_ana, syy_per_ana, sxy_per_ana, sxx_imp_ana, syy_imp_ana, sxy_imp_ana = np.loadtxt(
        "AnalyticalSolution.txt", skiprows=1, unpack=True)

    #Visulization
    N1 = 1
    fsize = 32
    msize = 12
    lw = 8
    malpha = 0.6
    lalpha = 0.6

    fig, ax = plt.subplots(1, 3, figsize=(32, 16))
    cmap = plt.get_cmap("tab10")

    ax[0].plot(sxx_imp_ana, y_ana, color=cmap(1), lw=lw, alpha=lalpha, label='Analytical_imp')
    ax[0].plot(sxx_imp[1::N1], y_imp[1::N1], 'o', color=cmap(1), markersize=msize, alpha=malpha, label='GEOSX_imp')
    ax[0].plot(sxx_per_ana, y_ana, color=cmap(2), lw=lw, alpha=lalpha, label='Analytical_per')
    ax[0].plot(sxx_per[1::N1], y_per[1::N1], 'o', color=cmap(2), markersize=msize, alpha=malpha, label='GEOSX_per')
    ax[0].set_xlim(-22.0, 8.00)
    ax[0].set_ylim(-300.0, 300.0)
    ax[0].xaxis.set_major_locator(MultipleLocator(10))
    ax[0].set_xlabel(r'$\Delta \sigma_{xx}$ (MPa)', size=fsize, weight="bold")
    ax[0].set_ylabel(r'y (m)', size=fsize, weight="bold")
    ax[0].legend(loc='upper left', fontsize=fsize * 0.6)
    ax[0].grid(True)
    ax[0].xaxis.set_tick_params(labelsize=fsize)
    ax[0].yaxis.set_tick_params(labelsize=fsize)

    ax[1].plot(syy_imp_ana, y_ana, color=cmap(1), lw=lw, alpha=lalpha, label='Analytical_imp')
    ax[1].plot(syy_imp[1::N1], y_imp[1::N1], 'o', color=cmap(1), markersize=msize, alpha=malpha, label='GEOSX_imp')
    ax[1].plot(syy_per_ana, y_ana, color=cmap(2), lw=lw, alpha=lalpha, label='Analytical_per')
    ax[1].plot(syy_per[1::N1], y_per[1::N1], 'o', color=cmap(2), markersize=msize, alpha=malpha, label='GEOSX_per')
    ax[1].set_xlim(-12.0, 12.00)
    ax[1].set_ylim(-300.0, 300.0)
    ax[1].xaxis.set_major_locator(MultipleLocator(4))
    ax[1].set_xlabel(r'$\Delta \sigma_{yy}$ (MPa)', size=fsize, weight="bold")
    ax[1].set_ylabel(r'y (m)', size=fsize, weight="bold")
    #ax[1].legend(loc='upper right',fontsize=fsize*0.6)
    ax[1].grid(True)
    ax[1].xaxis.set_tick_params(labelsize=fsize)
    ax[1].yaxis.set_tick_params(labelsize=fsize)

    ax[2].plot(sxy_imp_ana, y_ana, color=cmap(1), lw=lw, alpha=lalpha, label='Analytical_imp')
    ax[2].plot(sxy_imp[1::N1], y_imp[1::N1], 'o', color=cmap(1), markersize=msize, alpha=malpha, label='GEOSX_imp')
    ax[2].plot(sxy_per_ana, y_ana, color=cmap(2), lw=lw, alpha=lalpha, label='Analytical_per')
    ax[2].plot(sxy_per[1::N1], y_per[1::N1], 'o', color=cmap(2), markersize=msize, alpha=malpha, label='GEOSX_per')
    ax[2].set_xlim(-14.0, 14.00)
    ax[2].set_ylim(-300.0, 300.0)
    ax[2].xaxis.set_major_locator(MultipleLocator(7))
    ax[2].set_xlabel(r'$\Delta \sigma_{xy}$ (MPa)', size=fsize, weight="bold")
    ax[2].set_ylabel(r'y (m)', size=fsize, weight="bold")
    #ax[2].legend(loc='upper right',fontsize=fsize*0.6)
    ax[2].grid(True)
    ax[2].xaxis.set_tick_params(labelsize=fsize)
    ax[2].yaxis.set_tick_params(labelsize=fsize)

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

    plt.show()


if __name__ == "__main__":
    main()
