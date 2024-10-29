import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import HydrofractureSolutions
import sys


def getParametersFromXML(xmlFilePath):
    prefix = "../../../../../../../inputFiles/hydraulicFracturing/"

    tree = ElementTree.parse(prefix + xmlFilePath + "_benchmark.xml")

    maxTime = float(tree.find('Events').get('maxTime'))

    param = tree.findall('Geometry/Box')
    found_source = False
    for elem in param:
        if elem.get("name") == "source":
            source = elem.get("xMax")
            source = [float(i) for i in source[1:-1].split(",")]
            x_source = round(source[0])
            found_source = True
        if found_source: break

    tree = ElementTree.parse(prefix + xmlFilePath + "_base.xml")

    elasticParam = tree.find('Constitutive/ElasticIsotropic')    

    viscosity = float(tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultViscosity'))    

    fluidDensity = float(tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultDensity'))

    if xmlFilePath == 'kgdToughnessDominated' or xmlFilePath == 'kgdViscosityDominated':

        youngModulus = float(elasticParam.get('defaultYoungModulus'))
        poissonRatio = float(elasticParam.get('defaultPoissonRatio'))
        injectionRate = -2.0 * float(tree.find('FieldSpecifications/SourceFlux').get('scale')) / fluidDensity
        tree = ElementTree.parse(prefix + xmlFilePath + "_base.xml")
        toughness = float(tree.find('Solvers/SurfaceGenerator').get('rockToughness'))

    elif xmlFilePath == 'pennyShapedToughnessDominated' or xmlFilePath == 'pennyShapedViscosityDominated':
        K = float(elasticParam.get('defaultBulkModulus'))
        G = float(elasticParam.get('defaultShearModulus'))
        youngModulus = (9.0 * K * G) / (3.0 * K + G)
        poissonRatio = youngModulus / (2.0 * G) - 1.0
        injectionRate = -4.0 * float(tree.find('FieldSpecifications/SourceFlux').get('scale')) / fluidDensity
        tree = ElementTree.parse(prefix + xmlFilePath + "_benchmark.xml")
        toughness = float(tree.find('Solvers/SurfaceGenerator').get('rockToughness'))

    elif xmlFilePath == 'pknViscosityDominated':
        K = float(elasticParam.get('defaultBulkModulus'))
        G = float(elasticParam.get('defaultShearModulus'))
        youngModulus = (9.0 * K * G) / (3.0 * K + G)
        poissonRatio = youngModulus / (2.0 * G) - 1.0
        injectionRate = -4.0 * float(tree.find('FieldSpecifications/SourceFlux').get('scale')) / fluidDensity
        tree = ElementTree.parse(prefix + xmlFilePath + "_benchmark.xml")
        toughness = float(tree.find('Solvers/SurfaceGenerator').get('rockToughness'))
        found_core = False
        for elem in param:
            if elem.get("name") == "core":
                source = elem.get("xMax")
                source = [float(i) for i in source[1:-1].split(",")]
                x_source = round(source[1]) * 2.0
                found_core = True
            if found_core: break

    return [maxTime, youngModulus, poissonRatio, toughness, viscosity, injectionRate, x_source]



def main(xmlFilePathPrefix=''):
    if not xmlFilePathPrefix:
        xmlFilePathPrefix = sys.argv[1]

    tMax, E, nu, KIC, mu, Q0, xSource = getParametersFromXML(xmlFilePathPrefix)
    Ep = E / (1.0 - nu**2.0)

    t = np.arange(0.01 * tMax, tMax, 0.01 * tMax)
    radTimes = np.array([tMax])

    if xmlFilePathPrefix == 'kgdToughnessDominated':
        hfsolns = HydrofractureSolutions.KGDSolutions()
        kgdFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
        inletPressure = kgdFrac[5]
        halfLength = kgdFrac[6]
        inletAperture = kgdFrac[7]
        lablelist = ['Asymptotic ( $\mu$ => 0, $C_{L}$ => 0 )', 'GEOSX ( $\mu$ => 0, $C_{L}$ => 0 )']

    elif xmlFilePathPrefix == 'kgdViscosityDominated':
        hfsolns = HydrofractureSolutions.KGDSolutions()
        kgdFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
        inletPressure = kgdFrac[8]
        halfLength = kgdFrac[9]
        inletAperture = kgdFrac[10]
        lablelist = ['Asymptotic ( $K_{IC}$ => 0, $C_{L}$ => 0 )', 'GEOSX ( $K_{IC}$ => 0, $C_{L}$ => 0 )']

    elif xmlFilePathPrefix == 'pennyShapedToughnessDominated':
        hfsolns = HydrofractureSolutions.PennySolutions()
        pennyFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
        inletPressure = pennyFrac[5]
        halfLength = pennyFrac[6]
        inletAperture = pennyFrac[7]
        lablelist = ['Asymptotic ( $\mu$ => 0, $C_{L}$ => 0 )', 'GEOSX ( $\mu$ => 0, $C_{L}$ => 0 )']

    elif xmlFilePathPrefix == 'pennyShapedViscosityDominated':
        hfsolns = HydrofractureSolutions.PennySolutions()
        pennyFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
        inletPressure = pennyFrac[8]
        halfLength = pennyFrac[9]
        inletAperture = pennyFrac[10]
        lablelist = ['Asymptotic ( $K_{IC}$ => 0, $C_{L}$ => 0 )', 'GEOSX ( $K_{IC}$ => 0, $C_{L}$ => 0 )']

    elif xmlFilePathPrefix == 'pknViscosityDominated':
        pknFrac = HydrofractureSolutions.PKN_viscosityStorageDominated(E, nu, KIC, mu, Q0, t, xSource)
        halfLength, inletAperture, inletPressure = pknFrac.analyticalSolution()
        lablelist = ['Asymptotic ( $K_{IC}$ => 0, $C_{L}$ => 0 )', 'GEOSX ( $K_{IC}$ => 0, $C_{L}$ => 0 )']

    # GEOSX results   
    t_geosx, inletPressure_geosx, inletAperture_geosx, halfLength_geosx = np.loadtxt("model-results.txt", skiprows=1, unpack=True)

    # Visulization
    N1 = 1
    fsize = 30
    msize = 12
    lw = 8
    mew = 2
    malpha = 1.0

    fig, ax = plt.subplots(2, 2, figsize=(24, 18))
    cmap = plt.get_cmap("tab10")

    ax[0, 0].plot(t, halfLength, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[0, 0].plot(t_geosx[0::N1], halfLength_geosx[0::N1], 'o', color=cmap(0), mec=cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha)  
    ax[0, 0].set_xlim([min(t_geosx), max(t_geosx)])  
    ax[0, 0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0, 0].set_ylabel(r'Fracture Half Length (m)', size=fsize, weight="bold")
    ax[0, 0].legend(loc='lower right', fontsize=fsize * 0.7)
    ax[0, 0].grid(True)
    ax[0, 0].xaxis.set_tick_params(labelsize=fsize)
    ax[0, 0].yaxis.set_tick_params(labelsize=fsize)

    ax[0, 1].plot(t, inletAperture * 1000, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[0, 1].plot(t_geosx[0::N1], inletAperture_geosx[0::N1] * 1e3, 'o', color=cmap(0), mec=cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha)   
    ax[0, 1].set_xlim([min(t_geosx), max(t_geosx)]) 
    ax[0, 1].set_ylim(0, 2.0)
    ax[0, 1].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0, 1].set_ylabel(r'Fracture Mouth Opening (mm)', size=fsize, weight="bold")    
    ax[0, 1].grid(True)
    ax[0, 1].xaxis.set_tick_params(labelsize=fsize)
    ax[0, 1].yaxis.set_tick_params(labelsize=fsize)

    ax[1, 0].plot(t, inletPressure / 1.0e6, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[1, 0].plot(t_geosx[0::N1], inletPressure_geosx[0::N1] / 1e6, 'o', color=cmap(0), mec=cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha)    
    ax[1, 0].set_xlim([min(t_geosx), max(t_geosx)])
    ax[1, 0].set_ylim(0, 1.5)
    ax[1, 0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[1, 0].set_ylabel(r'Net Pressure at Well (MPa)', size=fsize, weight="bold")    
    ax[1, 0].grid(True)
    ax[1, 0].xaxis.set_tick_params(labelsize=fsize)
    ax[1, 0].yaxis.set_tick_params(labelsize=fsize)

    ax[1, 1].axis('off')
   
    plt.show()


if __name__ == "__main__":
    main()
