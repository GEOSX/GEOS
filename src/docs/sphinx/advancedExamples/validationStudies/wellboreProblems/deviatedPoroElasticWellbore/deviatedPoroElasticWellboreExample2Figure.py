import sys
sys.path.append('../')
import numpy as np
import matplotlib.pyplot as plt
import wellboreAnalyticalSolutions as analytic
import xml.etree.ElementTree as ElementTree
import os
import argparse


# Rotate stress from local coordinates of an inclined borehole to the global coordinates system
# See the description in fig.1 in Abousleiman and Cui 1998
def stressRotation(stress, phi_x, phi_z):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])
    rotz = np.array([[np.cos(phi_z), 0., np.sin(phi_z)], [0., 1., 0.], [-np.sin(phi_z), 0., np.cos(phi_z)]])

    return np.dot(np.dot(np.transpose(rotz), np.dot(np.dot(np.transpose(rotx), stress), rotx)), rotz)


# Rotate stress from global coordinates system to the local coordinates of an inclined borehole
# See the description in fig.1 in Abousleiman and Cui 1998
def stressRotationInv(stress, phi_x, phi_z):
    rotx = np.array([[np.cos(phi_x), np.sin(phi_x), 0.], [-np.sin(phi_x), np.cos(phi_x), 0.], [0., 0., 1.]])
    rotz = np.array([[np.cos(phi_z), 0., np.sin(phi_z)], [0., 1., 0.], [-np.sin(phi_z), 0., np.cos(phi_z)]])

    return np.dot(np.dot(np.transpose(rotx), np.dot(np.dot(np.transpose(rotz), stress), rotz)), rotx)


def analyticalResults(t, ri, theta, phi_x, phi_z, E, nu, M, Ks, kappa, bBiot, nE, nnu, nkappa, pi, pw, p0, Shmax, Shmin,
                      Sv):

    # For inclined borehole, the in-situ stress must be rotated to the local coordinates of the borehole
    # The solutions of Abousleiman and Cui 1998 are restricted to the case where the borehole is oriented in the direction of the material anisotropy
    S = analytic.stressRotation(Shmax, Shmin, Sv, phi_x, phi_z)
    Sx = S[0][0]
    Sxy = S[0][1]
    Sxz = S[0][2]
    Sy = S[1][1]
    Syz = S[1][2]
    Sz = S[2][2]

    if (Sx != Sy):
        theta_r = 0.5 * np.arctan(2. * Sxy / (Sx - Sy))
    else:
        theta_r = 0.

    PP0 = (Sx + Sy) / 2.
    S0 = -(((Sx - Sy) / 2.)**2. + Sxy**2.)**0.5

    r = np.arange(ri, 10. * ri, 0.005 * ri)

    # Elastic stiffnesses: see Eq.2 in Abousleiman and Cui 1998
    E_p = E / nE
    nu_p = nu / nnu
    kappa_p = kappa / nkappa

    M11 = E * (E_p - E * nu_p**2.) / (1. + nu) / (E_p - E_p * nu - 2. * E * nu_p**2.)
    M12 = E * (E_p * nu + E * nu_p**2.) / (1. + nu) / (E_p - E_p * nu - 2. * E * nu_p**2.)
    M13 = E * E_p * nu_p / (E_p - E_p * nu - 2. * E * nu_p**2.)
    M33 = E_p**2. * (1. - nu) / (E_p - E_p * nu - 2. * E * nu_p**2.)
    M44 = E / 2. / (1. + nu)
    #M55 = G_p
    G = M44

    # Anisotropic Biot's coefficients
    alpha = 1. - (M11 + M12 + M13) / (3. * Ks)
    alpha_p = 1. - (2. * M13 + M33) / (3. * Ks)

    # Fluid diffusion coefficient
    c = kappa * M * M11 / (M11 + alpha**2. * M)

    p, sig_rr, sig_tt, sig_rt = analytic.inTime_mode123(t, r, ri, PP0, pw, p0, pi, S0, theta, theta_r, c, alpha, M, G,
                                                        M11, M12, kappa)
    #sig_zz,tau_rz,tau_tz = analytic.inTime_outPlane(r, ri,theta,p0,Sx,Sy,Sz,Sxz,Syz,sig_rr,sig_tt,p,nu_p,alpha,alpha_p)
    return [r, sig_rr / 1e6, sig_tt / 1e6, p / 1e6]


def getParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    elasticParam = tree.find('Constitutive/ElasticIsotropic')

    bulkModulus = float(elasticParam.get('defaultBulkModulus'))
    shearModulus = float(elasticParam.get('defaultShearModulus'))

    maxTime = float(tree.find('Events').get('maxTime'))

    fsParams = tree.findall('FieldSpecifications/FieldSpecification')
    for fsParam in fsParams:
        if ((fsParam.get('fieldName') == "pressure") & (fsParam.get('initialCondition') == "1")):
            p0 = float(fsParam.get('scale'))
        if ((fsParam.get('fieldName') == "rock_stress") & (fsParam.get('initialCondition') == "1") &
            (fsParam.get('component') == "0")):
            ShmaxEffective = float(fsParam.get('scale'))
        if ((fsParam.get('fieldName') == "rock_stress") & (fsParam.get('initialCondition') == "1") &
            (fsParam.get('component') == "1")):
            ShminEffective = float(fsParam.get('scale'))
        if ((fsParam.get('fieldName') == "rock_stress") & (fsParam.get('initialCondition') == "1") &
            (fsParam.get('component') == "2")):
            SvEffective = float(fsParam.get('scale'))
        if ((fsParam.get('fieldName') == "pressure") & (fsParam.get('initialCondition') != "1")):
            pi = float(fsParam.get('scale'))

    porosity = float(tree.find('Constitutive/BiotPorosity').get('defaultReferencePorosity'))

    skeletonBulkModulus = float(tree.find('Constitutive/BiotPorosity').get('defaultGrainBulkModulus'))
    fluidCompressibility = float(tree.find('Constitutive/CompressibleSinglePhaseFluid').get('compressibility'))

    bBiot = 1.0 - bulkModulus / skeletonBulkModulus
    MBiot = 1.0 / (porosity * fluidCompressibility + (bBiot - porosity) / skeletonBulkModulus)

    permParam = tree.find('Constitutive/ConstantPermeability').get('permeabilityComponents')
    permeability = float(permParam.replace('{', '').replace('}', '').strip().split(',')[0])

    viscosity = float(tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultViscosity'))

    return [
        maxTime, MBiot, bBiot, bulkModulus, shearModulus, permeability, viscosity, pi, p0, bBiot * p0 - ShmaxEffective,
        bBiot * p0 - ShminEffective, bBiot * p0 - SvEffective
    ]


def getWellboreGeometryFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    meshParam = tree.find('Mesh/InternalWellbore')
    radius = float(meshParam.get("radius").replace('{', '').replace('}', '').strip().split(',')[0])

    # Wellbore deviation
    trajectoryParam = tree.find('Mesh/InternalWellbore').get('trajectory').replace(' ', '').split('},')
    top = trajectoryParam[0].replace('{', '').replace('}', '').strip().split(',')
    bottom = trajectoryParam[1].replace('{', '').replace('}', '').strip().split(',')

    dx = float(top[0]) - float(bottom[0])
    dy = float(top[1]) - float(bottom[1])
    dz = float(top[2]) - float(bottom[2])
    dl = np.sqrt(dx * dx + dy * dy)

    phi_x = np.arctan(dy / dx)
    phi_z = np.arctan(-dl / dz)

    return [radius, phi_x, phi_z]


def main():

   # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Script to generate figure from tutorial.")

    # Add arguments to accept individual file paths
    parser.add_argument('--geosDir', help='Path to the GEOS repository ', default='../../../../../../..')
    parser.add_argument('--outputDir', help='Path to output directory', default='.')

    # Parse the command-line arguments
    args = parser.parse_args()

    outputDir = args.outputDir
    geosDir = args.geosDir

    xmlFilePathPrefix = geosDir + "/inputFiles/wellbore/DeviatedPoroElasticWellbore_Drilling"

    geometry = getWellboreGeometryFromXML(xmlFilePathPrefix + "_benchmark.xml")
    parameters = getParametersFromXML(xmlFilePathPrefix + "_base.xml")

    # Time
    t = parameters[0]

    # Geometry
    ri = geometry[0]    #borehole radius
    phi_x = geometry[1]    #Azimuth angle
    phi_z = geometry[2]    #Inclination angle

    # Tangent angle to x-axis where data are extracted, converted degree to radial
    theta = 90. * np.pi / 180.

    # Poroelastic properties
    M = parameters[1]    # Biot's modulus
    bBiot = parameters[2]    # Biot's coefficient
    K = parameters[3]
    G = parameters[4]
    permeability = parameters[5]
    fluidViscosity = parameters[6]

    E = 1.0 / (1.0 / 9.0 / K + 1.0 / 3.0 / G)    # in-plane Young's modulus
    nu = (3.0 * K - 2.0 * G) / (6.0 * K + 2.0 * G)    # in-plane Poisson's ratio
    Ks = K / (1.0 - bBiot)    # Bulk modulus of the solid phase
    kappa = permeability / fluidViscosity    # (m2/Pa/s) ratio between the intrinsic permeability and the dynamic viscosity of fluid

    # Anisotropy: three additional parameters are needed for a transversly isotropic problem
    # Out-plane shear modulus is not required because only stresses are calculated
    nE = 1.0    # =1 for the isotropic case
    nnu = 1.0    # =1 for the isotropic case
    nkappa = 1.0    # ratio between the in-plane and out-plane permeability, =1 for the isotropic case

    # Loading: in-situ stress, in-situ pore pressure, borehole pore pressure and mud pressure
    pi = parameters[7]    # borehole pore pressure
    pw = pi    # cake effect is ignored
    p0 = parameters[8]    # in-situ pore pressure

    Shmax = parameters[9]
    Shmin = parameters[10]
    Sv = parameters[11]

    r_anal, sig_rr_anal, sig_tt_anal, pPore_anal = analyticalResults(t, ri, -theta, phi_x, phi_z, E, nu, M, Ks, kappa,
                                                                     bBiot, nE, nnu, nkappa, pi, pw, p0, Shmax, Shmin,
                                                                     Sv)

    fig = plt.figure(figsize=[13, 10])

    # Get radial coordinate and compute analytical results
    r = []
    for line in open( outputDir + '/stress_11_drilling.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            rval = values[0]
            r.append(rval)

    # Get stress_ij and pore pressure
    # These data are extracted along the y-axis from the well center (theta angle = 90°)
    stress_11, stress_12, stress_13, stress_22, stress_23, stress_33, pPore = [], [], [], [], [], [], []
    for line in open( outputDir + '/stress_11_drilling.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            sigVal = values[1] * 1e-6    # convert to MPa
            stress_11.append(sigVal)

    for line in open( outputDir + '/stress_12_drilling.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            sigVal = values[1] * 1e-6    # convert to MPa
            stress_12.append(sigVal)

    for line in open( outputDir + '/stress_13_drilling.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            sigVal = values[1] * 1e-6    # convert to MPa
            stress_13.append(sigVal)

    for line in open( outputDir + '/stress_22_drilling.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            sigVal = values[1] * 1e-6    # convert to MPa
            stress_22.append(sigVal)

    for line in open( outputDir + '/stress_23_drilling.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            sigVal = values[1] * 1e-6    # convert to MPa
            stress_23.append(sigVal)

    for line in open( outputDir + '/stress_33_drilling.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            sigVal = values[1] * 1e-6    # convert to MPa
            stress_33.append(sigVal)

    for line in open( outputDir + '/pressure_drilling.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            pPoreVal = values[1] * 1e-6    # convert to MPa
            pPore.append(pPoreVal)

    #Compute sig_rr, sig_tt
    sig_rr, sig_tt = [], []
    for i in range(len(stress_11)):
        stress = np.array([[stress_11[i],stress_12[i],stress_13[i]],\
                           [stress_12[i],stress_22[i],stress_23[i]],\
                           [stress_13[i],stress_23[i],stress_33[i]]])

        stressLocal = stressRotationInv(stress, theta + phi_x, phi_z)
        sig_rr.append(stressLocal[0][0])
        sig_tt.append(stressLocal[1][1])

    plt.subplot(221)
    plt.plot(r, sig_rr, 'ko', label='GEOSX result')
    plt.plot(r_anal, sig_rr_anal + bBiot * pPore_anal, 'k', linewidth=2, label='Analytic')
    plt.ylabel('Effective radial stress (MPa)')
    plt.xlabel('r (m)')
    plt.xlim(ri, 10 * ri)
    plt.legend()

    plt.subplot(222)
    plt.plot(r, sig_tt, 'ko', label='GEOSX result')
    plt.plot(r_anal, sig_tt_anal + bBiot * pPore_anal, 'k', linewidth=2, label='Analytic')
    plt.ylabel('Effective tangent stress (MPa)')
    plt.xlabel('r (m)')
    plt.xlim(ri, 10 * ri)

    plt.subplot(223)
    plt.plot(r, pPore, 'ko', label='GEOSX result')
    plt.plot(r_anal, pPore_anal, 'k', linewidth=2, label='Analytic')
    plt.ylabel('Pore pressure (MPa)')
    plt.xlabel('r (m)')
    plt.xlim(ri, 10 * ri)

    plt.show()


if __name__ == "__main__":
    main()
