import numpy as np
import os
import sys

search_root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(search_root)

import HydrofractureSolutions


def get_youngs_modulus_poisson_ratio(bulk_modulus, shear_modulus):
    E = 9.0 * bulk_modulus * shear_modulus / (3.0 * bulk_modulus + shear_modulus)
    nu = (1.5 * bulk_modulus - shear_modulus) / (3.0 + bulk_modulus * shear_modulus)
    return E, nu


def get_plane_strain_modulus(E, nu):
    return E / (1.0 - nu**2.0)


def kgd_toughness_dominated_solutions(**kwargs):
    E = 30.0e9
    nu = 0.25
    KIC = 1.0e6
    mu = 1.0e-6
    Q0 = -5e-2 * -2.0 / 1000.0
    xSource = 0.0

    Ep = get_plane_strain_modulus(E, nu)
    t = np.squeeze(kwargs['hydraulicAperture Time'][:, 0])
    radTimes = t[-1:]
    hfsolns = HydrofractureSolutions.KGDSolutions()
    kgdFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
    return kgdFrac[5], kgdFrac[7]


def kgd_toughness_dominated_pressure_curve(**kwargs):
    return kgd_toughness_dominated_solutions(**kwargs)[0]


def kgd_toughness_dominated_aperture_curve(**kwargs):
    return kgd_toughness_dominated_solutions(**kwargs)[1]


def kgd_viscosity_dominated_solutions(**kwargs):
    E = 30.0e9
    nu = 0.25
    KIC = 1.0e4
    mu = 1.0e-3
    Q0 = -5e-2 * -2.0 / 1000.0
    xSource = 0.0

    Ep = get_plane_strain_modulus(E, nu)
    t = np.squeeze(kwargs['hydraulicAperture Time'][:, 0])
    radTimes = t[-1:]
    hfsolns = HydrofractureSolutions.KGDSolutions()
    kgdFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
    return kgdFrac[8], kgdFrac[10]


def kgd_viscosity_dominated_pressure_curve(**kwargs):
    return kgd_viscosity_dominated_solutions(**kwargs)[0]


def kgd_viscosity_dominated_aperture_curve(**kwargs):
    return kgd_viscosity_dominated_solutions(**kwargs)[1]


def penny_shaped_toughness_dominated_solutions(**kwargs):
    bulk_modulus = 20.0e9
    shear_modulus = 12.0e9
    KIC = 3e6
    mu = 1.0e-6
    Q0 = -6.625 * -2.0 / 1000.0
    xSource = 0.0

    E, nu = get_youngs_modulus_poisson_ratio(bulk_modulus, shear_modulus)
    Ep = get_plane_strain_modulus(E, nu)
    t = np.squeeze(kwargs['hydraulicAperture Time'][:, 0])
    radTimes = t[-1:]
    hfsolns = HydrofractureSolutions.PennySolutions()
    pennyFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
    return pennyFrac[5], pennyFrac[7]


def penny_shaped_toughness_dominated_pressure_curve(**kwargs):
    return penny_shaped_toughness_dominated_solutions(**kwargs)[0]


def penny_shaped_toughness_dominated_aperture_curve(**kwargs):
    return penny_shaped_toughness_dominated_solutions(**kwargs)[1]


def penny_shaped_viscosity_dominated_solutions(**kwargs):
    bulk_modulus = 20.0e9
    shear_modulus = 12.0e9
    KIC = 0.3e6
    mu = 1.0e-3
    Q0 = -6.625 * -2.0 / 1000.0
    xSource = 0.0

    E, nu = get_youngs_modulus_poisson_ratio(bulk_modulus, shear_modulus)
    Ep = get_plane_strain_modulus(E, nu)
    t = np.squeeze(kwargs['hydraulicAperture Time'][:, 0])
    radTimes = t[-1:]
    hfsolns = HydrofractureSolutions.PennySolutions()
    pennyFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
    return pennyFrac[8], pennyFrac[10]


def penny_shaped_viscosity_dominated_pressure_curve(**kwargs):
    return penny_shaped_viscosity_dominated_solutions(**kwargs)[0]


def penny_shaped_viscosity_dominated_aperture_curve(**kwargs):
    return penny_shaped_viscosity_dominated_solutions(**kwargs)[1]


def pkn_viscosity_dominated_solutions(**kwargs):
    bulk_modulus = 20.0e9
    shear_modulus = 12.0e9
    KIC = 0.1e6
    mu = 1.0e-3
    Q0 = -6.625 * -2.0 / 1000.0
    xSource = 0.0
    height = 6.0

    E, nu = get_youngs_modulus_poisson_ratio(bulk_modulus, shear_modulus)
    t = np.squeeze(kwargs['hydraulicAperture Time'][:, 0])
    pknFrac = HydrofractureSolutions.PKN_viscosityStorageDominated(E, nu, KIC, mu, Q0, t, height)
    halfLength, inletAperture, inletPressure = pknFrac.analyticalSolution()
    return inletPressure, inletAperture


def pkn_viscosity_dominated_pressure_curve(**kwargs):
    return pkn_viscosity_dominated_solutions(**kwargs)[0]


def pkn_viscosity_dominated_aperture_curve(**kwargs):
    return pkn_viscosity_dominated_solutions(**kwargs)[1]
