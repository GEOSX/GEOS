/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file multiFluidSelector.hpp
 */
#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/fluid/multifluid/blackOil/DeadOilFluid.hpp"
#include "constitutive/fluid/multifluid/blackOil/BlackOilFluid.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/CO2BrineFluid.hpp"

#include "common/GeosxConfig.hpp"
#ifdef GEOS_USE_PVTPackage
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluidPVTPackage.hpp"
#endif
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"

namespace geos
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiFluidBase const & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
#ifdef GEOS_USE_PVTPackage
                               CompositionalMultiphaseFluidPVTPackage,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid,
// Including these in a CUDA build will lead to compiler segfault.
// Need to split compilation units for all the options
#if !defined(GEOS_DEVICE_COMPILE)
                               CO2BrineEzrokhiThermalFluid,
                               CompositionalTwoPhaseLohrenzBrayClarkViscosity,
#endif
                               CompositionalTwoPhaseConstantViscosity
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( MultiFluidBase & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
#ifdef GEOS_USE_PVTPackage
                               CompositionalMultiphaseFluidPVTPackage,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid,
// Including these in a CUDA build will lead to compiler segfault.
// Need to split compilation units for all the options"
#if !defined(GEOS_DEVICE_COMPILE)
                               CO2BrineEzrokhiThermalFluid,
                               CompositionalTwoPhaseLohrenzBrayClarkViscosity,
#endif
                               CompositionalTwoPhaseConstantViscosity
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
