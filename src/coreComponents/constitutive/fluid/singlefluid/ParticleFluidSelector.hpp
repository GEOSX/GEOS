/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file particleFluidSelector.hpp
 */
#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PARTICLEFLUIDSELECTOR_HPP_
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PARTICLEFLUIDSELECTOR_HPP_

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/fluid/singlefluid/ParticleFluid.hpp"

namespace geos
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( ParticleFluidBase const & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< ParticleFluid >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( ParticleFluidBase & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< ParticleFluid >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_PARTICLEFLUIDSELECTOR_HPP_
