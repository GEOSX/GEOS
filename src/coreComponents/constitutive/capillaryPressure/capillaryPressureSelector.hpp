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
 * @file capillaryPressureSelector.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP
#define GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/capillaryPressure/BrooksCoreyCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/JFunctionCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/TableCapillaryPressure.hpp"
#include "constitutive/capillaryPressure/VanGenuchtenCapillaryPressure.hpp"

namespace geos
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( CapillaryPressureBase const & capPres,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< BrooksCoreyCapillaryPressure,
                               JFunctionCapillaryPressure,
                               TableCapillaryPressure,
                               VanGenuchtenCapillaryPressure >::execute( capPres, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( CapillaryPressureBase & capPres,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< BrooksCoreyCapillaryPressure,
                               JFunctionCapillaryPressure,
                               TableCapillaryPressure,
                               VanGenuchtenCapillaryPressure >::execute( capPres, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_CAPILLARYPRESSURE_CAPILLARYPRESSURESELECTOR_HPP
