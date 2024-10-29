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
 * @file singlePhaseThermalConductivitySelector.hpp
 */

#ifndef GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP
#define GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/thermalConductivity/SinglePhaseConstantThermalConductivity.hpp"

namespace geos
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( SinglePhaseThermalConductivityBase const & thermalConductivity,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< SinglePhaseConstantThermalConductivity >::execute( thermalConductivity, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( SinglePhaseThermalConductivityBase & thermalConductivity,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< SinglePhaseConstantThermalConductivity >::execute( thermalConductivity, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_SINGLEPHASE_THERMALCONDUCTIVITY_THERMALCONDUCTIVITYSELECTOR_HPP
