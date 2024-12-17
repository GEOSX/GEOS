/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file relativePermeabilitySelector.hpp
 */

#ifndef GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP
#define GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp"
#include "constitutive/relativePermeability/BrooksCoreyBakerRelativePermeability.hpp"
#include "constitutive/relativePermeability/TableRelativePermeability.hpp"
#include "constitutive/relativePermeability/TableRelativePermeabilityHysteresis.hpp"
#include "constitutive/relativePermeability/BrooksCoreyStone2RelativePermeability.hpp"
#include "constitutive/relativePermeability/VanGenuchtenBakerRelativePermeability.hpp"
#include "constitutive/relativePermeability/VanGenuchtenStone2RelativePermeability.hpp"

namespace geos
{

namespace constitutive
{

#ifndef PASSTHROUGH_HANDLE_CASE
#define PASSTHROUGH_HANDLE_CASE( MODEL ) \
  if( dynamicCast< MODEL * >( &relPerm ) ) \
  { \
    lambda( static_cast< MODEL & >( relPerm ) ); \
    return true; \
  }

template< typename LAMBDA >
void constitutiveUpdatePassThru( RelativePermeabilityBase const & relPerm,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< BrooksCoreyRelativePermeability,
                               BrooksCoreyBakerRelativePermeability,
                               TableRelativePermeability,
                               TableRelativePermeabilityHysteresis,
                               VanGenuchtenBakerRelativePermeability,
                               BrooksCoreyStone2RelativePermeability,
                               VanGenuchtenBakerRelativePermeability,
                               VanGenuchtenStone2RelativePermeability >::execute( relPerm, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( RelativePermeabilityBase & relPerm,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< BrooksCoreyRelativePermeability,
                               BrooksCoreyBakerRelativePermeability,
                               BrooksCoreyStone2RelativePermeability,
                               VanGenuchtenBakerRelativePermeability,
                               VanGenuchtenStone2RelativePermeability,
                               TableRelativePermeability,
                               TableRelativePermeabilityHysteresis,
                               VanGenuchtenBakerRelativePermeability >::execute( relPerm, std::forward< LAMBDA >( lambda ) );
}

#undef PASSTHROUGH_HANDLE_CASE
#endif

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP
