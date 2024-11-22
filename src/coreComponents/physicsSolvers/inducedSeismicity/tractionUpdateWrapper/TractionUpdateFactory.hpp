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

#ifndef GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_TRACTIONUPDATEBUILDER_HPP
#define GEOS_PHYSICSSOLVERS_INDUCED_SEISMICITY_TRACTIONUPDATEBUILDER_HPP

#include "SpringSliderTractionUpdate.hpp"
#include "OneWayCoupledTractionUpdate.hpp"

namespace geos
{

namespace inducedSeismicity
{  

enum class TractionUpdateType : integer
{
    SpringSlider,
    OneWayCoupled
};

std::unique_ptr< FaultTractionUpdateBase > tractionUpdateFactory ( TractionUpdateType const & tractionUpdateType )
{
    switch (tractionUpdateType)
    {
    case TractionUpdateType::SpringSlider:
        return std::make_unique<SpringSliderTractionUpdate>();
        break;
    case TractionUpdateType::OneWayCoupled:
        return std::make_unique<OneWayCoupledTractionUpdate>();    
    default:
        GEOS_ERROR( GEOS_FMT( "Traction update type {} not recognized", tractionUpdateType ) );
        return nullptr;
        break;
    } 
}

ENUM_STRINGS( TractionUpdateType,
              "SpringSlider",
              "OneWayCoupled" );
}
}