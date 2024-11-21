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
 * @file WellFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace well
{

DECLARE_FIELD( energyPerforationFlux,
               "energyPerforationFlux",
               array1d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Energy perforation flux" );

DECLARE_FIELD( dEnergyPerforationFlux,
               "dEnergyPerforationFlux",
               array3d< real64 >,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of energy perforation flux with respect to pressure temperature and global component density" );


}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLFIELDS_HPP_
