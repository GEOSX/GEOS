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
 * @file PoromechanicsFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSFIELDS_HPP_

#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace poromechanics
{

DECLARE_FIELD( bulkDensity,
               "bulkDensity",
               array2d< real64 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Bulk density" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSFIELDS_HPP_
