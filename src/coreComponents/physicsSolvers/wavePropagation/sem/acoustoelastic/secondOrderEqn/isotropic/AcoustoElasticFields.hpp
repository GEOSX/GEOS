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
 * @file AcoustoElasticFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTOELASTICFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTOELASTICFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"


namespace geos
{

namespace fields
{

namespace acoustoelasticfields
{

DECLARE_FIELD( CouplingVectorx,
               "couplingVectorx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Coupling term on x." );

DECLARE_FIELD( CouplingVectory,
               "couplingVectory",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Coupling term on y." );

DECLARE_FIELD( CouplingVectorz,
               "couplingVectorz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Coupling term on z." );

}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_ACOUSTOELASTICFIELDS */
