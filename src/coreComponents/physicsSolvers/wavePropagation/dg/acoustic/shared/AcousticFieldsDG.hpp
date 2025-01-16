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
 * @file AcousticFieldsDG.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIELDSDG_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICFIELDSDG_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"


namespace geos
{

namespace fields
{

namespace acousticfieldsdg
{

DECLARE_FIELD( Pressure_nm1,
               "pressure_nm1",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n-1." );

DECLARE_FIELD( Pressure_n,
               "pressure_n",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n." );


DECLARE_FIELD( Pressure_np1,
               "pressure_np1",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );

//DECLARE_FIELD( ForcingRHS,
//               "rhs",
//               array1d< real32 >,
//               0,
//               NOPLOT,
//               WRITE_AND_READ,
//               "RHS" );

//DECLARE_FIELD( AcousticMassMatrix,
//               "acousticMassVector",
//               array2d< real32 >,
//               0,
//               NOPLOT,
//               WRITE_AND_READ,
//               "Diagonal of the Mass Matrix." );
//
//DECLARE_FIELD( StiffnessVector,
//               "stiffnessVector",
//               array2d< real32 >,
//               0,
//               NOPLOT,
//               WRITE_AND_READ,
//               "Stiffness vector contains R_h*Pressure_n." );

//DECLARE_FIELD( DampingVector,
//               "dampingVector",
//               array1d< real32 >,
//               0,
//               NOPLOT,
//               WRITE_AND_READ,
//               "Diagonal of the Damping Matrix." );

DECLARE_FIELD( AcousticVelocity,
               "acousticVelocity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium velocity of the cell" );

DECLARE_FIELD( AcousticDensity,
               "acousticDensity",
               array1d< real32 >,
               1,
               NOPLOT,
               WRITE_AND_READ,
               "Medium density of the cell" );

DECLARE_FIELD( AcousticFreeSurfaceFaceIndicator,
               "acousticFreeSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

DECLARE_FIELD( ElementToOpposite,
               "elementToOpposite",
               array2d< localIndex >,
               -1,
               NOPLOT,
               WRITE_AND_READ,
               "Map from elements to the neighbor opposite to each vertex. -1 for boundary, -2 for free surface" );

DECLARE_FIELD( ElementToOppositePermutation,
               "elementToOppositePermutation",
               array2d< integer >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Map from elements to the permutation of the neighboring element, opposite to each vertex." );

DECLARE_FIELD( CharacteristicSize
               "charactersticSize",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Characteristic size of every given element, used for penalty term computation. Often this is just the radius of the inscribed sphere.");

DECLARE_FIELD( MassPlusDampingInvIndex,
               "massPlusDampingInvIndex",
               array1d< localIndex >,
               -1,
               NOPLOT,
               WRITE_AND_READ,
               "Index in the list of the pre-computed mass+damping inverses, or -1 if not a boundary element" );

//DECLARE_FIELD( AcousticFreeSurfaceNodeIndicator,
//               "acousticFreeSurfaceNodeIndicator",
//               array1d< localIndex >,
//               0,
//               NOPLOT,
//               WRITE_AND_READ,
//               "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_HPP_ACOUSTICFIELDSDG */
