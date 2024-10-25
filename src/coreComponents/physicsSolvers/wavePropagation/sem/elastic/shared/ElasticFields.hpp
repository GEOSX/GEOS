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
 * @file ElasticFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{

namespace fields
{

namespace elasticfields
{

DECLARE_FIELD( Displacementx_nm1,
               "displacementx_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "x-component of displacement at time n-1." );

DECLARE_FIELD( Displacementy_nm1,
               "displacementy_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "y-component of displacement at time n-1." );

DECLARE_FIELD( Displacementz_nm1,
               "displacementz_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "z-component of displacement at time n-1." );

DECLARE_FIELD( Displacementx_n,
               "displacementx_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "x-component of displacement at time n." );

DECLARE_FIELD( Displacementy_n,
               "displacementy_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "y-component of displacement at time n." );

DECLARE_FIELD( Displacementz_n,
               "displacementz_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "z-component of displacement at time n." );

DECLARE_FIELD( Displacementx_np1,
               "displacementx_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "x-component of displacement at time n+1." );

DECLARE_FIELD( Displacementy_np1,
               "displacementy_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "y-component of displacement at time n+1." );

DECLARE_FIELD( Displacementz_np1,
               "displacementz_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "z-component of displacement at time n+1." );

DECLARE_FIELD( DivPsix,
               "divpsix",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "x-component of memory variables for attenuation." );

DECLARE_FIELD( DivPsiy,
               "divpsiy",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "y-component of memory variables for attenuation." );

DECLARE_FIELD( DivPsiz,
               "divpsiz",
               array2d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "z-component of memory variables for attenuation." );

DECLARE_FIELD( Stresstensorxx,
               "stresstensorxx",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xx-components of the stress tensor." );

DECLARE_FIELD( Stresstensoryy,
               "stresstensoryy",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "yy-components of the stress tensor." );

DECLARE_FIELD( Stresstensorzz,
               "stresstensorzz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "zz-components of the stress tensor." );

DECLARE_FIELD( Stresstensorxy,
               "stresstensorxy",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xy-components of the stress tensor (symetric of yx-component)." );

DECLARE_FIELD( Stresstensorxz,
               "stresstensorxz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "xz-components of the stress tensor (symetric of zx-component)." );

DECLARE_FIELD( Stresstensoryz,
               "stresstensoryz",
               array2d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "yz-components of the stress tensor (symetric of zy-component)." );

DECLARE_FIELD( ForcingRHS,
               "rhs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS" );

DECLARE_FIELD( ForcingRHSx,
               "rhsx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS for x-direction" );

DECLARE_FIELD( ForcingRHSy,
               "rhsy",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS for y-direction" );

DECLARE_FIELD( ForcingRHSz,
               "rhsz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "RHS for z-direction" );

DECLARE_FIELD( ElasticMassVector,
               "elasticMassVector",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Mass Matrix." );

DECLARE_FIELD( StiffnessVectorx,
               "stiffnessVectorx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "x-component of stiffness vector." );

DECLARE_FIELD( StiffnessVectory,
               "stiffnessVectory",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "y-component of stiffness vector." );

DECLARE_FIELD( StiffnessVectorz,
               "stiffnessVectorz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "z-component of stiffness vector." );

DECLARE_FIELD( StiffnessVectorAx,
               "stiffnessVectorAx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "x-component of attenuation stiffness vector." );

DECLARE_FIELD( StiffnessVectorAy,
               "stiffnessVectorAy",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "y-component of attenuation stiffness vector." );

DECLARE_FIELD( StiffnessVectorAz,
               "stiffnessVectorAz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "z-component of attenuation stiffness vector." );

DECLARE_FIELD( DampingVectorx,
               "dampingVectorx",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in x-direction." );

DECLARE_FIELD( DampingVectory,
               "dampingVectory",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in y-direction." );

DECLARE_FIELD( DampingVectorz,
               "dampingVectorz",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal Damping Matrix in z-direction." );

DECLARE_FIELD( ElasticVelocityVp,
               "elasticVelocityVp",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "P-waves speed in the cell" );

DECLARE_FIELD( ElasticVelocityVs,
               "elasticVelocityVs",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "S-waves speed in the cell" );

DECLARE_FIELD( ElasticDensity,
               "elasticDensity",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Medium density of the cell" );

DECLARE_FIELD( ElasticQualityFactorP,
               "elasticQualityFactorP",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Quality factor for P-wave attenuation in the cell" );

DECLARE_FIELD( ElasticQualityFactorS,
               "elasticQualityFactorS",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Quality factor for S-wave attenuation in the cell" );

DECLARE_FIELD( ElasticFreeSurfaceFaceIndicator,
               "elasticFreeSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

DECLARE_FIELD( ElasticFreeSurfaceNodeIndicator,
               "elasticFreeSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a node is on free surface 0 otherwise." );

DECLARE_FIELD( Lambda,
               "lambda",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "First Lame parameter: lambda" );

DECLARE_FIELD( Mu,
               "mu",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Second Lame parameter: mu" );


}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_ELASTICFIELDS */
