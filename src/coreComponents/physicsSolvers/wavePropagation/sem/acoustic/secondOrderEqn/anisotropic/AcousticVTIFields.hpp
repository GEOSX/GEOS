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
 * @file AcousticVTIFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"


namespace geos
{

namespace fields
{

namespace acousticvtifields
{

DECLARE_FIELD( Delta,
               "delta",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Delta thomsen anisotropy parameter" );

DECLARE_FIELD( Epsilon,
               "epsilon",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Epsilon thomsen anisotropy parameter" );

DECLARE_FIELD( F,
               "f",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "f quantity in VTI/TTI Fletcher's equations" );

DECLARE_FIELD( StiffnessVector_p,
               "stiffnessVector_p",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );

DECLARE_FIELD( StiffnessVector_q,
               "stiffnessVector_q",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Stiffness vector contains R_h*Pressure_n." );

DECLARE_FIELD( Pressure_p_nm1,
               "pressure_p_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n-1." );

DECLARE_FIELD( Pressure_p_n,
               "pressure_p_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar pressure at time n." );

DECLARE_FIELD( Pressure_p_np1,
               "pressure_p_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar pressure at time n+1." );

DECLARE_FIELD( Pressure_q_nm1,
               "pressure_q_nm1",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n-1." );

DECLARE_FIELD( Pressure_q_n,
               "pressure_q_n",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n." );

DECLARE_FIELD( Pressure_q_np1,
               "pressure_q_np1",
               array1d< real32 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Scalar auxiliary pressure q at time n+1." );

DECLARE_FIELD( DampingVector_p,
               "dampingVector_p",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for p terms in p equation." );

DECLARE_FIELD( DampingVector_pq,
               "dampingVector_pq",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for q terms in p equation." );

DECLARE_FIELD( DampingVector_q,
               "dampingVector_q",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for q terms in q equation." );

DECLARE_FIELD( DampingVector_qp,
               "dampingVector_qp",
               array1d< real32 >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Diagonal of the Damping Matrix for p terms in q equation." );

DECLARE_FIELD( LateralSurfaceFaceIndicator,
               "lateralSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Free surface indicator, 1 if a face is on a lateral surface 0 otherwise." );

DECLARE_FIELD( LateralSurfaceNodeIndicator,
               "lateralSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Lateral surface indicator, 1 if a face is on a lateral surface 0 otherwise." );

DECLARE_FIELD( BottomSurfaceFaceIndicator,
               "bottomSurfaceFaceIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Bottom surface indicator, 1 if a face is on the bottom surface 0 otherwise." );

DECLARE_FIELD( BottomSurfaceNodeIndicator,
               "bottomSurfaceNodeIndicator",
               array1d< localIndex >,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Bottom surface indicator, 1 if a face is on the bottom surface 0 otherwise." );
}

}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION__HPP_ACOUSTICVTIFIELDS */
