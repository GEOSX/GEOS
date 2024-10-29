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
 * @file SolidMechanicsLagrangianFEMKernels.hpp
 */

#pragma once

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace solidMechanicsLagrangianFEMKernels
{

inline void velocityUpdate( arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
                            arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
                            real64 const dt )
{
  GEOS_MARK_FUNCTION;

  localIndex const N = acceleration.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOS_DEVICE ( localIndex const i )
  {
    LvArray::tensorOps::scaledAdd< 3 >( velocity[ i ], acceleration[ i ], dt );
    LvArray::tensorOps::fill< 3 >( acceleration[ i ], 0 );
  } );
}

inline void velocityUpdate( arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
                            arrayView1d< real64 const > const & mass,
                            arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
                            real64 const dt,
                            SortedArrayView< localIndex const > const & indices )
{
  GEOS_MARK_FUNCTION;

  forAll< parallelDevicePolicy<> >( indices.size(), [=] GEOS_DEVICE ( localIndex const i )
  {
    localIndex const a = indices[ i ];
    LvArray::tensorOps::scale< 3 >( acceleration[ a ], 1.0 / mass[ a ] );
    LvArray::tensorOps::scaledAdd< 3 >( velocity[ a ], acceleration[ a ], dt );
  } );
}

inline void displacementUpdate( arrayView2d< real64 const, nodes::VELOCITY_USD > const & velocity,
                                arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat,
                                arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u,
                                real64 const dt )
{
  GEOS_MARK_FUNCTION;

  localIndex const N = velocity.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOS_DEVICE ( localIndex const i )
  {
    LvArray::tensorOps::scaledCopy< 3 >( uhat[ i ], velocity[ i ], dt );
    LvArray::tensorOps::add< 3 >( u[ i ], uhat[ i ] );
  } );
}


/**
 * @struct Structure to wrap templated function that implements the explicit time integration kernel.
 */
struct ExplicitKernel
{

  static inline real64
  calculateSingleNodalForce( localIndex const k,
                             localIndex const targetNode,
                             localIndex const numQuadraturePoints,
                             arrayView4d< real64 const > const & dNdX,
                             arrayView2d< real64 const > const & detJ,
                             arrayView3d< real64 const, solid::STRESS_USD > const & stress,
                             real64 ( & force )[ 3 ] )
  {
    GEOS_MARK_FUNCTION;
    localIndex const & a = targetNode;

    //Compute Quadrature
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      force[ 0 ] -= ( stress( k, q, 0 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 5 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 4 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );
      force[ 1 ] -= ( stress( k, q, 5 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 1 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 3 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );
      force[ 2 ] -= ( stress( k, q, 4 ) * dNdX( k, q, a, 0 ) +
                      stress( k, q, 3 ) * dNdX( k, q, a, 1 ) +
                      stress( k, q, 2 ) * dNdX( k, q, a, 2 ) ) * detJ( k, q );

    }//quadrature loop

    return 0;
  }

};


} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos
