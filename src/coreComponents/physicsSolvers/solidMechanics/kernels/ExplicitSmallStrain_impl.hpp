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
 * @file ExplictSmallStrain_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITSMALLTRAIN_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITSMALLTRAIN_IMPL_HPP_

//#include "ExplicitFiniteStrain.hpp"
#include "ExplicitSmallStrain.hpp"


namespace geos
{

/// Namespace to contain the solid mechanics kernels.
namespace solidMechanicsLagrangianFEMKernels
{


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
ExplicitSmallStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::ExplicitSmallStrain( NodeManager & nodeManager,
                                                                                        EdgeManager const & edgeManager,
                                                                                        FaceManager const & faceManager,
                                                                                        localIndex const targetRegionIndex,
                                                                                        SUBREGION_TYPE const & elementSubRegion,
                                                                                        FE_TYPE const & finiteElementSpace,
                                                                                        CONSTITUTIVE_TYPE & inputConstitutiveType,
                                                                                        real64 const dt,
                                                                                        string const elementListName ):
  Base( elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType ),
  m_X( nodeManager.referencePosition()),
  m_u( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
  m_vel( nodeManager.getField< fields::solidMechanics::velocity >() ),
  m_acc( nodeManager.getField< fields::solidMechanics::acceleration >() ),
  m_dt( dt ),
  m_elementList( elementSubRegion.template getReference< SortedArray< localIndex > >( elementListName ).toViewConst() )
{
  GEOS_UNUSED_VAR( edgeManager );
  GEOS_UNUSED_VAR( faceManager );
  GEOS_UNUSED_VAR( targetRegionIndex );
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
inline
void ExplicitSmallStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::setup( localIndex const k,
                                                                               StackVariables & stack ) const
{
  for( localIndex a=0; a< numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    for( int i=0; i<numDofPerTrialSupportPoint; ++i )
    {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
#endif

#if UPDATE_STRESS==2
      stack.varLocal[ a ][ i ] = m_vel[ nodeIndex ][ i ] * m_dt;
#else
      stack.varLocal[ a ][ i ] = m_u[ nodeIndex ][ i ];
#endif
    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
inline
void ExplicitSmallStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                               localIndex const q,
                                                                                               StackVariables & stack ) const
{
//#define USE_JACOBIAN
#if !defined( USE_JACOBIAN )
  real64 dNdX[ numNodesPerElem ][ 3 ];
  real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );
  /// Macro to substitute in the shape function derivatives.
  real64 strain[6] = {0};
  //real64 timeIncrement = 0.0;
  FE_TYPE::symmetricGradient( dNdX, stack.varLocal, strain );

  real64 stressLocal[ 6 ] = {0};
#if UPDATE_STRESS == 2
  m_constitutiveUpdate.smallStrainUpdate_StressOnly( k, q, m_dt, strain, stressLocal );
#else
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );
#endif

  for( localIndex c = 0; c < 6; ++c )
  {
#if UPDATE_STRESS == 2
    stressLocal[ c ] *= -detJ;
#elif UPDATE_STRESS == 1
    stressLocal[ c ] = -( stressLocal[ c ] + m_constitutiveUpdate.m_newStress( k, q, c ) ) * detJ;   // TODO: decide on
                                                                                                     // initial stress
                                                                                                     // strategy
#else
    stressLocal[ c ] *= -detJ;
#endif
  }

  FE_TYPE::plusGradNajAij( dNdX, stressLocal, stack.fLocal );

#else
  real64 invJ[3][3];
  real64 const detJ = FE_TYPE::inverseJacobianTransformation( q, stack.xLocal, invJ );

  real64 strain[6] = {0};
  FE_TYPE::symmetricGradient( q, invJ, stack.varLocal, strain );

  real64 stressLocal[ 6 ] = {0};
#if UPDATE_STRESS == 2
  m_constitutiveUpdate.smallStrainUpdate_StressOnly( k, q, m_dt, strain, stressLocal );
#else
  m_constitutiveUpdate.smallStrainNoStateUpdate_StressOnly( k, q, strain, stressLocal );
#endif

  for( localIndex c = 0; c < 6; ++c )
  {
#if UPDATE_STRESS == 2
    stressLocal[ c ] *= detJ;
#elif UPDATE_STRESS == 1
    stressLocal[ c ] = ( stressLocal[ c ] + m_constitutiveUpdate.m_newStress( k, q, c ) ) * DETJ;   // TODO: decide on
                                                                                                    // initial stress
                                                                                                    // strategy
#else
    stressLocal[ c ] *= DETJ;
#endif
  }

  FE_TYPE::plusGradNajAij( q, invJ, stressLocal, stack.fLocal );
#endif
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
inline
real64 ExplicitSmallStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::complete( localIndex const k,
                                                                                    StackVariables const & stack ) const
{
  for( localIndex a = 0; a < numNodesPerElem; ++a )
  {
    localIndex const nodeIndex = m_elemsToNodes( k, a );
    for( int b = 0; b < numDofPerTestSupportPoint; ++b )
    {
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_acc( nodeIndex, b ), stack.fLocal[ a ][ b ] );
    }
  }
  return 0;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
ExplicitSmallStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::kernelLaunch( localIndex const numElems,
                                                                                 KERNEL_TYPE const & kernelComponent )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( numElems );

  localIndex const numProcElems = kernelComponent.m_elementList.size();
  forAll< POLICY >( numProcElems,
                    [=] GEOS_DEVICE ( localIndex const index )
  {
    localIndex const k = kernelComponent.m_elementList[ index ];

    typename KERNEL_TYPE::StackVariables stack;

    kernelComponent.setup( k, stack );
    for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
    {
      kernelComponent.quadraturePointKernel( k, q, stack );
    }
    kernelComponent.complete( k, stack );
  } );
  return 0;
}


#undef UPDATE_STRESS

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITSMALLTRAIN_IMPL_HPP_
