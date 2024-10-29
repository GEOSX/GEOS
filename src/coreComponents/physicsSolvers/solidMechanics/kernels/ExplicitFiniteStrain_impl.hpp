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
 * @file ExplicitFiniteStrain_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_IMPL_HPP_

#include "constitutive/solid/SolidUtilities.hpp"
#include "ExplicitFiniteStrain.hpp"
#include "ExplicitSmallStrain_impl.hpp"
#include "finiteElement/Kinematics.h"

namespace geos
{

namespace solidMechanicsLagrangianFEMKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
ExplicitFiniteStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::ExplicitFiniteStrain( NodeManager & nodeManager,
                                                                                          EdgeManager const & edgeManager,
                                                                                          FaceManager const & faceManager,
                                                                                          localIndex const targetRegionIndex,
                                                                                          SUBREGION_TYPE const & elementSubRegion,
                                                                                          FE_TYPE const & finiteElementSpace,
                                                                                          CONSTITUTIVE_TYPE & inputConstitutiveType,
                                                                                          real64 const dt,
                                                                                          string const elementListName ):
  Base( nodeManager,
        edgeManager,
        faceManager,
        targetRegionIndex,
        elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType,
        dt,
        elementListName )
{}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
inline
void ExplicitFiniteStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::setup( localIndex const k,
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
      stack.uLocal[ a ][ i ] = m_u[ nodeIndex ][ i ];
      stack.varLocal[ a ][ i ] = m_vel[ nodeIndex ][ i ];
    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
inline
void ExplicitFiniteStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                                localIndex const q,
                                                                                                StackVariables & stack ) const
{
  real64 dNdX[ numNodesPerElem ][ 3 ];
  real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

  real64 dUhatdX[3][3] = { {0} };
  real64 dUdX[3][3] = { {0} };
  real64 F[3][3] = { {0} };
  real64 Ldt[3][3] = { {0} };
  real64 fInv[3][3] = { {0} };

  FE_TYPE::gradient( dNdX, stack.varLocal, dUhatdX );
  FE_TYPE::gradient( dNdX, stack.uLocal, dUdX );

  LvArray::tensorOps::scale< 3, 3 >( dUhatdX, m_dt );

  // calculate du/dX
  LvArray::tensorOps::scaledCopy< 3, 3 >( F, dUhatdX, 0.5 );
  LvArray::tensorOps::add< 3, 3 >( F, dUdX );
  LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
  LvArray::tensorOps::invert< 3 >( fInv, F );

  // chain rule: calculate dv/dx^(n+1/2) = dv/dX * dX/dx^(n+1/2)
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( Ldt, dUhatdX, fInv );

  // calculate gradient (end of step)
  LvArray::tensorOps::copy< 3, 3 >( F, dUhatdX );
  LvArray::tensorOps::add< 3, 3 >( F, dUdX );
  LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
  real64 const detF = LvArray::tensorOps::invert< 3 >( fInv, F );

  real64 Rot[ 3 ][ 3 ]{};
  real64 Dadt[ 6 ]{};
  HughesWinget( Rot, Dadt, Ldt );

  real64 stress[ 6 ]{};
  constitutive::SolidUtilities::
    hypoUpdate_StressOnly( m_constitutiveUpdate,
                           k, q, m_dt, Dadt, Rot, stress );

  real64 P[ 3 ][ 3 ]{};
  LvArray::tensorOps::Rij_eq_symAikBjk< 3 >( P, stress, fInv );
  LvArray::tensorOps::scale< 3, 3 >( P, -detJ * detF );

  FE_TYPE::plusGradNajAij( dNdX, P, stack.fLocal );
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
ExplicitFiniteStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::kernelLaunch( localIndex const numElems,
                                                                                  KERNEL_TYPE const & kernelComponent )
{
  return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
}

#undef UPDATE_STRESS


} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_IMPL_HPP_
