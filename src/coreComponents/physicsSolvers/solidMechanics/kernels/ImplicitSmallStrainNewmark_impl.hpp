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
 * @file ImplicitSmallStrainNewmark_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINNEWMARK_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINNEWMARK_IMPL_HPP_

#include "ImplicitSmallStrainNewmark.hpp"
#include "ImplicitSmallStrainQuasiStatic_impl.hpp"

namespace geos
{

namespace solidMechanicsLagrangianFEMKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
ImplicitSmallStrainNewmark< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
ImplicitSmallStrainNewmark( NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            localIndex const targetRegionIndex,
                            SUBREGION_TYPE const & elementSubRegion,
                            FE_TYPE const & finiteElementSpace,
                            CONSTITUTIVE_TYPE & inputConstitutiveType,
                            arrayView1d< globalIndex const > const & inputDofNumber,
                            globalIndex const rankOffset,
                            CRSMatrixView< real64, globalIndex const > const inputMatrix,
                            arrayView1d< real64 > const inputRhs,
                            real64 const inputDt,
                            real64 const (&inputGravityVector)[3],
                            real64 const inputNewmarkGamma,
                            real64 const inputNewmarkBeta,
                            real64 const inputMassDamping,
                            real64 const inputStiffnessDamping ):
  // real64 const inputDt ):
  Base( nodeManager,
        edgeManager,
        faceManager,
        targetRegionIndex,
        elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType,
        inputDofNumber,
        rankOffset,
        inputMatrix,
        inputRhs,
        inputDt,
        inputGravityVector ),
  m_vtilde( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
  m_uhattilde( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
  m_newmarkGamma( inputNewmarkGamma ),
  m_newmarkBeta( inputNewmarkBeta ),
  m_massDamping( inputMassDamping ),
  m_stiffnessDamping( inputStiffnessDamping )
  //m_dt( inputDt )
{}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
inline
void ImplicitSmallStrainNewmark< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  for( localIndex a=0; a<numNodesPerElem; ++a )
  {
    localIndex const localNodeIndex = m_elemsToNodes( k, a );
    for( localIndex i=0; i<numDofPerTrialSupportPoint; ++i )
    {
      stack.vtilde_local[ a ][ i ] = m_vtilde[ localNodeIndex ][ i ];
      stack.uhattilde_local[ a ][ i ] = m_uhattilde[ localNodeIndex ][ i ];
    }
  }
  Base::setup( k, stack );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
inline
void ImplicitSmallStrainNewmark< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack ) const
{

  Base::quadraturePointKernel( k, q, stack );
  real64 detJ=0;

  real64 N[numNodesPerElem];
  FE_TYPE::calcN( q, N );

  for( int a=0; a<numNodesPerElem; ++a )
  {
    for( int b=a; b<numNodesPerElem; ++b )
    {
      real64 const integrationFactor = m_density( k, q ) * N[a] * N[b] * detJ;
      real64 const temp1 = ( m_massDamping * m_newmarkGamma/( m_newmarkBeta * m_dt )
                             + 1.0 / ( m_newmarkBeta * m_dt * m_dt ) )* integrationFactor;

      constexpr int nsdof = numDofPerTestSupportPoint;
      for( int i=0; i<nsdof; ++i )
      {
        real64 const acc = 1.0 / ( m_newmarkBeta * m_dt * m_dt ) * ( stack.uhat_local[b][i] - stack.uhattilde_local[b][i] );
        real64 const vel = stack.vtilde_local[b][i] +
                           m_newmarkGamma/( m_newmarkBeta * m_dt ) *( stack.uhat_local[b][i]
                                                                      - stack.uhattilde_local[b][i] );

        stack.dRdU_InertiaMassDamping[ a*nsdof+i][ b*nsdof+i ] -= temp1;
        stack.localResidual[ a*nsdof+i ] -= ( m_massDamping * vel + acc ) * integrationFactor;
      }
    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
inline
real64 ImplicitSmallStrainNewmark< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{

  for( int a=0; a<numNodesPerElem; ++a )
  {
    for( int b=0; b<numNodesPerElem; ++b )
    {
      for( int i=0; i<numDofPerTestSupportPoint; ++i )
      {
        for( int j=0; j<numDofPerTrialSupportPoint; ++j )
        {
          stack.localResidual[ a*numDofPerTestSupportPoint+i ] =
            stack.localResidual[ a*numDofPerTestSupportPoint+i ] +
            m_stiffnessDamping * stack.localJacobian[ a*numDofPerTestSupportPoint+i][ b*numDofPerTrialSupportPoint+j ] *
            ( stack.vtilde_local[b][j] + m_newmarkGamma/(m_newmarkBeta * m_dt)*(stack.uhat_local[b][j]-stack.uhattilde_local[b][j]) );

          stack.localJacobian[a*numDofPerTestSupportPoint+i][b*numDofPerTrialSupportPoint+j] =
            stack.localJacobian[a*numDofPerTestSupportPoint+i][b*numDofPerTrialSupportPoint+j] +
            stack.localJacobian[a][b] * (1.0 + m_stiffnessDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) ) +
            stack.dRdU_InertiaMassDamping[ a*numDofPerTestSupportPoint+i ][ b*numDofPerTrialSupportPoint+j ];
        }
      }
    }
  }

  for( int a=0; a<stack.maxNumRows; ++a )
  {
    for( int b=0; b<stack.maxNumCols; ++b )
    {
      stack.localJacobian[a][b] += stack.localJacobian[a][b] * (1.0 + m_stiffnessDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) )
                                   + stack.dRdU_InertiaMassDamping[ a ][ b ];
    }
  }

  return Base::complete( k, stack );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
ImplicitSmallStrainNewmark< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::kernelLaunch( localIndex const numElems,
                                                                                        KERNEL_TYPE const & kernelComponent )
{
  return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
}


} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINNEWMARK_IMPL_HPP_
