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
 * @file ImplicitSmallStrainQuasiStatic_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_IMPL_HPP_

#include "ImplicitSmallStrainQuasiStatic.hpp"

namespace geos
{

namespace solidMechanicsLagrangianFEMKernels
{


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
ImplicitSmallStrainQuasiStatic( NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager const & faceManager,
                                localIndex const targetRegionIndex,
                                SUBREGION_TYPE const & elementSubRegion,
                                FE_TYPE const & finiteElementSpace,
                                CONSTITUTIVE_TYPE & inputConstitutiveType,
                                arrayView1d< globalIndex const > const inputDofNumber,
                                globalIndex const rankOffset,
                                CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                arrayView1d< real64 > const inputRhs,
                                real64 const inputDt,
                                real64 const (&inputGravityVector)[3] ):
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
        inputDt ),
  m_X( nodeManager.referencePosition()),
  m_disp( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
  m_uhat( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
  m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
  m_density( inputConstitutiveType.getDensity() )
{}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );

  localIndex const numSupportPoints = m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );

  stack.numRows =  3 * numSupportPoints;
  stack.numCols = stack.numRows;

  // #pragma unroll
  for( localIndex a = 0; a < numSupportPoints; ++a )
  {
    localIndex const localNodeIndex = m_elemsToNodes( k, a );

    // #pragma unroll
    for( int i = 0; i < numDofPerTestSupportPoint; ++i )
    {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
#endif
      stack.u_local[ a ][i] = m_disp[ localNodeIndex ][i];
      stack.uhat_local[ a ][i] = m_uhat[ localNodeIndex ][i];
      stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
    }
  }
  // Add stabilization to block diagonal parts of the local jacobian
  // (this is a no-operation with FEM classes)
  real64 const stabilizationScaling = computeStabilizationScaling( k );
  m_finiteElementSpace.template addGradGradStabilizationMatrix
  < FE_TYPE, numDofPerTrialSupportPoint, true >( stack.feStack,
                                                 stack.localJacobian,
                                                 -stabilizationScaling );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename STRESS_MODIFIER >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::quadraturePointKernel( localIndex const k,
                                                                                                          localIndex const q,
                                                                                                          StackVariables & stack,
                                                                                                          STRESS_MODIFIER && stressModifier ) const
{
  real64 dNdX[ numNodesPerElem ][ 3 ];
  real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, stack.feStack, dNdX );

  real64 strainInc[6] = {0};
  real64 stress[6] = {0};

  typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;

  FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainInc );

  m_constitutiveUpdate.smallStrainUpdate( k, q, m_dt, strainInc, stress, stiffness );

  stressModifier( stress );
  // #pragma unroll
  for( localIndex i=0; i<6; ++i )
  {
    stress[i] *= -detJxW;
  }

  real64 const gravityForce[3] = { m_gravityVector[0] * m_density( k, q )* detJxW,
                                   m_gravityVector[1] * m_density( k, q )* detJxW,
                                   m_gravityVector[2] * m_density( k, q )* detJxW };

  real64 N[numNodesPerElem];
  FE_TYPE::calcN( q, stack.feStack, N );
  FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                   stress,
                                   N,
                                   gravityForce,
                                   reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );
  real64 const stabilizationScaling = computeStabilizationScaling( k );
  m_finiteElementSpace.template addEvaluatedGradGradStabilizationVector< FE_TYPE, numDofPerTrialSupportPoint >( stack.feStack,
                                                                                                                stack.uhat_local,
                                                                                                                reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual),
                                                                                                                -stabilizationScaling );
#if !defined( GEOS_USE_HIP )
  stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJxW, stack.localJacobian );
# else
  stiffness.template BTDB< numNodesPerElem >( dNdX, -detJxW, stack.localJacobian ); // need to use full BTDB compute for hip
#endif
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::complete( localIndex const k,
                                                                                               StackVariables & stack ) const
{
  GEOS_UNUSED_VAR( k );
  real64 maxForce = 0;

#if !defined( GEOS_USE_HIP )
  // TODO: Does this work if BTDB is non-symmetric?
  CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );
#endif
  localIndex const numSupportPoints = m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );

  // #pragma unroll
  for( int localNode = 0; localNode < numSupportPoints; ++localNode )
  {
    // #pragma unroll
    for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ numDofPerTestSupportPoint * localNode + dim ] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() )
        continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localRowDofIndex,
                                                                              stack.localJacobian[ numDofPerTestSupportPoint * localNode + dim ],
                                                                              stack.numRows );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[ dof ], stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] );
      maxForce = fmax( maxForce, fabs( stack.localResidual[ numDofPerTestSupportPoint * localNode + dim ] ) );
    }
  }


  return maxForce;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
GEOS_FORCE_INLINE
real64
ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::kernelLaunch( localIndex const numElems,
                                                                                            KERNEL_TYPE const & kernelComponent )
{
  return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
}


} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_IMPL_HPP_
