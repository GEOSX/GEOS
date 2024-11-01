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
 * @file FixedStressThermoPoromechanics_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_FIXEDSTRESSTHERMOPOROMECHANICS_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_FIXEDSTRESSTHERMOPOROMECHANICS_IMPL_HPP_

#include "FixedStressThermoPoromechanics.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{

namespace solidMechanicsLagrangianFEMKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >

FixedStressThermoPoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
FixedStressThermoPoromechanics( NodeManager const & nodeManager,
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
  m_bulkDensity( elementSubRegion.template getField< fields::poromechanics::bulkDensity >() ),
  m_pressure( elementSubRegion.template getField< fields::flow::pressure >() ),
  m_pressure_n( elementSubRegion.template getField< fields::flow::pressure_n >() ),
  m_initialTemperature( elementSubRegion.template getField< fields::flow::initialTemperature >() ),
  m_temperature( elementSubRegion.template getField< fields::flow::temperature >() ),
  m_temperature_n( elementSubRegion.template getField< fields::flow::temperature_n >() )
{}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FixedStressThermoPoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );
  localIndex const numSupportPoints =
    m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
  stack.numRows = 3 * numSupportPoints;
  stack.numCols = stack.numRows;

  for( localIndex a = 0; a < numSupportPoints; ++a )
  {
    localIndex const localNodeIndex = m_elemsToNodes( k, a );

    for( int i = 0; i < 3; ++i )
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
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FixedStressThermoPoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack ) const
{
  real64 dNdX[ numNodesPerElem ][ 3 ];
  real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal,
                                                                           stack.feStack, dNdX );

  real64 strainInc[6] = {0};
  real64 totalStress[6] = {0};

  typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;

  FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainInc );

  // Evaluate total stress and its derivatives
  // TODO: allow for a customization of the kernel to pass the average pressure to the small strain update (to account for cap pressure
  // later)
  m_constitutiveUpdate.smallStrainUpdatePoromechanicsFixedStress( k, q,
                                                                  m_dt,
                                                                  m_pressure[k],
                                                                  m_pressure_n[k],
                                                                  m_temperature[k],
                                                                  m_temperature_n[k],
                                                                  strainInc,
                                                                  totalStress,
                                                                  stiffness );

  for( localIndex i=0; i<6; ++i )
  {
    totalStress[i] *= -detJxW;
  }

  // Here we consider the bodyForce is purely from the solid
  // Warning: here, we lag (in iteration) the displacement dependence of bulkDensity
  real64 const gravityForce[3] = { m_gravityVector[0] * m_bulkDensity( k, q )* detJxW,
                                   m_gravityVector[1] * m_bulkDensity( k, q )* detJxW,
                                   m_gravityVector[2] * m_bulkDensity( k, q )* detJxW };

  real64 N[numNodesPerElem];
  FE_TYPE::calcN( q, stack.feStack, N );
  FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                   totalStress,
                                   N,
                                   gravityForce,
                                   reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );
  real64 const stabilizationScaling = computeStabilizationScaling( k );
  m_finiteElementSpace.template
  addEvaluatedGradGradStabilizationVector< FE_TYPE,
                                           numDofPerTrialSupportPoint >( stack.feStack,
                                                                         stack.uhat_local,
                                                                         reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual),
                                                                         -stabilizationScaling );
  stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJxW, stack.localJacobian );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 FixedStressThermoPoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  GEOS_UNUSED_VAR( k );
  real64 maxForce = 0;

  // TODO: Does this work if BTDB is non-symmetric?
  CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );
  localIndex const numSupportPoints =
    m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
  for( int localNode = 0; localNode < numSupportPoints; ++localNode )
  {
    for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
    {
      localIndex const dof =
        LvArray::integerConversion< localIndex >( stack.localRowDofIndex[ numDofPerTestSupportPoint * localNode + dim ] - m_dofRankOffset );
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
real64
FixedStressThermoPoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::kernelLaunch( localIndex const numElems,
                                                                                            KERNEL_TYPE const & kernelComponent )
{
  return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
}

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_FIXEDSTRESSTHERMOPOROMECHANICS_IMPL_HPP_
