/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsDamage_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_IMPL_HPP_

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsDamage.hpp"

namespace geos
{

namespace poromechanicsDamageKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
SinglePhasePoromechanicsDamage< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
SinglePhasePoromechanicsDamage( NodeManager const & nodeManager,
                                EdgeManager const & edgeManager,
                                FaceManager const & faceManager,
                                localIndex const targetRegionIndex,
                                SUBREGION_TYPE const & elementSubRegion,
                                FE_TYPE const & finiteElementSpace,
                                CONSTITUTIVE_TYPE & inputConstitutiveType,
                                arrayView1d< globalIndex const > const inputDispDofNumber,
                                globalIndex const rankOffset,
                                CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                arrayView1d< real64 > const inputRhs,
                                real64 const inputDt,
                                real64 const (&gravityVector)[3],
                                string const inputFlowDofKey,
                                integer const performStressInitialization,
                                string const fluidModelKey ):
  Base( nodeManager,
        edgeManager,
        faceManager,
        targetRegionIndex,
        elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType,
        inputDispDofNumber,
        rankOffset,
        inputMatrix,
        inputRhs,
        inputDt,
        gravityVector,
        inputFlowDofKey,
        performStressInitialization,
        fluidModelKey ),
  m_fluidPressureGradient( elementSubRegion.template getReference< array2d< real64 > >( "pressureGradient" ) )
{}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanicsDamage< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
smallStrainUpdate( localIndex const k,
                   localIndex const q,
                   StackVariables & stack ) const
{
  real64 porosity = 0.0;
  real64 porosity_n = 0.0;
  real64 dPorosity_dVolStrain = 0.0;
  real64 dPorosity_dPressure = 0.0;
  real64 dPorosity_dTemperature = 0.0;
  real64 dSolidDensity_dPressure = 0.0;

  real64 fluidPressureGradient[3]{};

  for( integer i=0; i<3; ++i )
  {
    fluidPressureGradient[i] = m_fluidPressureGradient( k, i );
  }

  // Step 1: call the constitutive model to evaluate the total stress and compute porosity
  m_constitutiveUpdate.smallStrainUpdatePoromechanics( k, q,
                                                       m_dt,
                                                       m_pressure_n[k],
                                                       m_pressure[k],
                                                       stack.temperature,
                                                       stack.deltaTemperatureFromLastStep,
                                                       stack.strainIncrement,
                                                       stack.totalStress,
                                                       stack.dTotalStress_dPressure,
                                                       stack.dTotalStress_dTemperature,
                                                       stack.stiffness,
                                                       m_performStressInitialization,
                                                       porosity,
                                                       porosity_n,
                                                       dPorosity_dVolStrain,
                                                       dPorosity_dPressure,
                                                       dPorosity_dTemperature,
                                                       dSolidDensity_dPressure );

  // Step 2: compute fracture flow term and its derivative w.r.t pressure
  m_constitutiveUpdate.computeFractureFlowTerm( k, q,
                                                fluidPressureGradient,
                                                stack.fractureFlowTerm,
                                                stack.dFractureFlowTerm_dPressure );

  // Step 3: compute the body force
  Base::computeBodyForce( k, q,
                          porosity,
                          dPorosity_dVolStrain,
                          dPorosity_dPressure,
                          dPorosity_dTemperature,
                          dSolidDensity_dPressure,
                          stack );

  // Step 3: compute fluid mass increment
  Base::computeFluidIncrement( k, q,
                               porosity,
                               porosity_n,
                               dPorosity_dVolStrain,
                               dPorosity_dPressure,
                               dPorosity_dTemperature,
                               stack );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanicsDamage< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
assembleMomentumBalanceTerms( real64 const ( &N )[numNodesPerElem],
                              real64 const ( &dNdX )[numNodesPerElem][3],
                              real64 const & detJxW,
                              StackVariables & stack ) const
{
  using namespace PDEUtilities;

  constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
  constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;
  constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;

  // Step 1: compute local linear momentum balance residual
  LinearFormUtilities::compute< displacementTestSpace,
                                DifferentialOperator::SymmetricGradient >
  (
    stack.localResidualMomentum,
    dNdX,
    stack.totalStress,
    -detJxW );

  if( m_gravityAcceleration > 0.0 )
  {
    LinearFormUtilities::compute< displacementTestSpace,
                                  DifferentialOperator::Identity >
    (
      stack.localResidualMomentum,
      N,
      stack.bodyForce,
      detJxW );
  }

  LinearFormUtilities::compute< displacementTestSpace,
                                DifferentialOperator::Identity >
  (
    stack.localResidualMomentum,
    N,
    stack.fractureFlowTerm,
    detJxW );

  // Step 2: compute local linear momentum balance residual derivatives with respect to displacement
  BilinearFormUtilities::compute< displacementTestSpace,
                                  displacementTrialSpace,
                                  DifferentialOperator::SymmetricGradient,
                                  DifferentialOperator::SymmetricGradient >
  (
    stack.dLocalResidualMomentum_dDisplacement,
    dNdX,
    stack.stiffness, // fourth-order tensor handled via DiscretizationOps
    dNdX,
    -detJxW );

  if( m_gravityAcceleration > 0.0 )
  {
    BilinearFormUtilities::compute< displacementTestSpace,
                                    displacementTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Divergence >
    (
      stack.dLocalResidualMomentum_dDisplacement,
      N,
      stack.dBodyForce_dVolStrainIncrement,
      dNdX,
      detJxW );
  }

  // Step 3: compute local linear momentum balance residual derivatives with respect to pressure
  BilinearFormUtilities::compute< displacementTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::SymmetricGradient,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMomentum_dPressure,
    dNdX,
    stack.dTotalStress_dPressure,
    1.0,
    -detJxW );

  if( m_gravityAcceleration > 0.0 )
  {
    BilinearFormUtilities::compute< displacementTestSpace,
                                    pressureTrialSpace,
                                    DifferentialOperator::Identity,
                                    DifferentialOperator::Identity >
    (
      stack.dLocalResidualMomentum_dPressure,
      N,
      stack.dBodyForce_dPressure,
      1.0,
      detJxW );
  }

  BilinearFormUtilities::compute< displacementTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMomentum_dPressure,
    N,
    stack.dFractureFlowTerm_dPressure,
    1.0,
    detJxW );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanicsDamage< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack ) const
{
  // Step 1: compute displacement finite element basis functions (N), basis function derivatives (dNdX), and
  // determinant of the Jacobian transformation matrix times the quadrature weight (detJxW)
  real64 N[numNodesPerElem]{};
  real64 dNdX[numNodesPerElem][3]{};
  FE_TYPE::calcN( q, stack.feStack, N );
  real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal,
                                                                           stack.feStack, dNdX );

  // Step 2: compute strain increment
  LvArray::tensorOps::fill< 6 >( stack.strainIncrement, 0.0 );
  FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, stack.strainIncrement );

  // Step 3: compute 1) the total stress, 2) the body force terms, and 3) the fluidMassIncrement
  // using quantities returned by the PorousSolid constitutive model.
  // This function also computes the derivatives of these three quantities wrt primary variables
  smallStrainUpdate( k, q, stack );

  // Step 4: use the total stress and the body force to increment the local momentum balance residual
  // This function also fills the local Jacobian rows corresponding to the momentum balance.
  assembleMomentumBalanceTerms( N, dNdX, detJxW, stack );

  // Step 5: use the fluid mass increment to increment the local mass balance residual
  // This function also fills the local Jacobian rows corresponding to the mass balance.
  Base::assembleElementBasedFlowTerms( dNdX, detJxW, stack );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 SinglePhasePoromechanicsDamage< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  real64 const maxForce = Base::complete( k, stack );

  return maxForce;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64 SinglePhasePoromechanicsDamage< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
kernelLaunch( localIndex const numElems,
              KERNEL_TYPE const & kernelComponent )
{
  GEOS_MARK_FUNCTION;

  // Define a RAJA reduction variable to get the maximum residual contribution.
  RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

  forAll< POLICY >( numElems,
                    [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    typename KERNEL_TYPE::StackVariables stack;

    kernelComponent.setup( k, stack );
    for( integer q=0; q<KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
    {
      kernelComponent.quadraturePointKernel( k, q, stack );
    }
    maxResidual.max( kernelComponent.complete( k, stack ) );
  } );
  return maxResidual.get();
}

} // namespace poromechanicsDamageKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_IMPL_HPP_
