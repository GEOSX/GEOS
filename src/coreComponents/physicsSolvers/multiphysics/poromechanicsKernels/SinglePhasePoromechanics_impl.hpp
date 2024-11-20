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
 * @file SinglePhasePoromechanics_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_IMPL_HPP_

#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"

namespace geos
{

namespace poromechanicsKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
SinglePhasePoromechanics( NodeManager const & nodeManager,
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
        fluidModelKey ),
  m_fluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density() ),
  m_fluidDensity_n( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density_n() ),
  m_dFluidDensity_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).dDensity_dPressure() ),
  m_performStressInitialization( performStressInitialization )
{}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

  // Step 1: call the constitutive model to evaluate the total stress and compute porosity
  m_constitutiveUpdate.smallStrainUpdatePoromechanics( k, q,
                                                       m_dt,
                                                       m_pressure[k],
                                                       m_pressure_n[k],
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

  // Step 2: compute the body force
  computeBodyForce( k, q,
                    porosity,
                    dPorosity_dVolStrain,
                    dPorosity_dPressure,
                    dPorosity_dTemperature,
                    dSolidDensity_dPressure,
                    stack );

  // Step 3: compute fluid mass increment
  computeFluidIncrement( k, q,
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
void SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computeBodyForce( localIndex const k,
                  localIndex const q,
                  real64 const & porosity,
                  real64 const & dPorosity_dVolStrain,
                  real64 const & dPorosity_dPressure,
                  real64 const & dPorosity_dTemperature,
                  real64 const & dSolidDensity_dPressure,
                  StackVariables & stack ) const
{
  GEOS_UNUSED_VAR( dPorosity_dTemperature );

  real64 const mixtureDensity = ( 1.0 - porosity ) * m_solidDensity( k, q ) + porosity * m_fluidDensity( k, q );
  real64 const dMixtureDens_dVolStrainIncrement = dPorosity_dVolStrain * ( -m_solidDensity( k, q ) + m_fluidDensity( k, q ) );
  real64 const dMixtureDens_dPressure = dPorosity_dPressure * ( -m_solidDensity( k, q ) + m_fluidDensity( k, q ) )
                                        + ( 1.0 - porosity ) * dSolidDensity_dPressure
                                        + porosity * m_dFluidDensity_dPressure( k, q );

  LvArray::tensorOps::scaledCopy< 3 >( stack.bodyForce, m_gravityVector, mixtureDensity );
  LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dVolStrainIncrement, m_gravityVector, dMixtureDens_dVolStrainIncrement );
  LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dPressure, m_gravityVector, dMixtureDens_dPressure );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computeFluidIncrement( localIndex const k,
                       localIndex const q,
                       real64 const & porosity,
                       real64 const & porosity_n,
                       real64 const & dPorosity_dVolStrain,
                       real64 const & dPorosity_dPressure,
                       real64 const & dPorosity_dTemperature,
                       StackVariables & stack ) const
{
  GEOS_UNUSED_VAR( dPorosity_dTemperature );

  stack.fluidMassIncrement = porosity * m_fluidDensity( k, q ) - porosity_n * m_fluidDensity_n( k, q );
  stack.dFluidMassIncrement_dVolStrainIncrement = dPorosity_dVolStrain * m_fluidDensity( k, q );
  stack.dFluidMassIncrement_dPressure = dPorosity_dPressure * m_fluidDensity( k, q ) + porosity * m_dFluidDensity_dPressure( k, q );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

  LinearFormUtilities::compute< displacementTestSpace,
                                DifferentialOperator::Identity >
  (
    stack.localResidualMomentum,
    N,
    stack.bodyForce,
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

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
assembleElementBasedFlowTerms( real64 const ( &dNdX )[numNodesPerElem][3],
                               real64 const & detJxW,
                               StackVariables & stack ) const
{
  using namespace PDEUtilities;

  constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
  constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;
  constexpr FunctionSpace pressureTestSpace = pressureTrialSpace;

  // Step 1: compute local mass balance residual
  LinearFormUtilities::compute< pressureTestSpace,
                                DifferentialOperator::Identity >
  (
    stack.localResidualMass,
    1.0,
    stack.fluidMassIncrement,
    detJxW );

  // Step 2: compute local mass balance residual derivatives with respect to displacement
  BilinearFormUtilities::compute< pressureTestSpace,
                                  displacementTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Divergence >
  (
    stack.dLocalResidualMass_dDisplacement,
    1.0,
    stack.dFluidMassIncrement_dVolStrainIncrement,
    dNdX,
    detJxW );

  // Step 3: compute local mass balance residual derivatives with respect to pressure
  BilinearFormUtilities::compute< pressureTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMass_dPressure,
    1.0,
    stack.dFluidMassIncrement_dPressure,
    1.0,
    detJxW );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack ) const
{
  // Governing equations (strong form)
  // ---------------------------------
  //
  //   divergence( totalStress ) + bodyForce = 0                (quasi-static linear momentum balance)
  //   dFluidMass_dTime + divergence( fluidMassFlux ) = source  (fluid phase mass balance)
  //
  // with currently the following dependencies on the strainIncrement tensor and pressure
  //
  //   totalStress      = totalStress( strainIncrement, pressure)
  //   bodyForce        = bodyForce( strainIncrement, pressure)
  //   fluidMass        = fluidMass( strainIncrement, pressure)
  //   fluidMassFlux    = fludiMassFlux( pressure)
  //
  // Note that the fluidMassFlux will depend on the strainIncrement if a stress-dependent constitutive
  // model is assumed. A dependency on pressure can also occur in the source term of the mass
  // balance equation, e.g. if a Peaceman well model is used.
  //
  // In this kernel cell-based contributions to Jacobian matrix and residual
  // vector are computed. The face-based contributions, namely associated with the fluid
  // mass flux term are computed in a different kernel. The source term in the mass balance
  // equation is also treated elsewhere.
  //
  // Integration in time is performed using a backward Euler scheme (in the mass balance
  // equation LHS and RHS are multiplied by the timestep).
  //
  // Using a weak formulation of the governing equation the following terms are assembled in this kernel
  //
  //   Rmom = - \int symmetricGradient( \eta ) : totalStress + \int \eta \cdot bodyForce = 0
  //   Rmas = \int \chi ( fluidMass - fluidMass_n) = 0
  //
  //   dRmom_dVolStrain = - \int_Omega symmetricGradient( \eta ) : dTotalStress_dVolStrain
  //                      + \int \eta \cdot dBodyForce_dVolStrain
  //   dRmom_dPressure  = - \int_Omega symmetricGradient( \eta ) : dTotalStress_dPressure
  //                      + \int \eta \cdot dBodyForce_dPressure
  //   dRmas_dVolStrain = \int \chi dFluidMass_dVolStrain
  //   dRmas_dPressure  = \int \chi dFluidMass_dPressure
  //
  // with \eta and \chi test basis functions for the displacement and pressure field, respectively.
  // A continuous interpolation is used for the displacement, with \eta continuous finite element
  // basis functions. A piecewise-constant approximation is used for the pressure.

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
  assembleElementBasedFlowTerms( dNdX, detJxW, stack );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  GEOS_UNUSED_VAR( k );
  real64 maxForce = 0;
  localIndex const numSupportPoints =
    m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
  integer numDisplacementDofs = numSupportPoints * numDofPerTestSupportPoint;

  for( int localNode = 0; localNode < numSupportPoints; ++localNode )
  {
    for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
    {

      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint*localNode + dim] - m_dofRankOffset );

      // we need this check to filter out ghost nodes in the assembly
      if( dof < 0 || dof >= m_matrix.numRows() )
      {
        continue;
      }
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localRowDofIndex,
                                                                              stack.dLocalResidualMomentum_dDisplacement[numDofPerTestSupportPoint * localNode + dim],
                                                                              numDisplacementDofs );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] );
      maxForce = fmax( maxForce, fabs( stack.localResidualMomentum[numDofPerTestSupportPoint * localNode + dim] ) );
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              &stack.localPressureDofIndex,
                                                                              stack.dLocalResidualMomentum_dPressure[numDofPerTestSupportPoint * localNode + dim],
                                                                              1 );
    }
  }

  localIndex const dof = LvArray::integerConversion< localIndex >( stack.localPressureDofIndex - m_dofRankOffset );

  // we need this check to filter out ghost cells in the assembly
  if( 0 <= dof && dof < m_matrix.numRows() )
  {
    m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( dof,
                                                                    stack.localRowDofIndex,
                                                                    stack.dLocalResidualMass_dDisplacement[0],
                                                                    numDisplacementDofs );
    m_matrix.template addToRow< serialAtomic >( dof,
                                                &stack.localPressureDofIndex,
                                                stack.dLocalResidualMass_dPressure[0],
                                                1 );
    RAJA::atomicAdd< serialAtomic >( &m_rhs[dof], stack.localResidualMass[0] );
  }
  return maxForce;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
SinglePhasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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


} // namespace poromechanicsKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_IMPL_HPP_
