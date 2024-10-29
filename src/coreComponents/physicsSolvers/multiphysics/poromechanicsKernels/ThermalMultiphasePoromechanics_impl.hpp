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
 * @file ThermalMultiphasePoromechanics_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALMULTIPHASEPOROMECHANICS_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALMULTIPHASEPOROMECHANICS_IMPL_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalMultiphasePoromechanics.hpp"

namespace geos
{

namespace thermalPoromechanicsKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
ThermalMultiphasePoromechanics( NodeManager const & nodeManager,
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
                                localIndex const numComponents,
                                localIndex const numPhases,
                                integer const useTotalMassEquation,
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
        numComponents,
        numPhases,
        false, // do not use simple accumulation form
        useTotalMassEquation,
        performStressInitialization,
        fluidModelKey ),
  m_rockInternalEnergy_n( inputConstitutiveType.getInternalEnergy_n() ),
  m_rockInternalEnergy( inputConstitutiveType.getInternalEnergy() ),
  m_dRockInternalEnergy_dTemperature( inputConstitutiveType.getDinternalEnergy_dTemperature() ),
  m_temperature_n( elementSubRegion.template getField< fields::flow::temperature_n >() ),
  m_temperature( elementSubRegion.template getField< fields::flow::temperature >() )
{
  // extract fluid constitutive data views
  {
    string const fluidModelName = elementSubRegion.template getReference< string >( fluidModelKey );
    constitutive::MultiFluidBase const & fluid =
      elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelName );

    m_fluidPhaseInternalEnergy_n = fluid.phaseInternalEnergy_n();
    m_fluidPhaseInternalEnergy = fluid.phaseInternalEnergy();
    m_dFluidPhaseInternalEnergy = fluid.dPhaseInternalEnergy();
  }
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  // initialize displacement, pressure, temperature and fluid components degree of freedom
  Base::setup( k, stack );

  stack.localTemperatureDofIndex = stack.localPressureDofIndex + m_numComponents + 1;
  stack.temperature = m_temperature[k];
  stack.deltaTemperatureFromLastStep = m_temperature[k] - m_temperature_n[k];
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

  // Step 1: call the constitutive model to update the total stress, the porosity and their derivatives
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

  // Step 4: compute pore volume constraint
  computePoreVolumeConstraint( k,
                               porosity_n,
                               stack );
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computeBodyForce( localIndex const k,
                  localIndex const q,
                  real64 const & porosity,
                  real64 const & dPorosity_dVolStrain,
                  real64 const & dPorosity_dPressure,
                  real64 const & dPorosity_dTemperature,
                  real64 const & dSolidDensity_dPressure,
                  StackVariables & stack ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  Base::computeBodyForce( k, q,
                          porosity,
                          dPorosity_dVolStrain,
                          dPorosity_dPressure,
                          dPorosity_dTemperature,
                          dSolidDensity_dPressure,
                          stack, [&]( real64 const & totalMassDensity,
                                      real64 const & mixtureDensity )
  {
    GEOS_UNUSED_VAR( mixtureDensity );

    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseMassDensity = m_fluidPhaseMassDensity[k][q];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseMassDensity = m_dFluidPhaseMassDensity[k][q];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_fluidPhaseVolFrac[k];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];

    // Step 1: compute fluid total mass density and its derivatives

    real64 dTotalMassDensity_dTemperature = 0.0;

    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      dTotalMassDensity_dTemperature += dPhaseVolFrac( ip, Deriv::dT ) * phaseMassDensity( ip )
                                        + phaseVolFrac( ip ) * dPhaseMassDensity( ip, Deriv::dT );
    }

    // Step 2: compute the derivative of the bulk density (an average between total mass density and solid density) wrt temperature

    real64 const dMixtureDens_dTemperature = dPorosity_dTemperature * ( -m_solidDensity( k, q ) + totalMassDensity )
                                             + porosity * dTotalMassDensity_dTemperature;

    // Step 3: finally, get the body force

    LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dTemperature, m_gravityVector, dMixtureDens_dTemperature );

  } );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computeFluidIncrement( localIndex const k,
                       localIndex const q,
                       real64 const & porosity,
                       real64 const & porosity_n,
                       real64 const & dPorosity_dVolStrain,
                       real64 const & dPorosity_dPressure,
                       real64 const & dPorosity_dTemperature,
                       StackVariables & stack ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  stack.energyIncrement = 0.0;
  stack.dEnergyIncrement_dVolStrainIncrement = 0.0;
  stack.dEnergyIncrement_dPressure = 0.0;
  stack.dEnergyIncrement_dTemperature = 0.0;
  LvArray::tensorOps::fill< maxNumComponents >( stack.dEnergyIncrement_dComponents, 0.0 );
  LvArray::tensorOps::fill< maxNumComponents >( stack.dCompMassIncrement_dTemperature, 0.0 );

  Base::computeFluidIncrement( k, q,
                               porosity,
                               porosity_n,
                               dPorosity_dVolStrain,
                               dPorosity_dPressure,
                               dPorosity_dTemperature,
                               stack, [&]( real64 const & ip,
                                           real64 const & phaseAmount,
                                           real64 const & phaseAmount_n,
                                           real64 const & dPhaseAmount_dVolStrain,
                                           real64 const & dPhaseAmount_dP,
                                           real64 const (&dPhaseAmount_dC)[maxNumComponents] )
  {

    // We are in the loop over phases, ip provides the current phase index.
    // We have to do two things:
    //   1- Assemble the derivatives of the component mass balance equations with respect to temperature
    //   2- Assemble the phase-dependent part of the accumulation term of the energy equation

    real64 dPhaseInternalEnergy_dC[maxNumComponents]{};

    // construct the slices
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseInternalEnergy_n = m_fluidPhaseInternalEnergy_n[k][q];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseInternalEnergy = m_fluidPhaseInternalEnergy[k][q];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseInternalEnergy = m_dFluidPhaseInternalEnergy[k][q];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseDensity = m_fluidPhaseDensity[k][q];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseDensity = m_dFluidPhaseDensity[k][q];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const phaseCompFrac = m_fluidPhaseCompFrac[k][q];
    arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC -2 > const dPhaseCompFrac = m_dFluidPhaseCompFrac[k][q];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_fluidPhaseVolFrac[k];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dGlobalCompFrac_dGlobalCompDensity = m_dGlobalCompFraction_dGlobalCompDensity[k];

    // Step 1: assemble the derivatives of the component mass balance equations wrt temperature

    real64 const dPhaseAmount_dT = dPorosity_dTemperature * phaseVolFrac( ip ) * phaseDensity( ip )
                                   + porosity * ( dPhaseVolFrac( ip, Deriv::dT ) * phaseDensity( ip ) + phaseVolFrac( ip ) * dPhaseDensity( ip, Deriv::dT ) );

    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      stack.dCompMassIncrement_dTemperature[ic] += dPhaseAmount_dT * phaseCompFrac( ip, ic )
                                                   + phaseAmount * dPhaseCompFrac( ip, ic, Deriv::dT );
    }

    // Step 2: assemble the phase-dependent part of the accumulation term of the energy equation

    // local accumulation
    stack.energyIncrement += phaseAmount * phaseInternalEnergy( ip ) - phaseAmount_n * phaseInternalEnergy_n( ip );

    // derivatives w.r.t. vol strain increment, pressure and temperature
    stack.dEnergyIncrement_dVolStrainIncrement += dPhaseAmount_dVolStrain * phaseInternalEnergy( ip );
    stack.dEnergyIncrement_dPressure += dPhaseAmount_dP * phaseInternalEnergy( ip )
                                        + phaseAmount * dPhaseInternalEnergy( ip, Deriv::dP );
    stack.dEnergyIncrement_dTemperature += dPhaseAmount_dT * phaseInternalEnergy( ip )
                                           + phaseAmount * dPhaseInternalEnergy( ip, Deriv::dT );

    // derivatives w.r.t. component densities
    applyChainRule( m_numComponents, dGlobalCompFrac_dGlobalCompDensity, dPhaseInternalEnergy[ip], dPhaseInternalEnergy_dC, Deriv::dC );
    for( integer jc = 0; jc < m_numComponents; ++jc )
    {
      stack.dEnergyIncrement_dComponents[jc] += dPhaseAmount_dC[jc] * phaseInternalEnergy( ip )
                                                + phaseAmount * dPhaseInternalEnergy_dC[jc];
    }
  } );

  // Step 3: assemble the solid part of the accumulation term
  real64 const oneMinusPoro = 1 - porosity;

  // local accumulation and derivatives w.r.t. pressure and temperature
  stack.energyIncrement += oneMinusPoro * m_rockInternalEnergy( k, 0 ) - ( 1 - porosity_n ) * m_rockInternalEnergy_n( k, 0 );
  stack.dEnergyIncrement_dVolStrainIncrement += -dPorosity_dVolStrain * m_rockInternalEnergy( k, 0 );
  stack.dEnergyIncrement_dPressure += -dPorosity_dPressure * m_rockInternalEnergy( k, 0 );
  stack.dEnergyIncrement_dTemperature += -dPorosity_dTemperature * m_rockInternalEnergy( k, 0 ) + oneMinusPoro * m_dRockInternalEnergy_dTemperature( k, 0 );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computePoreVolumeConstraint( localIndex const k,
                             real64 const & porosity_n,
                             StackVariables & stack ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  Base::computePoreVolumeConstraint( k,
                                     porosity_n,
                                     stack );

  arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];

  stack.dPoreVolConstraint_dTemperature = 0.0;
  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    stack.dPoreVolConstraint_dTemperature -= dPhaseVolFrac( ip, Deriv::dT ) * porosity_n;
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
assembleMomentumBalanceTerms( real64 const ( &N )[numNodesPerElem],
                              real64 const ( &dNdX )[numNodesPerElem][3],
                              real64 const & detJxW,
                              StackVariables & stack ) const
{
  using namespace PDEUtilities;

  Base::assembleMomentumBalanceTerms( N, dNdX, detJxW, stack );

  constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
  constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;
  constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;

  // compute local linear momentum balance residual derivatives with respect to temperature

  BilinearFormUtilities::compute< displacementTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::SymmetricGradient,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMomentum_dTemperature,
    dNdX,
    stack.dTotalStress_dTemperature,
    1.0,
    -detJxW );

  BilinearFormUtilities::compute< displacementTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMomentum_dTemperature,
    N,
    stack.dBodyForce_dTemperature,
    1.0,
    detJxW );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
assembleElementBasedFlowTerms( real64 const ( &dNdX )[numNodesPerElem][3],
                               real64 const & detJxW,
                               StackVariables & stack ) const
{
  using namespace PDEUtilities;

  Base::assembleElementBasedFlowTerms( dNdX, detJxW, stack );

  constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
  constexpr FunctionSpace pressureTrialSpace = FunctionSpace::P0;
  constexpr FunctionSpace pressureTestSpace = pressureTrialSpace;

  // Step 1: compute local mass balance residual derivatives with respect to temperature

  BilinearFormUtilities::compute< pressureTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMass_dTemperature,
    1.0,
    stack.dCompMassIncrement_dTemperature,
    1.0,
    detJxW );

  // Step 2: compute local pore volume constraint residual derivatives with respect to temperature

  BilinearFormUtilities::compute< pressureTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualPoreVolConstraint_dTemperature,
    1.0,
    stack.dPoreVolConstraint_dTemperature,
    1.0,
    detJxW );

  // Step 3: compute local energy balance residual and its derivatives

  LinearFormUtilities::compute< pressureTestSpace,
                                DifferentialOperator::Identity >
  (
    stack.localResidualEnergy,
    1.0,
    stack.energyIncrement,
    detJxW );

  BilinearFormUtilities::compute< pressureTestSpace,
                                  displacementTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Divergence >
  (
    stack.dLocalResidualEnergy_dDisplacement,
    1.0,
    stack.dEnergyIncrement_dVolStrainIncrement,
    dNdX,
    detJxW );

  BilinearFormUtilities::compute< pressureTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualEnergy_dPressure,
    1.0,
    stack.dEnergyIncrement_dPressure,
    1.0,
    detJxW );

  BilinearFormUtilities::compute< pressureTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualEnergy_dTemperature,
    1.0,
    stack.dEnergyIncrement_dTemperature,
    1.0,
    detJxW );

  BilinearFormUtilities::compute< pressureTestSpace,
                                  pressureTrialSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualEnergy_dComponents,
    1.0,
    stack.dEnergyIncrement_dComponents,
    1.0,
    detJxW );

}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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
  assembleElementBasedFlowTerms( dNdX, detJxW, stack );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  using namespace compositionalMultiphaseUtilities;

  real64 const maxForce = Base::complete( k, stack );

  localIndex const numSupportPoints =
    m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
  integer numDisplacementDofs = numSupportPoints * numDofPerTestSupportPoint;

  if( m_useTotalMassEquation > 0 )
  {
    // Apply equation/variable change transformation(s)
    real64 work[maxNumComponents + 1]{};
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, 1, stack.dLocalResidualMass_dTemperature, work );
  }

  // Step 1: assemble the derivatives of linear momentum balance wrt temperature into the global matrix

  for( int localNode = 0; localNode < numSupportPoints; ++localNode )
  {
    for( integer dim = 0; dim < numDofPerTestSupportPoint; ++dim )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint*localNode + dim] - m_dofRankOffset );

      // we need this check to filter out ghost nodes in the assembly
      if( dof < 0 || dof >= m_matrix.numRows() )
      {
        continue;
      }

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              &stack.localTemperatureDofIndex,
                                                                              stack.dLocalResidualMomentum_dTemperature[numDofPerTestSupportPoint * localNode + dim],
                                                                              1 );
    }
  }

  localIndex const massDof = LvArray::integerConversion< localIndex >( stack.localPressureDofIndex - m_dofRankOffset );

  // we need this check to filter out ghost nodes in the assembly
  if( 0 <= massDof && massDof < m_matrix.numRows() )
  {

    // Step 2: assemble the derivatives of mass balance residual wrt temperature into the global matrix

    for( localIndex i = 0; i < m_numComponents; ++i )
    {
      m_matrix.template addToRow< serialAtomic >( massDof + i,
                                                  &stack.localTemperatureDofIndex,
                                                  stack.dLocalResidualMass_dTemperature[i],
                                                  1 );
    }

    // Step 3: assemble the derivatives of pore-volume constraint residual wrt temperature into the global matrix

    m_matrix.template addToRow< serialAtomic >( massDof + m_numComponents,
                                                &stack.localTemperatureDofIndex,
                                                stack.dLocalResidualPoreVolConstraint_dTemperature[0],
                                                1 );
  }

  // Step 4: assemble the energy balance and its derivatives into the global matrix

  localIndex const energyDof = LvArray::integerConversion< localIndex >( stack.localTemperatureDofIndex - m_dofRankOffset );

  // we need this check to filter out ghost nodes in the assembly
  if( 0 <= energyDof && energyDof < m_matrix.numRows() )
  {
    m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( energyDof,
                                                                    stack.localRowDofIndex,
                                                                    stack.dLocalResidualEnergy_dDisplacement[0],
                                                                    numDisplacementDofs );
    m_matrix.template addToRow< serialAtomic >( energyDof,
                                                &stack.localPressureDofIndex,
                                                stack.dLocalResidualEnergy_dPressure[0],
                                                1 );
    m_matrix.template addToRow< serialAtomic >( energyDof,
                                                &stack.localTemperatureDofIndex,
                                                stack.dLocalResidualEnergy_dTemperature[0],
                                                1 );
    m_matrix.template addToRow< serialAtomic >( energyDof,
                                                stack.localComponentDofIndices,
                                                stack.dLocalResidualEnergy_dComponents[0],
                                                m_numComponents );

    RAJA::atomicAdd< serialAtomic >( &m_rhs[energyDof], stack.localResidualEnergy[0] );
  }

  return maxForce;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64 ThermalMultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

} // namespace thermalporomechanicsKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALMULTIPHASEPOROMECHANICS_IMPL_HPP_
