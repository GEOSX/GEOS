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
 * @file MultiphasePoromechanics_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICS_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICS_IMPL_HPP_

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp"

namespace geos
{

namespace poromechanicsKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
MultiphasePoromechanics( NodeManager const & nodeManager,
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
                         integer const useSimpleAccumulation,
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
        fluidModelKey ),
  m_fluidPhaseVolFrac( elementSubRegion.template getField< fields::flow::phaseVolumeFraction >() ),
  m_fluidPhaseVolFrac_n( elementSubRegion.template getField< fields::flow::phaseVolumeFraction_n >() ),
  m_dFluidPhaseVolFrac( elementSubRegion.template getField< fields::flow::dPhaseVolumeFraction >() ),
  m_dGlobalCompFraction_dGlobalCompDensity( elementSubRegion.template getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
  m_compDens( elementSubRegion.template getField< fields::flow::globalCompDensity >() ),
  m_compDens_n( elementSubRegion.template getField< fields::flow::globalCompDensity_n >() ),
  m_numComponents( numComponents ),
  m_numPhases( numPhases ),
  m_useSimpleAccumulation( useSimpleAccumulation ),
  m_useTotalMassEquation( useTotalMassEquation ),
  m_performStressInitialization( performStressInitialization )
{
  GEOS_ERROR_IF_GT_MSG( m_numComponents, maxNumComponents,
                        "MultiphasePoromechanics solver allows at most " <<
                        maxNumComponents << " components at the moment" );

  // extract fluid constitutive data views
  {
    string const fluidModelName = elementSubRegion.template getReference< string >( fluidModelKey );
    constitutive::MultiFluidBase const & fluid =
      elementSubRegion.template getConstitutiveModel< constitutive::MultiFluidBase >( fluidModelName );

    m_fluidPhaseDensity = fluid.phaseDensity();
    m_fluidPhaseDensity_n = fluid.phaseDensity_n();
    m_dFluidPhaseDensity = fluid.dPhaseDensity();

    m_fluidPhaseCompFrac = fluid.phaseCompFraction();
    m_fluidPhaseCompFrac_n = fluid.phaseCompFraction_n();
    m_dFluidPhaseCompFrac = fluid.dPhaseCompFraction();

    m_fluidPhaseMassDensity = fluid.phaseMassDensity();
    m_dFluidPhaseMassDensity = fluid.dPhaseMassDensity();
  }
}

/**
 * @brief Copy global values from primary field to a local stack array.
 * @copydoc ::geos::finiteElement::ImplicitKernelBase::setup
 *
 * For the MultiphasePoromechanics implementation, global values from the displacement,
 * incremental displacement, and degree of freedom numbers are placed into
 * element local stack storage.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  // initialize displacement dof and pressure dofs
  Base::setup( k, stack );

  // setup component dofs
  // for now, maxNumComponents > m_numComponents, so we pad localComponentDofIndices with -1
  LvArray::tensorOps::fill< maxNumComponents >( stack.localComponentDofIndices, -1.0 );
  for( integer flowDofIndex=0; flowDofIndex < m_numComponents; ++flowDofIndex )
  {
    stack.localComponentDofIndices[flowDofIndex] = stack.localPressureDofIndex + flowDofIndex + 1;
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

  // Step 4: compute pore volume constraint
  computePoreVolumeConstraint( k,
                               porosity_n,
                               stack );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computeBodyForce( localIndex const k,
                  localIndex const q,
                  real64 const & porosity,
                  real64 const & dPorosity_dVolStrain,
                  real64 const & dPorosity_dPressure,
                  real64 const & dPorosity_dTemperature,
                  real64 const & dSolidDensity_dPressure,
                  StackVariables & stack,
                  FUNC && bodyForceKernelOp ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  GEOS_UNUSED_VAR( dPorosity_dTemperature );

  arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const phaseMassDensity = m_fluidPhaseMassDensity[k][q];
  arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const dPhaseMassDensity = m_dFluidPhaseMassDensity[k][q];
  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_fluidPhaseVolFrac[k];
  arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];
  arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const dGlobalCompFrac_dGlobalCompDensity = m_dGlobalCompFraction_dGlobalCompDensity[k];

  // Step 1: compute fluid total mass density and its derivatives

  real64 totalMassDensity = 0.0;
  real64 dTotalMassDensity_dPressure = 0.0;
  real64 dTotalMassDensity_dComponents[maxNumComponents]{};
  real64 dPhaseMassDensity_dComponents[maxNumComponents]{};

  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    totalMassDensity += phaseVolFrac( ip ) * phaseMassDensity( ip );
    dTotalMassDensity_dPressure += dPhaseVolFrac( ip, Deriv::dP ) * phaseMassDensity( ip )
                                   + phaseVolFrac( ip ) * dPhaseMassDensity( ip, Deriv::dP );

    applyChainRule( m_numComponents,
                    dGlobalCompFrac_dGlobalCompDensity,
                    dPhaseMassDensity[ip],
                    dPhaseMassDensity_dComponents,
                    Deriv::dC );
    for( integer jc = 0; jc < m_numComponents; ++jc )
    {
      dTotalMassDensity_dComponents[jc] += dPhaseVolFrac( ip, Deriv::dC+jc ) * phaseMassDensity( ip )
                                           + phaseVolFrac( ip ) * dPhaseMassDensity_dComponents[jc];
    }
  }

  // Step 2: compute mixture density as an average between total mass density and solid density

  real64 const mixtureDensity = ( 1.0 - porosity ) * m_solidDensity( k, q ) + porosity * totalMassDensity;
  real64 const dMixtureDens_dVolStrainIncrement = dPorosity_dVolStrain * ( -m_solidDensity( k, q ) + totalMassDensity );
  real64 const dMixtureDens_dPressure = dPorosity_dPressure * ( -m_solidDensity( k, q ) + totalMassDensity )
                                        + ( 1.0 - porosity ) * dSolidDensity_dPressure
                                        + porosity * dTotalMassDensity_dPressure;
  LvArray::tensorOps::scale< maxNumComponents >( dTotalMassDensity_dComponents, porosity );

  // Step 3: finally, get the body force

  LvArray::tensorOps::scaledCopy< 3 >( stack.bodyForce, m_gravityVector, mixtureDensity );
  LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dVolStrainIncrement, m_gravityVector, dMixtureDens_dVolStrainIncrement );
  LvArray::tensorOps::scaledCopy< 3 >( stack.dBodyForce_dPressure, m_gravityVector, dMixtureDens_dPressure );
  LvArray::tensorOps::Rij_eq_AiBj< 3, maxNumComponents >( stack.dBodyForce_dComponents, m_gravityVector, dTotalMassDensity_dComponents );

  // Step 4: customize the kernel (for instance, to add thermal derivatives)
  bodyForceKernelOp( totalMassDensity, mixtureDensity );

}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computeFluidIncrement( localIndex const k,
                       localIndex const q,
                       real64 const & porosity,
                       real64 const & porosity_n,
                       real64 const & dPorosity_dVolStrain,
                       real64 const & dPorosity_dPressure,
                       real64 const & dPorosity_dTemperature,
                       StackVariables & stack,
                       FUNC && fluidIncrementKernelOp ) const
{
  LvArray::tensorOps::fill< maxNumComponents >( stack.compMassIncrement, 0.0 );
  LvArray::tensorOps::fill< maxNumComponents >( stack.dCompMassIncrement_dVolStrainIncrement, 0.0 );
  LvArray::tensorOps::fill< maxNumComponents >( stack.dCompMassIncrement_dPressure, 0.0 );
  LvArray::tensorOps::fill< maxNumComponents, maxNumComponents >( stack.dCompMassIncrement_dComponents, 0.0 );

  if( m_useSimpleAccumulation )
  {
    arraySlice1d< real64 const, compflow::USD_COMP - 1 >
    const compDens = m_compDens[k];
    arraySlice1d< real64 const, compflow::USD_COMP - 1 >
    const compDens_n = m_compDens_n[k];

    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      stack.compMassIncrement[ic] += porosity * compDens[ic] - porosity_n * compDens_n[ic];
      stack.dCompMassIncrement_dPressure[ic] += dPorosity_dPressure * compDens[ic];
      stack.dCompMassIncrement_dVolStrainIncrement[ic] += dPorosity_dVolStrain * compDens[ic];
      stack.dCompMassIncrement_dComponents[ic][ic] += porosity;
    }
  }
  else
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    GEOS_UNUSED_VAR( dPorosity_dTemperature );

    // temporary work arrays and slices
    real64 dPhaseAmount_dC[maxNumComponents]{};
    real64 dPhaseCompFrac_dC[maxNumComponents]{};

    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 >
    const phaseDensity = m_fluidPhaseDensity[k][q];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 >
    const phaseDensity_n = m_fluidPhaseDensity_n[k][q];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 >
    const dPhaseDensity = m_dFluidPhaseDensity[k][q];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 >
    const phaseCompFrac = m_fluidPhaseCompFrac[k][q];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 >
    const phaseCompFrac_n = m_fluidPhaseCompFrac_n[k][q];
    arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 2 >
    const dPhaseCompFrac = m_dFluidPhaseCompFrac[k][q];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 >
    const phaseVolFrac = m_fluidPhaseVolFrac[k];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 >
    const phaseVolFrac_n = m_fluidPhaseVolFrac_n[k];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 >
    const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 >
    const dGlobalCompFrac_dGlobalCompDensity = m_dGlobalCompFraction_dGlobalCompDensity[k];

    // loop over the fluid phases
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {

      // compute the mass of the current phase
      real64 const phaseAmount = porosity * phaseVolFrac( ip ) * phaseDensity( ip );
      real64 const phaseAmount_n = porosity_n * phaseVolFrac_n( ip ) * phaseDensity_n( ip );

      real64 const dPhaseAmount_dVolStrain = dPorosity_dVolStrain * phaseVolFrac( ip ) * phaseDensity( ip );
      real64 const dPhaseAmount_dP = dPorosity_dPressure * phaseVolFrac( ip ) * phaseDensity( ip )
                                     + porosity * ( dPhaseVolFrac( ip, Deriv::dP ) * phaseDensity( ip )
                                                    + phaseVolFrac( ip ) * dPhaseDensity( ip, Deriv::dP ) );

      applyChainRule( m_numComponents,
                      dGlobalCompFrac_dGlobalCompDensity,
                      dPhaseDensity[ip],
                      dPhaseAmount_dC,
                      Deriv::dC );

      for( integer jc = 0; jc < m_numComponents; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac( ip )
                              + phaseDensity( ip ) * dPhaseVolFrac( ip, Deriv::dC + jc );
        dPhaseAmount_dC[jc] *= porosity;
      }

      // for each phase, compute the amount of each component transported by the phase
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        stack.compMassIncrement[ic] += phaseAmount * phaseCompFrac( ip, ic )
                                       - phaseAmount_n * phaseCompFrac_n( ip, ic );

        stack.dCompMassIncrement_dPressure[ic] += dPhaseAmount_dP * phaseCompFrac( ip, ic )
                                                  + phaseAmount * dPhaseCompFrac( ip, ic, Deriv::dP );
        stack.dCompMassIncrement_dVolStrainIncrement[ic] += dPhaseAmount_dVolStrain * phaseCompFrac( ip, ic );

        applyChainRule( m_numComponents,
                        dGlobalCompFrac_dGlobalCompDensity,
                        dPhaseCompFrac[ip][ic],
                        dPhaseCompFrac_dC,
                        Deriv::dC );

        for( integer jc = 0; jc < m_numComponents; ++jc )
        {
          stack.dCompMassIncrement_dComponents[ic][jc] += dPhaseAmount_dC[jc] * phaseCompFrac( ip, ic )
                                                          + phaseAmount * dPhaseCompFrac_dC[jc];
        }
      }

      // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
      // possible use: assemble the derivatives wrt temperature, and the accumulation term of the energy equation for this phase
      fluidIncrementKernelOp( ip, phaseAmount, phaseAmount_n, dPhaseAmount_dVolStrain, dPhaseAmount_dP, dPhaseAmount_dC );

    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
computePoreVolumeConstraint( localIndex const k,
                             real64 const & porosity_n,
                             StackVariables & stack ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const phaseVolFrac = m_fluidPhaseVolFrac[k];
  arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const dPhaseVolFrac = m_dFluidPhaseVolFrac[k];

  stack.poreVolConstraint = 1.0;
  stack.dPoreVolConstraint_dPressure = 0.0;
  LvArray::tensorOps::fill< 1, maxNumComponents >( stack.dPoreVolConstraint_dComponents, 0.0 );

  for( integer ip = 0; ip < m_numPhases; ++ip )
  {
    stack.poreVolConstraint -= phaseVolFrac( ip );
    stack.dPoreVolConstraint_dPressure -= dPhaseVolFrac( ip, Deriv::dP ) * porosity_n;

    for( integer jc = 0; jc < m_numComponents; ++jc )
    {
      stack.dPoreVolConstraint_dComponents[0][jc] -= dPhaseVolFrac( ip, Deriv::dC+jc ) * porosity_n;
    }
  }
  stack.poreVolConstraint *= porosity_n;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

  // Step 4: compute local linear momentum balance residual derivatives with respect to components
  BilinearFormUtilities::compute< displacementTestSpace,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMomentum_dComponents,
    N,
    stack.dBodyForce_dComponents,
    1.0,
    detJxW );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
assembleElementBasedFlowTerms( real64 const ( &dNdX )[numNodesPerElem][3],
                               real64 const & detJxW,
                               StackVariables & stack ) const
{
  using namespace PDEUtilities;

  constexpr FunctionSpace displacementTrialSpace = FE_TYPE::template getFunctionSpace< numDofPerTrialSupportPoint >();
  constexpr FunctionSpace displacementTestSpace = displacementTrialSpace;

  // Step 1: mass balance equations

  // compute local component mass balance residual
  LinearFormUtilities::compute< FunctionSpace::P0,
                                DifferentialOperator::Identity >
  (
    stack.localResidualMass,
    1.0,
    stack.compMassIncrement,
    detJxW );

  // compute local mass balance residual derivatives with respect to displacement
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  displacementTestSpace,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Divergence >
  (
    stack.dLocalResidualMass_dDisplacement,
    1.0,
    stack.dCompMassIncrement_dVolStrainIncrement,
    dNdX,
    detJxW );

  // compute local mass balance residual derivatives with respect to pressure
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMass_dPressure,
    1.0,
    stack.dCompMassIncrement_dPressure,
    1.0,
    detJxW );

  // compute local mass balance residual derivatives with respect to components
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualMass_dComponents,
    1.0,
    stack.dCompMassIncrement_dComponents,
    1.0,
    detJxW );


  // Step 2: pore volume constraint equation

  // compute local pore volume contraint residual
  LinearFormUtilities::compute< FunctionSpace::P0,
                                DifferentialOperator::Identity >
  (
    stack.localResidualPoreVolConstraint,
    1.0,
    stack.poreVolConstraint,
    detJxW );

  // compute local pore volume contraint residual derivatives with respect to pressure
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualPoreVolConstraint_dPressure,
    1.0,
    stack.dPoreVolConstraint_dPressure,
    1.0,
    detJxW );

  // compute local pore volume contraint residual derivatives with respect to components
  BilinearFormUtilities::compute< FunctionSpace::P0,
                                  FunctionSpace::P0,
                                  DifferentialOperator::Identity,
                                  DifferentialOperator::Identity >
  (
    stack.dLocalResidualPoreVolConstraint_dComponents,
    1.0,
    stack.dPoreVolConstraint_dComponents,
    1.0,
    detJxW );
}


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

/**
 * @copydoc geos::finiteElement::ImplicitKernelBase::complete
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  using namespace compositionalMultiphaseUtilities;

  GEOS_UNUSED_VAR( k );

  real64 maxForce = 0;
  localIndex const numSupportPoints =
    m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );
  integer numDisplacementDofs = numSupportPoints * numDofPerTestSupportPoint;
  constexpr integer maxNumDisplacementDofs = FE_TYPE::maxSupportPoints * numDofPerTestSupportPoint;

  if( m_useTotalMassEquation > 0 )
  {
    // Apply equation/variable change transformation(s)
    real64 work[maxNumDisplacementDofs > ( maxNumComponents + 1 ) ? maxNumDisplacementDofs : maxNumComponents + 1];
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, numDisplacementDofs, stack.dLocalResidualMass_dDisplacement, work );
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, 1, stack.dLocalResidualMass_dPressure, work );
    shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( m_numComponents, m_numComponents, stack.dLocalResidualMass_dComponents, work );
    shiftElementsAheadByOneAndReplaceFirstElementWithSum( m_numComponents, stack.localResidualMass );
  }

  for( int localNode = 0; localNode < numSupportPoints; ++localNode )
  {
    for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint * localNode + dim] - m_dofRankOffset );

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

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localComponentDofIndices,
                                                                              stack.dLocalResidualMomentum_dComponents[numDofPerTestSupportPoint * localNode + dim],
                                                                              m_numComponents );
    }
  }

  localIndex const dof = LvArray::integerConversion< localIndex >( stack.localPressureDofIndex - m_dofRankOffset );

  // we need this check to filter out ghost cells in the assembly
  if( 0 <= dof && dof < m_matrix.numRows() )
  {
    for( localIndex i = 0; i < m_numComponents; ++i )
    {
      m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( dof + i,
                                                                      stack.localRowDofIndex,
                                                                      stack.dLocalResidualMass_dDisplacement[i],
                                                                      numDisplacementDofs );
      m_matrix.template addToRow< serialAtomic >( dof + i,
                                                  &stack.localPressureDofIndex,
                                                  stack.dLocalResidualMass_dPressure[i],
                                                  1 );
      m_matrix.template addToRow< serialAtomic >( dof + i,
                                                  stack.localComponentDofIndices,
                                                  stack.dLocalResidualMass_dComponents[i],
                                                  m_numComponents );
      RAJA::atomicAdd< serialAtomic >( &m_rhs[dof+i], stack.localResidualMass[i] );
    }

    m_matrix.template addToRow< serialAtomic >( dof + m_numComponents,
                                                &stack.localPressureDofIndex,
                                                stack.dLocalResidualPoreVolConstraint_dPressure[0],
                                                1 );

    m_matrix.template addToRow< serialAtomic >( dof + m_numComponents,
                                                stack.localComponentDofIndices,
                                                stack.dLocalResidualPoreVolConstraint_dComponents[0],
                                                m_numComponents );

    RAJA::atomicAdd< serialAtomic >( &m_rhs[dof+m_numComponents], stack.localResidualPoreVolConstraint[0] );
  }
  return maxForce;
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64 MultiphasePoromechanics< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
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

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_MULTIPHASEPOROMECHANICS_IMPL_HPP_
