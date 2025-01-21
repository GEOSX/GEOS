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
 * @file SinglePhaseWellKernels.cpp
 */

#include "SinglePhaseWellKernels.hpp"

// TODO: move keys to WellControls
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellFields.hpp"
namespace geos
{

namespace singlePhaseWellKernels
{

/******************************** ControlEquationHelper ********************************/

GEOS_HOST_DEVICE
inline
void
ControlEquationHelper::
  switchControl( bool const isProducer,
                 WellControls::Control const & currentControl,
                 real64 const & targetBHP,
                 real64 const & targetRate,
                 real64 const & currentBHP,
                 real64 const & currentVolRate,
                 WellControls::Control & newControl )
{
  // if isViable is true at the end of the following checks, no need to switch
  bool controlIsViable = false;

  // The limiting flow rates are treated as upper limits, while the pressure limits
  // are treated as lower limits in production wells and upper limits in injectors.
  // The well changes its mode of control whenever the existing control mode would
  // violate one of these limits.
  // BHP control
  if( currentControl == WellControls::Control::BHP )
  {
    // the control is viable if the reference rate is below the max rate
    controlIsViable = ( LvArray::math::abs( currentVolRate ) <= LvArray::math::abs( targetRate ) + EPS );
  }
  else // rate control
  {
    // the control is viable if the reference pressure is below/above the max/min pressure
    if( isProducer )
    {
      // targetBHP specifies a min pressure here
      controlIsViable = ( currentBHP >= targetBHP - EPS );
    }
    else
    {
      // targetBHP specifies a max pressure here
      controlIsViable = ( currentBHP <= targetBHP + EPS );
    }
  }

  if( controlIsViable )
  {
    newControl = currentControl;
  }
  else
  {
    // Note: if BHP control is not viable, we switch to TOTALVOLRATE
    //       if TOTALVOLRATE are not viable, we switch to BHP
    newControl = ( currentControl == WellControls::Control::BHP )
               ? WellControls::Control::TOTALVOLRATE
               : WellControls::Control::BHP;
  }
}

template< integer IS_THERMAL >
GEOS_HOST_DEVICE
inline
void
ControlEquationHelper::
  compute( globalIndex const rankOffset,
           WellControls::Control const currentControl,
           real64 const & targetBHP,
           real64 const & targetRate,
           real64 const & currentBHP,
           real64 const & dCurrentBHP_dPres,
           real64 const & currentVolRate,
           real64 const & dCurrentVolRate_dPres,
           real64 const & dCurrentVolRate_dRate,
           globalIndex const dofNumber,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{
  localIndex const eqnRowIndex = dofNumber + ROFFSET::CONTROL - rankOffset;
  globalIndex const presDofColIndex = dofNumber + COFFSET::DPRES;
  globalIndex const rateDofColIndex = dofNumber + COFFSET::DRATE;

  real64 controlEqn = 0;
  real64 dControlEqn_dRate = 0;
  real64 dControlEqn_dPres = 0;

  // Note: We assume in the computation of currentBHP that the reference elevation
  //       is in the top well element. This is enforced by a check in the solver.
  //       If we wanted to allow the reference elevation to be outside the top
  //       well element, it would make more sense to check the BHP constraint in
  //       the well element that contains the reference elevation.

  // BHP control
  if( currentControl == WellControls::Control::BHP )
  {
    // control equation is a difference between current BHP and target BHP
    controlEqn = currentBHP - targetBHP;
    dControlEqn_dPres = dCurrentBHP_dPres;
  }
  // Total volumetric rate control
  else if( currentControl == WellControls::Control::TOTALVOLRATE )
  {
    // control equation is the difference between volumetric current rate and target rate
    controlEqn = currentVolRate - targetRate;
    dControlEqn_dRate = dCurrentVolRate_dRate;
    dControlEqn_dPres = dCurrentVolRate_dPres;
  }
  else
  {
    GEOS_ERROR( "This constraint is not supported in SinglePhaseWell" );
  }

  localRhs[eqnRowIndex] += controlEqn;
  localMatrix.addToRow< serialAtomic >( eqnRowIndex,
                                        &presDofColIndex,
                                        &dControlEqn_dPres,
                                        1 );
  localMatrix.addToRow< serialAtomic >( eqnRowIndex,
                                        &rateDofColIndex,
                                        &dControlEqn_dRate,
                                        1 );
}

/******************************** FluxKernel ********************************/
#define INST_FluxKernel( IS_THERMAL ) \
  template \
  void  \
  FluxKernel::  \
    launch< IS_THERMAL >( localIndex const size,  \
                          globalIndex const rankOffset,  \
                          arrayView1d< globalIndex const > const & wellElemDofNumber,  \
                          arrayView1d< localIndex const > const & nextWellElemIndex,  \
                          arrayView1d< real64 const > const & connRate,  \
                          real64 const & dt,  \
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,  \
                          arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 0 );
INST_FluxKernel( 1 );

template< integer IS_THERMAL >
void
FluxKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  // loop over the well elements to compute the fluxes between elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {

    // 1) Compute the flux and its derivatives

    /*  currentConnRate < 0 flow from iwelem to iwelemNext
     *  currentConnRate > 0 flow from iwelemNext to iwelem
     *  With this convention, currentConnRate < 0 at the last connection for a producer
     *                        currentConnRate > 0 at the last connection for a injector
     */

    // get next well element index
    localIndex const iwelemNext = nextWellElemIndex[iwelem];

    // there is nothing to upwind for single-phase flow
    real64 const currentConnRate = connRate[iwelem];
    real64 const flux = dt * currentConnRate;
    real64 const dFlux_dRate = dt;

    // 2) Assemble the flux into residual and Jacobian
    if( iwelemNext < 0 )
    {
      // flux terms
      real64 const oneSidedLocalFlux = -flux;
      real64 const oneSidedLocalFluxJacobian_dRate = -dFlux_dRate;

      // jacobian indices
      globalIndex const oneSidedEqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - rankOffset;
      globalIndex const oneSidedDofColIndex_dRate = wellElemDofNumber[iwelem] + COFFSET::DRATE;

      if( oneSidedEqnRowIndex >= 0 && oneSidedEqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndex,
                                                      &oneSidedDofColIndex_dRate,
                                                      &oneSidedLocalFluxJacobian_dRate,
                                                      1 );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[oneSidedEqnRowIndex], oneSidedLocalFlux );
      }
    }
    else
    {
      // local working variables and arrays
      globalIndex eqnRowIndices[2]{};

      real64 localFlux[2]{};
      real64 localFluxJacobian_dRate[2]{};

      // flux terms
      localFlux[TAG::NEXT]    =   flux;
      localFlux[TAG::CURRENT] = -flux;

      localFluxJacobian_dRate[TAG::NEXT]    =   dFlux_dRate;
      localFluxJacobian_dRate[TAG::CURRENT] = -dFlux_dRate;

      // indices
      eqnRowIndices[TAG::CURRENT] = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - rankOffset;
      eqnRowIndices[TAG::NEXT]    = wellElemDofNumber[iwelemNext] + ROFFSET::MASSBAL - rankOffset;
      globalIndex const dofColIndex_dRate = wellElemDofNumber[iwelem] + COFFSET::DRATE;

      for( localIndex i = 0; i < 2; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                        &dofColIndex_dRate,
                                                        &localFluxJacobian_dRate[i],
                                                        1 );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localFlux[i] );
        }
      }
    }
  } );
}

/******************************** PressureRelationKernel ********************************/

#define INST_PressureRelationKernel( IS_THERMAL ) \
  template \
  localIndex \
  PressureRelationKernel:: \
    launch< IS_THERMAL >( localIndex const size, \
                          globalIndex const rankOffset, \
                          bool const isLocallyOwned, \
                          localIndex const iwelemControl, \
                          WellControls const & wellControls, \
                          real64 const & timeAtEndOfStep, \
                          arrayView1d< globalIndex const > const & wellElemDofNumber, \
                          arrayView1d< real64 const > const & wellElemGravCoef, \
                          arrayView1d< localIndex const > const & nextWellElemIndex, \
                          arrayView1d< real64 const > const & wellElemPressure, \
                          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity, \
                          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity, \
                          CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                          arrayView1d< real64 > const & localRhs )

INST_PressureRelationKernel( 0 );
INST_PressureRelationKernel( 1 );

template< integer IS_THERMAL >
localIndex
PressureRelationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const iwelemControl,
          WellControls const & wellControls,
          real64 const & timeAtEndOfStep,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity,
          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  using Deriv = constitutive::singlefluid::DerivativeOffset;
  // static well control data
  bool const isProducer = wellControls.isProducer();
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const targetBHP = wellControls.getTargetBHP( timeAtEndOfStep );
  real64 const targetRate = wellControls.getTargetTotalRate( timeAtEndOfStep );

  // dynamic well control data
  real64 const & currentBHP =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentBHPString() );
  real64 const & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentBHP_dPresString() );

  real64 const & currentVolRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::currentVolRateString() );
  real64 const & dCurrentVolRate_dPres =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dPresString() );
  real64 const & dCurrentVolRate_dRate =
    wellControls.getReference< real64 >( SinglePhaseWell::viewKeyStruct::dCurrentVolRate_dRateString() );

  RAJA::ReduceMax< parallelDeviceReduce, localIndex > switchControl( 0 );

  // loop over the well elements to compute the pressure relations between well elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {

    localIndex const iwelemNext = nextWellElemIndex[iwelem];

    if( iwelemNext < 0 && isLocallyOwned ) // if iwelemNext < 0, form control equation
    {
      WellControls::Control newControl = currentControl;
      ControlEquationHelper::switchControl( isProducer,
                                            currentControl,
                                            targetBHP,
                                            targetRate,
                                            currentBHP,
                                            currentVolRate,
                                            newControl );
      if( currentControl != newControl )
      {
        switchControl.max( 1 );
      }

      ControlEquationHelper::compute< IS_THERMAL >( rankOffset,
                                                    newControl,
                                                    targetBHP,
                                                    targetRate,
                                                    currentBHP,
                                                    dCurrentBHP_dPres,
                                                    currentVolRate,
                                                    dCurrentVolRate_dPres,
                                                    dCurrentVolRate_dRate,
                                                    wellElemDofNumber[iwelemControl],
                                                    localMatrix,
                                                    localRhs );
    }
    else if( iwelemNext >= 0 )  // if iwelemNext >= 0, form momentum equation
    {

      // local working variables and arrays
      globalIndex dofColIndices[2]{};
      real64 localPresRelJacobian[2]{};

      // compute avg density
      real64 const avgDensity = 0.5 * ( wellElemDensity[iwelem][0] + wellElemDensity[iwelemNext][0] );
      real64 const dAvgDensity_dPresNext    = 0.5 * dWellElemDensity[iwelemNext][0][Deriv::dP];
      real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemDensity[iwelem][0][Deriv::dP];

      // compute depth diff times acceleration
      real64 const gravD = wellElemGravCoef[iwelemNext] - wellElemGravCoef[iwelem];

      // compute the current pressure in the two well elements
      real64 const pressureCurrent = wellElemPressure[iwelem];
      real64 const pressureNext    = wellElemPressure[iwelemNext];

      // compute momentum flux and derivatives
      real64 const localPresRel = pressureNext - pressureCurrent - avgDensity * gravD;
      localPresRelJacobian[TAG::NEXT]    =  1 - dAvgDensity_dPresNext * gravD;
      localPresRelJacobian[TAG::CURRENT] = -1 - dAvgDensity_dPresCurrent * gravD;

      // TODO: add friction and acceleration terms

      // jacobian indices
      globalIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::CONTROL - rankOffset;
      dofColIndices[TAG::NEXT]      = wellElemDofNumber[iwelemNext] + COFFSET::DPRES;
      dofColIndices[TAG::CURRENT]   = wellElemDofNumber[iwelem] + COFFSET::DPRES;

      if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                          &dofColIndices[0],
                                                                          &localPresRelJacobian[0],
                                                                          2 );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], localPresRel );
      }
    }
  } );
  return switchControl.get();
}

/******************************** PerforationKernel ********************************/

template< integer IS_THERMAL >
GEOS_HOST_DEVICE
inline
void
PerforationKernel::
  compute( real64 const & resPressure,
           real64 const & resDensity,
           arraySlice1d< real64 const > const & dResDensity,
           real64 const & dResDensity_dPres,
           real64 const & resViscosity,
           arraySlice1d< real64 const > const & dResViscosity,
           real64 const & dResViscosity_dPres,
           real64 const & wellElemGravCoef,
           real64 const & wellElemPressure,
           real64 const & wellElemDensity,
           arraySlice1d< real64 const > const & dWellElemDensity,
           real64 const & dWellElemDensity_dPres,
           real64 const & wellElemViscosity,
           arraySlice1d< real64 const > const & dWellElemViscosity,
           real64 const & dWellElemViscosity_dPres,
           real64 const & perfGravCoef,
           real64 const & trans,
           real64 & perfRate,
           arraySlice2d< real64 > const & dPerfRate,
           arraySlice1d< real64 > const & dPerfRate_dPres )
{

  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

  // local working variables and arrays
  real64 pressure[2]{};
  real64 dPressure[2][DerivOffset::nDer]{};

  real64 dPressure_dP[2]{};  // tjb delete
  real64 multiplier[2]{};

  perfRate = 0.0;
  for( integer i = 0; i < 2; ++i )
  {
    for( integer j=0; j<DerivOffset::nDer; j++ )
    {
      dPerfRate[i][j] = 0.0;
    }
  }

  dPerfRate_dPres[0] = 0.0;  // tjb delete
  dPerfRate_dPres[1] = 0.0;  // tjb delete

  // 1) Reservoir side
  // get reservoir variables
  pressure[TAG::RES] = resPressure;
  dPressure_dP[TAG::RES] = 1;  // tjb delete
  dPressure[TAG::RES][DerivOffset::dP] = 1.0;

  // TODO: add a buoyancy term for the reservoir side here

  // multiplier for reservoir side in the flux
  multiplier[TAG::RES] = 1;

  // 2) Well side


  // get well variables
  pressure[TAG::WELL] = wellElemPressure;
  dPressure[TAG::WELL][DerivOffset::dP] = 1.0;
  dPressure_dP[TAG::WELL] = 1.0; // tjb delete

  real64 const gravD = ( perfGravCoef - wellElemGravCoef );
  pressure[TAG::WELL]     += wellElemDensity * gravD;
  dPressure[TAG::WELL][DerivOffset::dP] += dWellElemDensity[DerivOffset::dP] * gravD;
  if constexpr (IS_THERMAL)
  {
    dPressure[TAG::WELL][DerivOffset::dT] = dWellElemDensity[DerivOffset::dT] * gravD;
  }
  dPressure_dP[TAG::WELL] += dWellElemDensity_dPres * gravD; // tjb remove

  // multiplier for well side in the flux
  multiplier[TAG::WELL] = -1;

  // compute potential difference
  real64 potDif = 0.0;
  for( localIndex i = 0; i < 2; ++i )
  {
    potDif += multiplier[i] * trans * pressure[i];
    dPerfRate_dPres[i] = multiplier[i] * trans * dPressure_dP[i];// tjb remove
    for( integer j=0; j<DerivOffset::nDer; j++ )
    {
      dPerfRate[i][j] = multiplier[i] * trans * dPressure[i][j];
    }
  }


  // choose upstream cell based on potential difference
  localIndex const k_up = (potDif >= 0) ? TAG::RES : TAG::WELL;

  // compute upstream density, viscosity, and mobility
  real64 densityUp       = 0.0;
  real64 viscosityUp     = 0.0;
  real64 dDensityUp[DerivOffset::nDer]{};
  real64 dViscosityUp[DerivOffset::nDer]{};
  for( integer j=0; j<DerivOffset::nDer; j++ )
  {
    dDensityUp[j]=0.0;
    dViscosityUp[j]=0.0;
  }
  real64 dDensityUp_dP   = 0.0;  // tjb remove
  real64 dViscosityUp_dP = 0.0; // tjb remove

  // upwinding the variables
  if( k_up == TAG::RES ) // use reservoir vars
  {
    densityUp     = resDensity;
    viscosityUp     = resViscosity;
    for( integer j=0; j<DerivOffset::nDer; j++ )
    {
      dDensityUp[j] = dResDensity[j];
      dViscosityUp[j] = dResViscosity[j];
    }
    dDensityUp_dP = dResDensity_dPres; // tjb remove
    dViscosityUp_dP = dResViscosity_dPres; // tjb remove
  }
  else // use well vars
  {
    densityUp = wellElemDensity;
    viscosityUp = wellElemViscosity;
    for( integer j=0; j<DerivOffset::nDer; j++ )
    {
      dDensityUp[j] = dWellElemDensity[j];
      dViscosityUp[j] = dWellElemViscosity[j];
    }
    dDensityUp_dP = dWellElemDensity_dPres;  // tjb remove
    dViscosityUp_dP = dWellElemViscosity_dPres;  // tjb remove
  }

  // compute mobility
  real64 const mobilityUp     = densityUp / viscosityUp;
  real64 dMobilityUp[DerivOffset::nDer]{};
  for( integer j=0; j<DerivOffset::nDer; j++ )
  {
    dMobilityUp[j] = dDensityUp[j] / viscosityUp
                     - mobilityUp / viscosityUp * dViscosityUp[j];
  }
  real64 const dMobilityUp_dP = dDensityUp_dP / viscosityUp
                                - mobilityUp / viscosityUp * dViscosityUp_dP;  // tjb remove

  perfRate = mobilityUp * potDif;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    for( integer j=0; j<DerivOffset::nDer; j++ )
    {
      dPerfRate[ke][j]  *= mobilityUp;
    }
    dPerfRate_dPres[ke] *= mobilityUp; // tjb remove
  }
  for( integer j=0; j<DerivOffset::nDer; j++ )
  {
    dPerfRate[k_up][j]  += dMobilityUp[j]* potDif;
  }
  dPerfRate_dPres[k_up] += dMobilityUp_dP * potDif; // tjb remove
  assert( fabs( dPerfRate_dPres[k_up] -dPerfRate[k_up][DerivOffset::dP] ) < FLT_EPSILON );
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    assert( fabs( dPerfRate_dPres[ke] -dPerfRate[ke][DerivOffset::dP] ) < FLT_EPSILON );
  }
}
#define INST_PerforationKernel( IS_THERMAL ) \
  template  \
  void \
  PerforationKernel:: \
    launch< IS_THERMAL >( localIndex const size, \
                          ElementViewConst< arrayView1d< real64 const > > const & resPressure, \
                          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & resDensity, \
                          ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const & dResDensity, \
                          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & resViscosity, \
                          ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const & dResViscosity, \
                          arrayView1d< real64 const > const & wellElemGravCoef, \
                          arrayView1d< real64 const > const & wellElemPressure, \
                          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity, \
                          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity, \
                          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemViscosity, \
                          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemViscosity, \
                          arrayView1d< real64 const > const & perfGravCoef, \
                          arrayView1d< localIndex const > const & perfWellElemIndex, \
                          arrayView1d< real64 const > const & perfTransmissibility, \
                          arrayView1d< localIndex const > const & resElementRegion, \
                          arrayView1d< localIndex const > const & resElementSubRegion, \
                          arrayView1d< localIndex const > const & resElementIndex, \
                          arrayView1d< real64 > const & perfRate, \
                          arrayView3d< real64 > const & dPerfRate, \
                          arrayView2d< real64 > const & dPerfRate_dPres )

INST_PerforationKernel( 0 );
INST_PerforationKernel( 1 );


template< integer IS_THERMAL >
void
PerforationKernel::
  launch( localIndex const size,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & resDensity,
          ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const & dResDensity,
          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & resViscosity,
          ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const & dResViscosity,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity,
          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemViscosity,
          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemViscosity,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< localIndex const > const & perfWellElemIndex,
          arrayView1d< real64 const > const & perfTransmissibility,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 > const & perfRate,
          arrayView3d< real64 > const & dPerfRate,
          arrayView2d< real64 > const & dPerfRate_dPres )
{
  using Deriv = constitutive::singlefluid::DerivativeOffset;
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
  {

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];

    // get the local index of the well element
    localIndex const iwelem = perfWellElemIndex[iperf];

    compute< IS_THERMAL >( resPressure[er][esr][ei],
                           resDensity[er][esr][ei][0],
                           dResDensity[er][esr][ei][0],
                           dResDensity[er][esr][ei][0][Deriv::dP],  // tjb remove
                           resViscosity[er][esr][ei][0],
                           dResViscosity[er][esr][ei][0],
                           dResViscosity[er][esr][ei][0][Deriv::dP],// tjb remove
                           wellElemGravCoef[iwelem],
                           wellElemPressure[iwelem],
                           wellElemDensity[iwelem][0],
                           dWellElemDensity[iwelem][0],
                           dWellElemDensity[iwelem][0][Deriv::dP],// tjb remove
                           wellElemViscosity[iwelem][0],
                           dWellElemViscosity[iwelem][0],
                           dWellElemViscosity[iwelem][0][Deriv::dP],// tjb remove
                           perfGravCoef[iperf],
                           perfTransmissibility[iperf],
                           perfRate[iperf],
                           dPerfRate[iperf],
                           dPerfRate_dPres[iperf] );


  } );
}


/******************************** AccumulationKernel ********************************/

void
AccumulationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity,
          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity_n,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  using Deriv = constitutive::singlefluid::DerivativeOffset;
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {

    if( wellElemGhostRank[iwelem] >= 0 )
    {
      return;
    }

    localIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - rankOffset;
    globalIndex const presDofColIndex = wellElemDofNumber[iwelem] + COFFSET::DPRES;

    real64 const localAccum = wellElemVolume[iwelem] * ( wellElemDensity[iwelem][0] - wellElemDensity_n[iwelem][0] );
    real64 const localAccumJacobian = wellElemVolume[iwelem] * dWellElemDensity[iwelem][0][Deriv::dP];

    // add contribution to global residual and jacobian (no need for atomics here)
    localMatrix.addToRow< serialAtomic >( eqnRowIndex, &presDofColIndex, &localAccumJacobian, 1 );
    localRhs[eqnRowIndex] += localAccum;
  } );
}


/******************************** PressureInitializationKernel ********************************/

void
PresTempInitializationKernel::
  launch( integer const isThermal,
          localIndex const perforationSize,
          localIndex const subRegionSize,
          localIndex const numPerforations,
          WellControls const & wellControls,
          real64 const & currentTime,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView1d< real64 const > > const & resTemp,
          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & resDensity,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPressure,
          arrayView1d< real64 > const & wellElemTemperature )
{
  real64 const targetBHP = wellControls.getTargetBHP( currentTime );
  real64 const refWellElemGravCoef = wellControls.getReferenceGravityCoef();
  real64 const initialPressureCoef = wellControls.getInitialPressureCoefficient();
  WellControls::Control const currentControl = wellControls.getControl();
  bool const isProducer = wellControls.isProducer();



  // Step 1: we loop over all the perforations on this rank to compute the following quantities:
  //   - Sum of densities over the perforated reservoir elements
  // In passing, we save the min gravCoef difference between the reference depth and the perforation depth
  // Note that we use gravCoef instead of depth for the (unlikely) case in which the gravityVector is not aligned with z

  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumDensity( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTemp( 0 );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > localMinGravCoefDiff( 1e9 );

  forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
  {
    // get the reservoir (sub)region and element indices
    localIndex const er = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei = resElementIndex[iperf];

    // save the min gravCoef difference between the reference depth and the perforation depth (times g)
    localMinGravCoefDiff.min( LvArray::math::abs( refWellElemGravCoef - perfGravCoef[iperf] ) );

    // increment the fluid density
    sumDensity += resDensity[er][esr][ei][0];

    // increment the temperature
    sumTemp += resTemp[er][esr][ei];
  } );

  real64 const minGravCoefDiff = MpiWrapper::min( localMinGravCoefDiff.get() );



  // Step 2: we assign average quantities over the well (i.e., over all the ranks)

  real64 const avgDensity = MpiWrapper::sum( sumDensity.get() ) / numPerforations;

  // Step 3: we compute the approximate pressure at the reference depth
  // We make a distinction between pressure-controlled wells and rate-controlled wells

  real64 refPres = 0.0;
  real64 avgTemp = 0;


  // for a producer, we use the temperature and component fractions from the reservoir
  if( isProducer )
  {
    // use average temperature from reservoir
    avgTemp = MpiWrapper::sum( sumTemp.get() ) / numPerforations;
  }
  // for an injector, we use the injection stream values
  else
  {
    // use temperature from injection stream
    avgTemp = wellControls.getInjectionTemperature();
  }
  // if the well is controlled by pressure, initialize the reference pressure at the target pressure
  if( currentControl == WellControls::Control::BHP )
  {
    refPres = targetBHP;
  }
  // if the well is controlled by rate, initialize the reference pressure using the pressure at the closest perforation
  else
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > localRefPres( 1e9 );
    real64 const alpha = ( isProducer ) ? 1 - initialPressureCoef : 1 + initialPressureCoef;

    forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      // get the perforation pressure and save the estimated reference pressure
      real64 const gravCoefDiff = LvArray::math::abs( refWellElemGravCoef - perfGravCoef[iperf] );
      if( isZero( gravCoefDiff - minGravCoefDiff ) )
      {
        localRefPres.min( alpha * resPressure[er][esr][ei] + avgDensity * ( refWellElemGravCoef - perfGravCoef[iperf] ) );
      }
    } );
    refPres = MpiWrapper::min( localRefPres.get() );
  }



  // Step 4: we are ready to assign the primary variables on the well elements:
  //  - pressure: hydrostatic pressure using our crude approximation of the total mass density

  RAJA::ReduceMax< parallelDeviceReduce, integer > foundNegativePressure( 0 );
  RAJA::ReduceMax< parallelDeviceReduce, integer > foundNegativeTemp( 0 );

  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    wellElemPressure[iwelem] = refPres + avgDensity * ( wellElemGravCoef[iwelem] - refWellElemGravCoef );
    wellElemTemperature[iwelem] = avgTemp;
    if( wellElemPressure[iwelem] <= 0 )
    {
      foundNegativePressure.max( 1 );
    }
    if( wellElemTemperature[iwelem] < 0 )
    {
      foundNegativeTemp.max( 1 );
    }
  } );

  GEOS_THROW_IF( foundNegativePressure.get() == 1,
                 wellControls.getDataContext() << ": Invalid well initialization, negative pressure was found.",
                 InputError );
  if( isThermal )   // tjb change  temp in isothermal cases shouldnt be an issue (also what if temp in fluid prop calcs like compo)
  {
    GEOS_THROW_IF( foundNegativeTemp.get() == 1,
                   wellControls.getDataContext() << "Invalid well initialization, negative temperature was found.",
                   InputError );
  }
}

/******************************** RateInitializationKernel ********************************/

void
RateInitializationKernel::
  launch( localIndex const subRegionSize,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDens,
          arrayView1d< real64 > const & connRate )
{
  real64 const targetRate = wellControls.getTargetTotalRate( currentTime );
  WellControls::Control const control = wellControls.getControl();
  bool const isProducer = wellControls.isProducer();

  // Estimate the connection rates
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    if( control == WellControls::Control::BHP )
    {
      // if BHP constraint set rate below the absolute max rate
      // with the appropriate sign (negative for prod, positive for inj)
      if( isProducer )
      {
        connRate[iwelem] = LvArray::math::max( 0.1 * targetRate * wellElemDens[iwelem][0], -1e3 );
      }
      else
      {
        connRate[iwelem] = LvArray::math::min( 0.1 * targetRate * wellElemDens[iwelem][0], 1e3 );
      }
    }
    else
    {
      connRate[iwelem] = targetRate * wellElemDens[iwelem][0];
    }
  } );
}

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */


template< integer IS_THERMAL >
class FaceBasedAssemblyKernel
{
public:

  using COFFSET = singlePhaseWellKernels::ColOffset;
  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using TAG = singlePhaseWellKernels::ElemTag;

  using FLUID_PROP_COFFSET = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;
  using WJ_COFFSET = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  using WJ_ROFFSET = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;

  using CP_Deriv = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;

  /// Number of Dof's set in this kernal
  static constexpr integer numDof = WJ_COFFSET::nDer;  // tjb revisit

  /// Compile time value for the number of equations except rate, momentum, energy
  static constexpr integer numEqn = WJ_COFFSET::nDer;// tjb revisit

  static constexpr integer maxNumElems = 2;
  static constexpr integer maxStencilSize = 2;
  /**
   * @brief Constructor for the kernel interface
   * @param[in] dt time step size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor
   * @param[in] wellControls well information
   * @param[in] subRegion  region containing well
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( real64 const dt,
                           globalIndex const rankOffset,
                           string const wellDofKey,
                           WellControls const & wellControls,
                           WellElementSubRegion const & subRegion,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    :
    m_dt( dt ),
    m_rankOffset( rankOffset ),
    m_wellElemDofNumber ( subRegion.getReference< array1d< globalIndex > >( wellDofKey ) ),
    m_nextWellElemIndex ( subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString()) ),
    m_connRate ( subRegion.getField< geos::fields::well::connectionRate >() ),
    m_localMatrix( localMatrix ),
    m_localRhs ( localRhs ),
    m_isProducer ( wellControls.isProducer() )
  {}


  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] ie the element index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iwelem,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {
    // 1) Compute the flux and its derivatives

    /*  currentConnRate < 0 flow from iwelem to iwelemNext
     *  currentConnRate > 0 flow from iwelemNext to iwelem
     *  With this convention, currentConnRate < 0 at the last connection for a producer
     *                        currentConnRate > 0 at the last connection for a injector
     */

    // get next well element index
    localIndex const iwelemNext = m_nextWellElemIndex[iwelem];

    // there is nothing to upwind for single-phase flow
    real64 const currentConnRate = m_connRate[iwelem];
    real64 const flux = m_dt * currentConnRate;
    real64 const dFlux_dRate = m_dt;

    // 2) Assemble the flux into residual and Jacobian
    if( iwelemNext < 0 )
    {
      // flux terms
      real64 const oneSidedLocalFlux = -flux;
      real64 const oneSidedLocalFluxJacobian_dRate = -dFlux_dRate;

      // jacobian indices
      globalIndex const oneSidedEqnRowIndex = m_wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - m_rankOffset;
      globalIndex const oneSidedDofColIndex_dRate = m_wellElemDofNumber[iwelem] + COFFSET::DRATE;

      if( oneSidedEqnRowIndex >= 0 && oneSidedEqnRowIndex < m_localMatrix.numRows() )
      {
        m_localMatrix.addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndex,
                                                        &oneSidedDofColIndex_dRate,
                                                        &oneSidedLocalFluxJacobian_dRate,
                                                        1 );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[oneSidedEqnRowIndex], oneSidedLocalFlux );
      }
    }
    else
    {
      // local working variables and arrays
      globalIndex eqnRowIndices[2]{};

      real64 localFlux[2]{};
      real64 localFluxJacobian_dRate[2]{};

      // flux terms
      localFlux[TAG::NEXT]    =   flux;
      localFlux[TAG::CURRENT] = -flux;

      localFluxJacobian_dRate[TAG::NEXT]    =   dFlux_dRate;
      localFluxJacobian_dRate[TAG::CURRENT] = -dFlux_dRate;

      // indices
      eqnRowIndices[TAG::CURRENT] = m_wellElemDofNumber[iwelem] + ROFFSET::MASSBAL - m_rankOffset;
      eqnRowIndices[TAG::NEXT]    = m_wellElemDofNumber[iwelemNext] + ROFFSET::MASSBAL - m_rankOffset;
      globalIndex const dofColIndex_dRate = m_wellElemDofNumber[iwelem] + COFFSET::DRATE;

      for( localIndex i = 0; i < 2; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < m_localMatrix.numRows() )
        {
          m_localMatrix.addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                          &dofColIndex_dRate,
                                                          &localFluxJacobian_dRate[i],
                                                          1 );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[eqnRowIndices[i]], localFlux[i] );
        }
      }
    }


  }


  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElements the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElements,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const ie )
    {
      kernelComponent.computeFlux( ie );
    } );
  }

protected:
  /// Time step size
  real64 const m_dt;
  /// Rank offset for calculating row/col Jacobian indices
  integer const m_rankOffset;

  /// Reference to the degree-of-freedom numbers
  arrayView1d< globalIndex const > const m_wellElemDofNumber;
  /// Next element index, needed since iterating over element nodes, not edges
  arrayView1d< localIndex const > const m_nextWellElemIndex;

  /// Connection rate
  arrayView1d< real64 const > const m_connRate;

  /// Element component fraction
  arrayView2d< real64 const, compflow::USD_COMP > const m_wellElemCompFrac;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  /// Well type
  bool const m_isProducer;

};

/**
 * @class FaceBasedAssemblyKernelFactory
 */
class FaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] dt time step size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] wellControls object holding well control/constraint information
   * @param[in] subregion well subregion
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComps,
                   real64 const dt,
                   globalIndex const rankOffset,
                   string const dofKey,
                   WellControls const & wellControls,
                   WellElementSubRegion const & subRegion,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer isThermal=0;
    geos::internal::kernelLaunchSelectorThermalSwitch( isThermal, [&]( auto IS_THERMAL )
    {

      integer constexpr istherm = IS_THERMAL();

      using kernelType = FaceBasedAssemblyKernel< istherm >;
      kernelType kernel( dt, rankOffset, dofKey, wellControls, subRegion, localMatrix, localRhs );
      kernelType::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }
};

} // end namespace singlePhaseWellKernels

} // end namespace geos
