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
 * @file CompositionalMultiphaseWellKernels.cpp
 */

#include "CompositionalMultiphaseWellKernels.hpp"

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
// TODO: move keys to WellControls
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"

namespace geos
{

using namespace constitutive;
namespace compositionalMultiphaseWellKernels
{

/******************************** ControlEquationHelper ********************************/

GEOS_HOST_DEVICE
inline
void
ControlEquationHelper::
  switchControl( bool const isProducer,
                 WellControls::Control const & currentControl,
                 integer const phasePhaseIndex,
                 real64 const & targetBHP,
                 real64 const & targetPhaseRate,
                 real64 const & targetTotalRate,
                 real64 const & targetMassRate,
                 real64 const & currentBHP,
                 arrayView1d< real64 const > const & currentPhaseVolRate,
                 real64 const & currentTotalVolRate,
                 WellControls::Control & newControl )
{
  // if isViable is true at the end of the following checks, no need to switch
  bool controlIsViable = false;

  // The limiting flow rates are treated as upper limits, while the pressure limits
  // are treated as lower limits in production wells and upper limits in injectors.
  // The well changes its mode of control whenever the existing control mode would
  // violate one of these limits.

  // Currently, the available constraints are:
  //   - Producer: BHP, PHASEVOLRATE
  //   - Injector: BHP, TOTALVOLRATE

  // TODO: support GRAT, WRAT, LIQUID for producers and check if any of the active constraint is violated

  // BHP control
  if( currentControl == WellControls::Control::BHP )
  {
    // the control is viable if the reference oil rate is below the max rate for producers
    if( isProducer )
    {
      controlIsViable = ( LvArray::math::abs( currentPhaseVolRate[phasePhaseIndex] ) <= LvArray::math::abs( targetPhaseRate ) );
    }
    // the control is viable if the reference total rate is below the max rate for injectors
    else
    {
      controlIsViable = ( LvArray::math::abs( currentTotalVolRate ) <= LvArray::math::abs( targetTotalRate ) );
    }
  }
  else // rate control
  {
    // the control is viable if the reference pressure is below/above the max/min pressure
    if( isProducer )
    {
      // targetBHP specifies a min pressure here
      controlIsViable = ( currentBHP >= targetBHP );
    }
    else
    {
      // targetBHP specifies a max pressure here
      controlIsViable = ( currentBHP <= targetBHP );
    }
  }

  if( controlIsViable )
  {
    newControl = currentControl;
  }
  else
  {
    if( isProducer )
    {
      newControl = ( currentControl == WellControls::Control::BHP )
           ? WellControls::Control::PHASEVOLRATE
           : WellControls::Control::BHP;
    }
    else
    {
      if( isZero( targetMassRate ) )
      {
        newControl = ( currentControl == WellControls::Control::BHP )
                  ? WellControls::Control::TOTALVOLRATE
                  : WellControls::Control::BHP;
      }
      else
      {
        newControl = ( currentControl == WellControls::Control::BHP )
                  ? WellControls::Control::MASSRATE
                  : WellControls::Control::BHP;
      }
    }
  }
}

template< integer NC, integer IS_THERMAL >
GEOS_HOST_DEVICE
inline
void
ControlEquationHelper::
  compute( globalIndex const rankOffset,
           WellControls::Control const currentControl,
           integer const targetPhaseIndex,
           real64 const & targetBHP,
           real64 const & targetPhaseRate,
           real64 const & targetTotalRate,
           real64 const & targetMassRate,
           real64 const & currentBHP,
           arrayView1d< real64 const > const & dCurrentBHP,
           arrayView1d< real64 const > const & currentPhaseVolRate,
           arrayView2d< real64 const > const & dCurrentPhaseVolRate,

           real64 const & currentTotalVolRate,
           arrayView1d< real64 const > const & dCurrentTotalVolRate,
           real64 const & massDensity,
           globalIndex const dofNumber,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{

  using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
  using Deriv = multifluid::DerivativeOffset;

  localIndex const eqnRowIndex      = dofNumber + ROFFSET::CONTROL - rankOffset;
  globalIndex dofColIndices[COFFSET_WJ::nDer]{};
  for( integer ic = 0; ic < COFFSET_WJ::nDer; ++ic )
  {
    dofColIndices[ ic ] = dofNumber + ic;
  }

  real64 controlEqn = 0;
  real64 dControlEqn[NC+2+IS_THERMAL]{};

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
    dControlEqn[COFFSET_WJ::dP] = dCurrentBHP[Deriv::dP];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn[COFFSET_WJ::dC+ic] = dCurrentBHP[Deriv::dC+ic];
    }
    if constexpr ( IS_THERMAL )

      dControlEqn[COFFSET_WJ::dT] = dCurrentBHP[Deriv::dT];

  }
  // Oil volumetric rate control
  else if( currentControl == WellControls::Control::PHASEVOLRATE )
  {
    controlEqn = currentPhaseVolRate[targetPhaseIndex] - targetPhaseRate;
    dControlEqn[COFFSET_WJ::dP] = dCurrentPhaseVolRate[targetPhaseIndex][COFFSET_WJ::dP];
    dControlEqn[COFFSET_WJ::dQ] = dCurrentPhaseVolRate[targetPhaseIndex][COFFSET_WJ::dQ];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn[COFFSET_WJ::dC+ic] = dCurrentPhaseVolRate[targetPhaseIndex][COFFSET_WJ::dC+ic];
    }
    if constexpr ( IS_THERMAL )
      dControlEqn[COFFSET_WJ::dT] = dCurrentBHP[Deriv::dT];
  }
  // Total volumetric rate control
  else if( currentControl == WellControls::Control::TOTALVOLRATE )
  {
    controlEqn = currentTotalVolRate - targetTotalRate;
    dControlEqn[COFFSET_WJ::dP] = dCurrentTotalVolRate[COFFSET_WJ::dP];
    dControlEqn[COFFSET_WJ::dQ] = dCurrentTotalVolRate[COFFSET_WJ::dQ];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn[COFFSET_WJ::dC+ic] = dCurrentTotalVolRate[COFFSET_WJ::dC+ic];
    }
    if constexpr ( IS_THERMAL )
      dControlEqn[COFFSET_WJ::dT] = dCurrentTotalVolRate[COFFSET_WJ::dT];
  }
  // Total mass rate control
  else if( currentControl == WellControls::Control::MASSRATE )
  {
    controlEqn = massDensity*currentTotalVolRate - targetMassRate;
    dControlEqn[COFFSET_WJ::dP] = massDensity*dCurrentTotalVolRate[COFFSET_WJ::dP];
    dControlEqn[COFFSET_WJ::dQ] = massDensity*dCurrentTotalVolRate[COFFSET_WJ::dQ];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn[COFFSET_WJ::dC+ic] = massDensity*dCurrentTotalVolRate[COFFSET_WJ::dC+ic];
    }
    if constexpr ( IS_THERMAL )
      dControlEqn[COFFSET_WJ::dT] = massDensity*dCurrentTotalVolRate[COFFSET_WJ::dT];
  }
  // Total mass rate control
  else if( currentControl == WellControls::Control::MASSRATE )
  {
    controlEqn = massDensity*currentTotalVolRate - targetMassRate;
    dControlEqn[COFFSET_WJ::dP] = massDensity*dCurrentTotalVolRate[COFFSET_WJ::dP];
    dControlEqn[COFFSET_WJ::dQ] = massDensity*dCurrentTotalVolRate[COFFSET_WJ::dQ];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn[COFFSET_WJ::dC+ic] = massDensity*dCurrentTotalVolRate[COFFSET_WJ::dC+ic];
    }
  }
  else
  {
    GEOS_ERROR( "This constraint is not supported in CompositionalMultiphaseWell" );
  }
  localRhs[eqnRowIndex] += controlEqn;

  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                            dofColIndices,
                                                            dControlEqn,
                                                            COFFSET_WJ::nDer );


}

/******************************** PressureRelationKernel ********************************/

template< integer NC, integer IS_THERMAL >
GEOS_HOST_DEVICE
void
PressureRelationKernel::
  compute( real64 const & gravCoef,
           real64 const & gravCoefNext,
           real64 const & pres,
           real64 const & presNext,
           real64 const & totalMassDens,
           real64 const & totalMassDensNext,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDensNext,
           real64 & localPresRel,
           real64 ( & localPresRelJacobian )[2*(NC+1 + IS_THERMAL)] )
{
  // local working variables and arrays
  real64 dAvgMassDens_dCompCurrent[NC]{};
  real64 dAvgMassDens_dCompNext[NC]{};

  // compute the average density at the interface between well elements
  real64 const avgMassDens = 0.5 * ( totalMassDensNext + totalMassDens );
  real64 const dAvgMassDens_dPresNext    = 0.5 * dTotalMassDensNext[Deriv::dP];
  real64 const dAvgMassDens_dPresCurrent = 0.5 * dTotalMassDens[Deriv::dP];
  for( integer ic = 0; ic < NC; ++ic )
  {
    dAvgMassDens_dCompNext[ic]    = 0.5 * dTotalMassDensNext[Deriv::dC+ic];
    dAvgMassDens_dCompCurrent[ic] = 0.5 * dTotalMassDens[Deriv::dC+ic];
  }

  // compute depth diff times acceleration
  real64 const gravD = gravCoefNext - gravCoef;

  // TODO: add friction and acceleration terms

  localPresRel = ( presNext - pres - avgMassDens * gravD );

  // localPresRelJacbain contains dP, dC and potentially dT derivatives for neighboring well elements
  // TAG::NEXT is 1, CURRENT is 0 , not sure why indexes are setup as below
  localPresRelJacobian[TAG::NEXT *(NC+1+IS_THERMAL)]    = ( 1 - dAvgMassDens_dPresNext * gravD );
  localPresRelJacobian[TAG::CURRENT *(NC+1+IS_THERMAL)] = ( -1 - dAvgMassDens_dPresCurrent * gravD );

  for( integer ic = 0; ic < NC; ++ic )
  {
    localPresRelJacobian[TAG::NEXT *(NC+1+IS_THERMAL) + ic+1]    = -dAvgMassDens_dCompNext[ic] * gravD;
    localPresRelJacobian[TAG::CURRENT *(NC+1+IS_THERMAL) + ic+1] = -dAvgMassDens_dCompCurrent[ic] * gravD;
  }
  if constexpr ( IS_THERMAL )
  {
    localPresRelJacobian[TAG::NEXT *(NC+1+IS_THERMAL)+NC+1]    =  0.5 * dTotalMassDensNext[Deriv::dT];
    localPresRelJacobian[TAG::CURRENT *(NC+1+IS_THERMAL)+NC+1] = 0.5 * dTotalMassDens[Deriv::dT];
  }
}

template< integer NC, integer IS_THERMAL >
void
PressureRelationKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const iwelemControl,
          integer const targetPhaseIndex,
          WellControls const & wellControls,
          real64 const & timeAtEndOfStep,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & wellElemTotalMassDens,
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens,
          bool & controlHasSwitched,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  using COFFSET_WJ = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
  // static well control data
  bool const isProducer = wellControls.isProducer();
  WellControls::Control const currentControl = wellControls.getControl();
  real64 const targetBHP = wellControls.getTargetBHP( timeAtEndOfStep );
  real64 const targetTotalRate = wellControls.getTargetTotalRate( timeAtEndOfStep );
  real64 const targetPhaseRate = wellControls.getTargetPhaseRate( timeAtEndOfStep );
  real64 const targetMassRate = wellControls.getTargetMassRate( timeAtEndOfStep );

  // dynamic well control data
  real64 const & currentBHP =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  arrayView1d< real64 const > const & dCurrentBHP =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHPString() );

  arrayView1d< real64 const > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  arrayView2d< real64 const > const & dCurrentPhaseVolRate =
    wellControls.getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRateString() );

  real64 const & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
  arrayView1d< real64 const > const & dCurrentTotalVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRateString() );

  real64 const & massDensity  =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::massDensityString() );

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
                                            targetPhaseIndex,
                                            targetBHP,
                                            targetPhaseRate,
                                            targetTotalRate,
                                            targetMassRate,
                                            currentBHP,
                                            currentPhaseVolRate,
                                            currentTotalVolRate,
                                            newControl );
      if( currentControl != newControl )
      {
        switchControl.max( 1 );
      }
      ControlEquationHelper::compute< NC, IS_THERMAL >( rankOffset,
                                                        newControl,
                                                        targetPhaseIndex,
                                                        targetBHP,
                                                        targetPhaseRate,
                                                        targetTotalRate,
                                                        targetMassRate,
                                                        currentBHP,
                                                        dCurrentBHP,
                                                        currentPhaseVolRate,
                                                        dCurrentPhaseVolRate,
                                                        currentTotalVolRate,
                                                        dCurrentTotalVolRate,
                                                        massDensity,
                                                        wellElemDofNumber[iwelemControl],
                                                        localMatrix,
                                                        localRhs );
      // TODO: for consistency, we should assemble here, not in compute...

    }
    else if( iwelemNext >= 0 ) // if iwelemNext >= 0, form momentum equation
    {

      real64 localPresRel = 0;
      real64 localPresRelJacobian[2*(NC+1+IS_THERMAL)]{};

      compute< NC, IS_THERMAL >(
        wellElemGravCoef[iwelem],
        wellElemGravCoef[iwelemNext],
        wellElemPressure[iwelem],
        wellElemPressure[iwelemNext],
        wellElemTotalMassDens[iwelem],
        wellElemTotalMassDens[iwelemNext],
        dWellElemTotalMassDens[iwelem],
        dWellElemTotalMassDens[iwelemNext],
        localPresRel,
        localPresRelJacobian );


      // local working variables and arrays
      globalIndex dofColIndices[2*(NC+1+IS_THERMAL)];

      globalIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::CONTROL - rankOffset;
      dofColIndices[TAG::NEXT *(NC+1+IS_THERMAL)]    = wellElemDofNumber[iwelemNext] + COFFSET_WJ::dP;
      dofColIndices[TAG::CURRENT *(NC+1+IS_THERMAL)] = wellElemDofNumber[iwelem] + COFFSET_WJ::dP;

      for( integer ic = 0; ic < NC; ++ic )
      {
        dofColIndices[TAG::NEXT *(NC+1+IS_THERMAL) + ic+1]    = wellElemDofNumber[iwelemNext] + COFFSET_WJ::dC + ic;
        dofColIndices[TAG::CURRENT *(NC+1+IS_THERMAL) + ic+1] = wellElemDofNumber[iwelem] + COFFSET_WJ::dC + ic;
      }
      if constexpr ( IS_THERMAL )
      {
        dofColIndices[TAG::NEXT *(NC+1+IS_THERMAL)+NC+1]    = wellElemDofNumber[iwelemNext] + COFFSET_WJ::dT;
        dofColIndices[TAG::CURRENT *(NC+1+IS_THERMAL)+NC+1] = wellElemDofNumber[iwelem] + COFFSET_WJ::dT;
      }
      if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                          dofColIndices,
                                                                          localPresRelJacobian,
                                                                          2 * (NC+1+IS_THERMAL) );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], localPresRel );
      }
    }
  } );
  controlHasSwitched = ( switchControl.get() == 1 );
}

#define INST_PressureRelationKernel( NC, IS_THERMAL ) \
  template \
  void PressureRelationKernel:: \
    launch< NC, IS_THERMAL >( localIndex const size, \
                              globalIndex const rankOffset, \
                              bool const isLocallyOwned, \
                              localIndex const iwelemControl, \
                              integer const targetPhaseIndex, \
                              WellControls const & wellControls, \
                              real64 const & timeAtEndOfStep, \
                              arrayView1d< globalIndex const > const & wellElemDofNumber, \
                              arrayView1d< real64 const > const & wellElemGravCoef, \
                              arrayView1d< localIndex const > const & nextWellElemIndex, \
                              arrayView1d< real64 const > const & wellElemPressure, \
                              arrayView1d< real64 const > const & wellElemTotalMassDens, \
                              arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens, \
                              bool & controlHasSwitched, \
                              CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                              arrayView1d< real64 > const & localRhs )

INST_PressureRelationKernel( 1, 0 );
INST_PressureRelationKernel( 1, 1 );
INST_PressureRelationKernel( 2, 0 );
INST_PressureRelationKernel( 2, 1 );
INST_PressureRelationKernel( 3, 0 );
INST_PressureRelationKernel( 3, 1 );
INST_PressureRelationKernel( 4, 0 );
INST_PressureRelationKernel( 4, 1 );
INST_PressureRelationKernel( 5, 0 );
INST_PressureRelationKernel( 5, 1 );

void
PresTempCompFracInitializationKernel::
  launch( localIndex const perforationSize,
          localIndex const subRegionSize,
          integer const numComps,
          integer const numPhases,
          localIndex const numPerforations,
          WellControls const & wellControls,
          real64 const & currentTime,
          ElementViewConst< arrayView1d< real64 const > > const & resPres,
          ElementViewConst< arrayView1d< real64 const > > const & resTemp,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const & resCompDens,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseMassDens,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPres,
          arrayView1d< real64 > const & wellElemTemp,
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac )
{
  integer constexpr MAX_NUM_COMP = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  real64 const targetBHP = wellControls.getTargetBHP( currentTime );
  real64 const refWellElemGravCoef = wellControls.getReferenceGravityCoef();
  real64 const initialPresCoef = wellControls.getInitialPressureCoefficient();
  WellControls::Control const currentControl = wellControls.getControl();
  bool const isProducer = wellControls.isProducer();



  // Step 1: we loop over all the perforations on this rank to compute the following quantities:
  //   - Sum of total mass densities over the perforated reservoir elements
  //   - Sum of the temperatures over the perforated reservoir elements
  //   - Sum of the component fractions over the perforated reservoir elements
  // In passing, we save the min gravCoef difference between the reference depth and the perforation depth
  // Note that we use gravCoef instead of depth for the (unlikely) case in which the gravityVector is not aligned with z

  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTotalMassDens( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumTemp( 0 );
  RAJA::ReduceSum< parallelDeviceReduce, real64 > sumCompFrac[MAX_NUM_COMP]{};
  RAJA::ReduceMin< parallelDeviceReduce, real64 > localMinGravCoefDiff( 1e9 );

  forAll< parallelDevicePolicy<> >( perforationSize, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
  {
    // get the reservoir (sub)region and element indices
    localIndex const er = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei = resElementIndex[iperf];

    // save the min gravCoef difference between the reference depth and the perforation depth (times g)
    localMinGravCoefDiff.min( LvArray::math::abs( refWellElemGravCoef - perfGravCoef[iperf] ) );

    // increment the temperature
    sumTemp += resTemp[er][esr][ei];

    // increment the total mass density
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      sumTotalMassDens += resPhaseVolFrac[er][esr][ei][ip] * resPhaseMassDens[er][esr][ei][0][ip];
    }

    // increment the component fractions
    real64 perfTotalDens = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      perfTotalDens += resCompDens[er][esr][ei][ic];
    }
    for( integer ic = 0; ic < numComps; ++ic )
    {
      sumCompFrac[ic] += resCompDens[er][esr][ei][ic] / perfTotalDens;
    }
  } );
  real64 const minGravCoefDiff = MpiWrapper::min( localMinGravCoefDiff.get() );



  // Step 2: we assign average quantities over the well (i.e., over all the ranks)
  // For composition and temperature, we make a distinction between injection and production

  // for total mass density, we always use the values of the perforated reservoir elements, even for injectors
  real64 const avgTotalMassDens = MpiWrapper::sum( sumTotalMassDens.get() ) / numPerforations;

//  stackArray1d< real64, MAX_NUM_COMP > avgCompFrac( numComps );
  array1d< real64 > aveCompFracArray( numComps );
  arrayView1d< real64 > const avgCompFrac = aveCompFracArray.toView();

  real64 avgTemp = 0;

  // for a producer, we use the temperature and component fractions from the reservoir
  if( isProducer )
  {
    // use average temperature from reservoir
    avgTemp = MpiWrapper::sum( sumTemp.get() ) / numPerforations;

    // use average comp frac from reservoir
    for( integer ic = 0; ic < numComps; ++ic )
    {
      avgCompFrac[ic] = MpiWrapper::sum( sumCompFrac[ic].get() ) / numPerforations;
    }
  }
  // for an injector, we use the injection stream values
  else
  {
    // use temperature from injection stream
    avgTemp = wellControls.getInjectionTemperature();

    // use comp frac from injection stream
    for( integer ic = 0; ic < numComps; ++ic )
    {
      avgCompFrac[ic] = wellControls.getInjectionStream()[ic];
    }
  }



  // Step 3: we compute the approximate pressure at the reference depth
  // We make a distinction between pressure-controlled wells and rate-controlled wells

  real64 refPres = 0.0;

  // if the well is controlled by pressure, initialize the reference pressure at the target pressure
  if( currentControl == WellControls::Control::BHP )
  {
    refPres = targetBHP;
  }
  // if the well is controlled by rate, initialize the reference pressure using the pressure at the closest perforation
  else
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > localRefPres( 1e9 );
    real64 const alpha = ( isProducer ) ? 1 - initialPresCoef : 1 + initialPresCoef;

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
        localRefPres.min( alpha * resPres[er][esr][ei] + avgTotalMassDens * ( refWellElemGravCoef - perfGravCoef[iperf] ) );
      }
    } );
    refPres = MpiWrapper::min( localRefPres.get() );
  }



  // Step 4: we are ready to assign the primary variables on the well elements:
  //  - pressure: hydrostatic pressure using our crude approximation of the total mass density
  //  - temperature: uniform, using the average temperature computed above
  //  - component fraction: uniform, using the average component fraction computed above

  RAJA::ReduceMax< parallelDeviceReduce, integer > foundNegativeTemp( 0 );
  RAJA::ReduceMax< parallelDeviceReduce, integer > foundNegativePres( 0 );
  RAJA::ReduceMax< parallelDeviceReduce, integer > foundInconsistentCompFrac( 0 );


  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    wellElemPres[iwelem] = refPres + avgTotalMassDens * ( wellElemGravCoef[iwelem] - refWellElemGravCoef );
    wellElemTemp[iwelem] = avgTemp;

    real64 sumCompFracForCheck = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      wellElemCompFrac[iwelem][ic] = avgCompFrac[ic];
      sumCompFracForCheck += wellElemCompFrac[iwelem][ic];
    }

    if( wellElemPres[iwelem] <= 0 )
    {
      foundNegativePres.max( 1 );
    }
    if( wellElemTemp[iwelem] <= 0 )
    {
      foundNegativeTemp.max( 1 );
    }
    if( !isZero( sumCompFracForCheck - 1.0, constitutive::MultiFluidConstants::minForSpeciesPresence ) )
    {
      foundInconsistentCompFrac.max( 1 );
    }

  } );


  GEOS_THROW_IF( foundNegativePres.get() == 1,
                 wellControls.getDataContext() << "Invalid well initialization, negative pressure was found.",
                 InputError );
  GEOS_THROW_IF( foundNegativeTemp.get() == 1,
                 wellControls.getDataContext() << "Invalid well initialization, negative temperature was found.",
                 InputError );
  GEOS_THROW_IF( foundInconsistentCompFrac.get() == 1,
                 wellControls.getDataContext() << "Invalid well initialization, inconsistent component fractions were found.",
                 InputError );


}

/******************************** CompDensInitializationKernel ********************************/

void
CompDensInitializationKernel::
  launch( localIndex const subRegionSize,
          integer const numComponents,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & wellElemTotalDens,
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens )
{
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    for( integer ic = 0; ic < numComponents; ++ic )
    {
      wellElemCompDens[iwelem][ic] = wellElemCompFrac[iwelem][ic] * wellElemTotalDens[iwelem][0];
    }
  } );
}

/******************************** RateInitializationKernel ********************************/

void
RateInitializationKernel::
  launch( localIndex const subRegionSize,
          integer const targetPhaseIndex,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens,
          arrayView1d< real64 > const & connRate )
{
  WellControls::Control const control = wellControls.getControl();
  bool const isProducer = wellControls.isProducer();
  real64 const targetTotalRate = wellControls.getTargetTotalRate( currentTime );
  real64 const targetPhaseRate = wellControls.getTargetPhaseRate( currentTime );
  real64 const targetMassRate = wellControls.getTargetMassRate( currentTime );

  // Estimate the connection rates
  forAll< parallelDevicePolicy<> >( subRegionSize, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    if( control == WellControls::Control::BHP )
    {
      // if BHP constraint set rate below the absolute max rate
      // with the appropriate sign (negative for prod, positive for inj)
      if( isProducer )
      {
        connRate[iwelem] = LvArray::math::max( 0.1 * targetPhaseRate * phaseDens[iwelem][0][targetPhaseIndex], -1e3 );
      }
      else
      {
        if( isZero( targetMassRate ) )
        {
          connRate[iwelem] = LvArray::math::min( 0.1 * targetTotalRate * totalDens[iwelem][0], 1e3 );
        }
        else
        {
          connRate[iwelem] = targetMassRate;
        }

      }
    }
    else if( control == WellControls::Control::MASSRATE )
    {
      connRate[iwelem] = targetMassRate;
      connRate[iwelem] = targetMassRate* totalDens[iwelem][0];
    }
    else
    {
      if( isProducer )
      {
        connRate[iwelem] = targetPhaseRate * phaseDens[iwelem][0][targetPhaseIndex];
      }
      else
      {
        connRate[iwelem] = targetTotalRate * totalDens[iwelem][0];
      }
    }
  } );
}


} // end namespace compositionalMultiphaseWellKernels

} // end namespace geos
