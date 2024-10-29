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
 * @file CompositionalMultiphaseWellKernels.cpp
 */

#include "CompositionalMultiphaseWellKernels.hpp"

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
// TODO: move keys to WellControls
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWell.hpp"

namespace geos
{

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

template< integer NC >
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
           real64 const & dCurrentBHP_dPres,
           arrayView1d< real64 const > const & dCurrentBHP_dCompDens,
           arrayView1d< real64 const > const & currentPhaseVolRate,
           arrayView1d< real64 const > const & dCurrentPhaseVolRate_dPres,
           arrayView2d< real64 const > const & dCurrentPhaseVolRate_dCompDens,
           arrayView1d< real64 const > const & dCurrentPhaseVolRate_dRate,
           real64 const & currentTotalVolRate,
           real64 const & dCurrentTotalVolRate_dPres,
           arrayView1d< real64 const > const & dCurrentTotalVolRate_dCompDens,
           real64 const & dCurrentTotalVolRate_dRate,
           real64 const & massDensity,
           globalIndex const dofNumber,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{
  localIndex const eqnRowIndex      = dofNumber + ROFFSET::CONTROL - rankOffset;
  globalIndex const presDofColIndex = dofNumber + COFFSET::DPRES;
  globalIndex const rateDofColIndex = dofNumber + COFFSET::DCOMP + NC;

  globalIndex compDofColIndices[NC]{};
  for( integer ic = 0; ic < NC; ++ic )
  {
    compDofColIndices[ ic ] = presDofColIndex + ic + 1;
  }

  real64 controlEqn = 0;
  real64 dControlEqn_dPres = 0;
  real64 dControlEqn_dRate = 0;
  real64 dControlEqn_dComp[NC]{};

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
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn_dComp[ic] = dCurrentBHP_dCompDens[ic];
    }
  }
  // Oil volumetric rate control
  else if( currentControl == WellControls::Control::PHASEVOLRATE )
  {
    controlEqn = currentPhaseVolRate[targetPhaseIndex] - targetPhaseRate;
    dControlEqn_dPres = dCurrentPhaseVolRate_dPres[targetPhaseIndex];
    dControlEqn_dRate = dCurrentPhaseVolRate_dRate[targetPhaseIndex];
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn_dComp[ic] = dCurrentPhaseVolRate_dCompDens[targetPhaseIndex][ic];
    }
  }
  // Total volumetric rate control
  else if( currentControl == WellControls::Control::TOTALVOLRATE )
  {
    controlEqn = currentTotalVolRate - targetTotalRate;
    dControlEqn_dPres = dCurrentTotalVolRate_dPres;
    dControlEqn_dRate = dCurrentTotalVolRate_dRate;
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn_dComp[ic] = dCurrentTotalVolRate_dCompDens[ic];
    }
  }
  // Total mass rate control
  else if( currentControl == WellControls::Control::MASSRATE )
  {
    controlEqn = massDensity*currentTotalVolRate - targetMassRate;
    dControlEqn_dPres = massDensity*dCurrentTotalVolRate_dPres;
    dControlEqn_dRate = massDensity*dCurrentTotalVolRate_dRate;
    for( integer ic = 0; ic < NC; ++ic )
    {
      dControlEqn_dComp[ic] = massDensity*dCurrentTotalVolRate_dCompDens[ic];
    }
  }
  else
  {
    GEOS_ERROR( "This constraint is not supported in CompositionalMultiphaseWell" );
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
  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                            compDofColIndices,
                                                            dControlEqn_dComp,
                                                            NC );
}

/******************************** FluxKernel ********************************/

template< integer NC >
GEOS_HOST_DEVICE
void
FluxKernel::
  computeExit( real64 const & dt,
               real64 const ( &compFlux )[NC],
               real64 const ( &dCompFlux_dRate )[NC],
               real64 const ( &dCompFlux_dPresUp )[NC],
               real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
               real64 ( & oneSidedFlux )[NC],
               real64 ( & oneSidedFluxJacobian_dRate )[NC][1],
               real64 ( & oneSidedFluxJacobian_dPresCompUp )[NC][NC + 1] )
{
  for( integer ic = 0; ic < NC; ++ic )
  {
    oneSidedFlux[ic] = -dt * compFlux[ic];

    // derivative with respect to rate
    oneSidedFluxJacobian_dRate[ic][0] = -dt * dCompFlux_dRate[ic];

    // derivative with respect to upstream pressure
    oneSidedFluxJacobian_dPresCompUp[ic][0] = -dt * dCompFlux_dPresUp[ic];

    // derivatives with respect to upstream component densities
    for( integer jdof = 0; jdof < NC; ++jdof )
    {
      oneSidedFluxJacobian_dPresCompUp[ic][jdof+1] = -dt * dCompFlux_dCompDensUp[ic][jdof];
    }
  }
}

template< integer NC >
GEOS_HOST_DEVICE
void
FluxKernel::
  compute( real64 const & dt,
           real64 const ( &compFlux )[NC],
           real64 const ( &dCompFlux_dRate )[NC],
           real64 const ( &dCompFlux_dPresUp )[NC],
           real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
           real64 ( & localFlux )[2*NC],
           real64 ( & localFluxJacobian_dRate )[2*NC][1],
           real64 ( & localFluxJacobian_dPresCompUp )[2*NC][NC + 1] )
{
  // flux terms
  for( integer ic = 0; ic < NC; ++ic )
  {
    localFlux[TAG::NEXT *NC+ic]    = dt * compFlux[ic];
    localFlux[TAG::CURRENT *NC+ic] = -dt * compFlux[ic];

    // derivative with respect to rate
    localFluxJacobian_dRate[TAG::NEXT *NC+ic][0]    = dt * dCompFlux_dRate[ic];
    localFluxJacobian_dRate[TAG::CURRENT *NC+ic][0] = -dt * dCompFlux_dRate[ic];

    // derivative with respect to upstream pressure
    localFluxJacobian_dPresCompUp[TAG::NEXT *NC+ic][0]    = dt * dCompFlux_dPresUp[ic];
    localFluxJacobian_dPresCompUp[TAG::CURRENT *NC+ic][0] = -dt * dCompFlux_dPresUp[ic];

    // derivatives with respect to upstream component densities
    for( integer jdof = 0; jdof < NC; ++jdof )
    {
      localFluxJacobian_dPresCompUp[TAG::NEXT *NC+ic][jdof+1]    =  dt * dCompFlux_dCompDensUp[ic][jdof];
      localFluxJacobian_dPresCompUp[TAG::CURRENT *NC+ic][jdof+1] = -dt * dCompFlux_dCompDensUp[ic][jdof];
    }
  }
}

template< integer NC >
void
FluxKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          integer const useTotalMassEquation,
          WellControls const & wellControls,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  using namespace compositionalMultiphaseUtilities;

  bool const isProducer = wellControls.isProducer();
  arrayView1d< real64 const > const & injection = wellControls.getInjectionStream();

  // loop over the well elements to compute the fluxes between elements
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {
    // create local work arrays
    real64 compFracUp[NC]{};
    real64 dCompFrac_dCompDensUp[NC][NC]{};

    real64 compFlux[NC]{};
    real64 dCompFlux_dRate[NC]{};
    real64 dCompFlux_dPresUp[NC]{};
    real64 dCompFlux_dCompDensUp[NC][NC]{};

    // Step 1) decide the upwind well element

    /*  currentConnRate < 0 flow from iwelem to iwelemNext
     *  currentConnRate > 0 flow from iwelemNext to iwelem
     *  With this convention, currentConnRate < 0 at the last connection for a producer
     *                        currentConnRate > 0 at the last connection for a injector
     */

    localIndex const iwelemNext = nextWellElemIndex[iwelem];
    real64 const currentConnRate = connRate[iwelem];
    localIndex iwelemUp = -1;

    if( iwelemNext < 0 && !isProducer ) // exit connection, injector
    {
      // we still need to define iwelemUp for Jacobian assembly
      iwelemUp = iwelem;

      // just copy the injection stream into compFrac
      for( integer ic = 0; ic < NC; ++ic )
      {
        compFracUp[ic] = injection[ic];
        for( integer jc = 0; jc < NC; ++jc )
        {
          dCompFrac_dCompDensUp[ic][jc] = 0.0;
        }
      }
    }
    else
    {
      // first set iwelemUp to the upstream cell
      if( ( iwelemNext < 0 && isProducer )  // exit connection, producer
          || currentConnRate < 0 ) // not an exit connection, iwelem is upstream
      {
        iwelemUp = iwelem;
      }
      else // not an exit connection, iwelemNext is upstream
      {
        iwelemUp = iwelemNext;
      }

      // copy the vars of iwelemUp into compFrac
      for( integer ic = 0; ic < NC; ++ic )
      {
        compFracUp[ic] = wellElemCompFrac[iwelemUp][ic];
        for( integer jc = 0; jc < NC; ++jc )
        {
          dCompFrac_dCompDensUp[ic][jc] = dWellElemCompFrac_dCompDens[iwelemUp][ic][jc];
        }
      }
    }

    // Step 2) compute upstream transport coefficient

    for( integer ic = 0; ic < NC; ++ic )
    {
      compFlux[ic]          = compFracUp[ic] * currentConnRate;
      dCompFlux_dRate[ic]   = compFracUp[ic];
      dCompFlux_dPresUp[ic] = 0.0; // none of these quantities depend on pressure
      for( integer jc = 0; jc < NC; ++jc )
      {
        dCompFlux_dCompDensUp[ic][jc] = dCompFrac_dCompDensUp[ic][jc] * currentConnRate;
      }
    }

    globalIndex const offsetUp = wellElemDofNumber[iwelemUp];
    globalIndex const offsetCurrent = wellElemDofNumber[iwelem];

    if( iwelemNext < 0 )  // exit connection
    {
      // for this case, we only need NC mass conservation equations
      // so we do not use the arrays initialized before the loop
      real64 oneSidedFlux[NC]{};
      real64 oneSidedFluxJacobian_dRate[NC][1]{};
      real64 oneSidedFluxJacobian_dPresCompUp[NC][NC+1]{};

      computeExit< NC >( dt,
                         compFlux,
                         dCompFlux_dRate,
                         dCompFlux_dPresUp,
                         dCompFlux_dCompDensUp,
                         oneSidedFlux,
                         oneSidedFluxJacobian_dRate,
                         oneSidedFluxJacobian_dPresCompUp );


      globalIndex oneSidedEqnRowIndices[NC]{};
      globalIndex oneSidedDofColIndices_dPresCompUp[NC+1]{};
      globalIndex oneSidedDofColIndices_dRate = 0;

      // jacobian indices
      for( integer ic = 0; ic < NC; ++ic )
      {
        // mass balance equations for all components
        oneSidedEqnRowIndices[ic] = offsetUp + ROFFSET::MASSBAL + ic - rankOffset;
      }

      // in the dof ordering used in this class, there are 1 pressure dofs
      // and NC compDens dofs before the rate dof in this block
      localIndex const dRateColOffset = COFFSET::DCOMP + NC;
      oneSidedDofColIndices_dRate = offsetCurrent + dRateColOffset;

      for( integer jdof = 0; jdof < NC+1; ++jdof )
      {
        // dofs are the **upstream** pressure and component densities
        oneSidedDofColIndices_dPresCompUp[jdof] = offsetUp + COFFSET::DPRES + jdof;
      }

      if( useTotalMassEquation > 0 )
      {
        // Apply equation/variable change transformation(s)
        real64 work[NC + 1]{};
        shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, 1, oneSidedFluxJacobian_dRate, work );
        shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC + 1, oneSidedFluxJacobian_dPresCompUp, work );
        shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, oneSidedFlux );
      }

      for( integer i = 0; i < NC; ++i )
      {
        if( oneSidedEqnRowIndices[i] >= 0 && oneSidedEqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                        &oneSidedDofColIndices_dRate,
                                                        oneSidedFluxJacobian_dRate[i],
                                                        1 );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                                            oneSidedDofColIndices_dPresCompUp,
                                                                            oneSidedFluxJacobian_dPresCompUp[i],
                                                                            NC+1 );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[oneSidedEqnRowIndices[i]], oneSidedFlux[i] );
        }
      }
    }
    else // not an exit connection
    {
      real64 localFlux[2*NC]{};
      real64 localFluxJacobian_dRate[2*NC][1]{};
      real64 localFluxJacobian_dPresCompUp[2*NC][NC+1]{};

      compute< NC >( dt,
                     compFlux,
                     dCompFlux_dRate,
                     dCompFlux_dPresUp,
                     dCompFlux_dCompDensUp,
                     localFlux,
                     localFluxJacobian_dRate,
                     localFluxJacobian_dPresCompUp );


      globalIndex eqnRowIndices[2*NC]{};
      globalIndex dofColIndices_dPresCompUp[NC+1]{};
      globalIndex dofColIndices_dRate = 0;

      globalIndex const offsetNext = wellElemDofNumber[iwelemNext];

      // jacobian indices
      for( integer ic = 0; ic < NC; ++ic )
      {
        // mass balance equations for all components
        eqnRowIndices[TAG::NEXT *NC+ic]    = offsetNext + ROFFSET::MASSBAL + ic - rankOffset;
        eqnRowIndices[TAG::CURRENT *NC+ic] = offsetCurrent + ROFFSET::MASSBAL + ic - rankOffset;
      }

      // in the dof ordering used in this class, there are 1 pressure dofs
      // and NC compDens dofs before the rate dof in this block
      localIndex const dRateColOffset = COFFSET::DCOMP + NC;
      dofColIndices_dRate = offsetCurrent + dRateColOffset;

      for( integer jdof = 0; jdof < NC+1; ++jdof )
      {
        // dofs are the **upstream** pressure and component densities
        dofColIndices_dPresCompUp[jdof] = offsetUp + COFFSET::DPRES + jdof;
      }

      if( useTotalMassEquation > 0 )
      {
        // Apply equation/variable change transformation(s)
        real64 work[NC + 1]{};
        shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC, 1, 2, localFluxJacobian_dRate, work );
        shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC, NC + 1, 2, localFluxJacobian_dPresCompUp, work );
        shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( NC, NC, 2, localFlux );
      }

      for( integer i = 0; i < 2*NC; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                        &dofColIndices_dRate,
                                                        localFluxJacobian_dRate[i],
                                                        1 );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                            dofColIndices_dPresCompUp,
                                                                            localFluxJacobian_dPresCompUp[i],
                                                                            NC+1 );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localFlux[i] );
        }
      }
    }
  } );
}

#define INST_FluxKernel( NC ) \
  template \
  void FluxKernel:: \
    launch< NC >( localIndex const size, \
                  globalIndex const rankOffset, \
                  integer const useTotalMassEquation, \
                  WellControls const & wellControls, \
                  arrayView1d< globalIndex const > const & wellElemDofNumber, \
                  arrayView1d< localIndex const > const & nextWellElemIndex, \
                  arrayView1d< real64 const > const & connRate, \
                  arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac, \
                  arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens, \
                  real64 const & dt, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 1 );
INST_FluxKernel( 2 );
INST_FluxKernel( 3 );
INST_FluxKernel( 4 );
INST_FluxKernel( 5 );

/******************************** PressureRelationKernel ********************************/

template< integer NC >
GEOS_HOST_DEVICE
void
PressureRelationKernel::
  compute( real64 const & gravCoef,
           real64 const & gravCoefNext,
           real64 const & pres,
           real64 const & presNext,
           real64 const & totalMassDens,
           real64 const & totalMassDensNext,
           real64 const & dTotalMassDens_dPres,
           real64 const & dTotalMassDens_dPresNext,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDens,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDensNext,
           real64 & localPresRel,
           real64 ( & localPresRelJacobian )[2*(NC+1)] )
{
  // local working variables and arrays
  real64 dAvgMassDens_dCompCurrent[NC]{};
  real64 dAvgMassDens_dCompNext[NC]{};

  // compute the average density at the interface between well elements
  real64 const avgMassDens = 0.5 * ( totalMassDensNext + totalMassDens );
  real64 const dAvgMassDens_dPresNext    = 0.5 * dTotalMassDens_dPresNext;
  real64 const dAvgMassDens_dPresCurrent = 0.5 * dTotalMassDens_dPres;
  for( integer ic = 0; ic < NC; ++ic )
  {
    dAvgMassDens_dCompNext[ic]    = 0.5 * dTotalMassDens_dCompDensNext[ic];
    dAvgMassDens_dCompCurrent[ic] = 0.5 * dTotalMassDens_dCompDens[ic];
  }

  // compute depth diff times acceleration
  real64 const gravD = gravCoefNext - gravCoef;

  // TODO: add friction and acceleration terms

  localPresRel = ( presNext - pres - avgMassDens * gravD );
  localPresRelJacobian[TAG::NEXT *(NC+1)]    = ( 1 - dAvgMassDens_dPresNext * gravD );
  localPresRelJacobian[TAG::CURRENT *(NC+1)] = ( -1 - dAvgMassDens_dPresCurrent * gravD );

  for( integer ic = 0; ic < NC; ++ic )
  {
    localPresRelJacobian[TAG::NEXT *(NC+1) + ic+1]    = -dAvgMassDens_dCompNext[ic] * gravD;
    localPresRelJacobian[TAG::CURRENT *(NC+1) + ic+1] = -dAvgMassDens_dCompCurrent[ic] * gravD;
  }
}

template< integer NC >
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
          arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres,
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens,
          bool & controlHasSwitched,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{

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
  real64 const & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dPresString() );
  arrayView1d< real64 const > const & dCurrentBHP_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dCompDensString() );

  arrayView1d< real64 const > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  arrayView1d< real64 const > const & dCurrentPhaseVolRate_dPres =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dPresString() );
  arrayView2d< real64 const > const & dCurrentPhaseVolRate_dCompDens =
    wellControls.getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dCompDensString() );
  arrayView1d< real64 const > const & dCurrentPhaseVolRate_dRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dRateString() );

  real64 const & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
  real64 const & dCurrentTotalVolRate_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dPresString() );
  arrayView1d< real64 const > const & dCurrentTotalVolRate_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dCompDensString() );
  real64 const & dCurrentTotalVolRate_dRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dRateString() );
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

      ControlEquationHelper::compute< NC >( rankOffset,
                                            newControl,
                                            targetPhaseIndex,
                                            targetBHP,
                                            targetPhaseRate,
                                            targetTotalRate,
                                            targetMassRate,
                                            currentBHP,
                                            dCurrentBHP_dPres,
                                            dCurrentBHP_dCompDens,
                                            currentPhaseVolRate,
                                            dCurrentPhaseVolRate_dPres,
                                            dCurrentPhaseVolRate_dCompDens,
                                            dCurrentPhaseVolRate_dRate,
                                            currentTotalVolRate,
                                            dCurrentTotalVolRate_dPres,
                                            dCurrentTotalVolRate_dCompDens,
                                            dCurrentTotalVolRate_dRate,
                                            massDensity,
                                            wellElemDofNumber[iwelemControl],
                                            localMatrix,
                                            localRhs );

      // TODO: for consistency, we should assemble here, not in compute...

    }
    else if( iwelemNext >= 0 ) // if iwelemNext >= 0, form momentum equation
    {

      real64 localPresRel = 0;
      real64 localPresRelJacobian[2*(NC+1)]{};

      compute< NC >( wellElemGravCoef[iwelem],
                     wellElemGravCoef[iwelemNext],
                     wellElemPressure[iwelem],
                     wellElemPressure[iwelemNext],
                     wellElemTotalMassDens[iwelem],
                     wellElemTotalMassDens[iwelemNext],
                     dWellElemTotalMassDens_dPres[iwelem],
                     dWellElemTotalMassDens_dPres[iwelemNext],
                     dWellElemTotalMassDens_dCompDens[iwelem],
                     dWellElemTotalMassDens_dCompDens[iwelemNext],
                     localPresRel,
                     localPresRelJacobian );


      // local working variables and arrays
      globalIndex dofColIndices[2*(NC+1)];

      globalIndex const eqnRowIndex = wellElemDofNumber[iwelem] + ROFFSET::CONTROL - rankOffset;
      dofColIndices[TAG::NEXT *(NC+1)]    = wellElemDofNumber[iwelemNext] + COFFSET::DPRES;
      dofColIndices[TAG::CURRENT *(NC+1)] = wellElemDofNumber[iwelem] + COFFSET::DPRES;

      for( integer ic = 0; ic < NC; ++ic )
      {
        dofColIndices[TAG::NEXT *(NC+1) + ic+1]    = wellElemDofNumber[iwelemNext] + COFFSET::DCOMP + ic;
        dofColIndices[TAG::CURRENT *(NC+1) + ic+1] = wellElemDofNumber[iwelem] + COFFSET::DCOMP + ic;
      }

      if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndex,
                                                                          dofColIndices,
                                                                          localPresRelJacobian,
                                                                          2 * (NC+1) );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndex], localPresRel );
      }
    }
  } );
  controlHasSwitched = ( switchControl.get() == 1 );
}

#define INST_PressureRelationKernel( NC ) \
  template \
  void PressureRelationKernel:: \
    launch< NC >( localIndex const size, \
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
                  arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres, \
                  arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens, \
                  bool & controlHasSwitched, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_PressureRelationKernel( 1 );
INST_PressureRelationKernel( 2 );
INST_PressureRelationKernel( 3 );
INST_PressureRelationKernel( 4 );
INST_PressureRelationKernel( 5 );


/******************************** PerforationKernel ********************************/

template< integer NC, integer NP >
GEOS_HOST_DEVICE
void
PerforationKernel::
  compute( bool const & disableReservoirToWellFlow,
           real64 const & resPres,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & resPhaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dResPhaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dResCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & resPhaseDens,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dResPhaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & resPhaseVisc,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dResPhaseVisc,
           arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > const & resPhaseCompFrac,
           arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > const & dResPhaseCompFrac,
           arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > const & resPhaseRelPerm,
           arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > const & dResPhaseRelPerm_dPhaseVolFrac,
           real64 const & wellElemGravCoef,
           real64 const & wellElemPres,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & wellElemCompDens,
           real64 const & wellElemTotalMassDens,
           real64 const & dWellElemTotalMassDens_dPres,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dWellElemTotalMassDens_dCompDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & wellElemCompFrac,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dWellElemCompFrac_dCompDens,
           real64 const & perfGravCoef,
           real64 const & trans,
           arraySlice1d< real64 > const & compPerfRate,
           arraySlice2d< real64 > const & dCompPerfRate_dPres,
           arraySlice3d< real64 > const & dCompPerfRate_dComp )
{
  using Deriv = multifluid::DerivativeOffset;

  // local working variables and arrays
  real64 pres[2]{};
  real64 dPres_dP[2]{};
  real64 dPres_dC[2][NC]{};
  real64 dFlux_dP[2]{};
  real64 dFlux_dC[2][NC]{};
  real64 dMult_dP[2]{};
  real64 dMult_dC[2][NC]{};
  real64 dPotDiff_dP[2]{};
  real64 dPotDiff_dC[2][NC]{};
  real64 multiplier[2]{};

  real64 dResTotalMob_dC[NC]{};
  real64 dDens_dC[NC]{};
  real64 dVisc_dC[NC]{};
  real64 dRelPerm_dC[NC]{};
  real64 dMob_dC[NC]{};
  real64 dCompFrac_dCompDens[NC]{};


  // Step 1: reset the perforation rates

  for( integer ic = 0; ic < NC; ++ic )
  {
    compPerfRate[ic] = 0.0;
    for( integer ke = 0; ke < 2; ++ke )
    {
      dCompPerfRate_dPres[ke][ic] = 0.0;
      for( integer jc = 0; jc < NC; ++jc )
      {
        dCompPerfRate_dComp[ke][ic][jc] = 0.0;
      }
    }
  }


  // Step 2: copy the variables from the reservoir and well element

  // a) get reservoir variables

  pres[TAG::RES] = resPres;
  dPres_dP[TAG::RES] = 1.0;
  multiplier[TAG::RES] = 1.0;

  // Here in the absence of a buoyancy term we assume that the reservoir cell is perforated at its center
  // TODO: add a buoyancy term for the reservoir side here


  // b) get well variables

  pres[TAG::WELL] = wellElemPres;
  dPres_dP[TAG::WELL] = 1.0;
  multiplier[TAG::WELL] = -1.0;

  real64 const gravD = ( perfGravCoef - wellElemGravCoef );

  pres[TAG::WELL] += wellElemTotalMassDens * gravD;
  dPres_dP[TAG::WELL] += dWellElemTotalMassDens_dPres * gravD;
  for( integer ic = 0; ic < NC; ++ic )
  {
    dPres_dC[TAG::WELL][ic] += dWellElemTotalMassDens_dCompDens[ic] * gravD;
  }


  // Step 3: compute potential difference

  real64 potDiff = 0.0;
  for( integer i = 0; i < 2; ++i )
  {
    potDiff += multiplier[i] * trans * pres[i];
    dPotDiff_dP[i] += multiplier[i] * trans * dPres_dP[i];

    for( integer ic = 0; ic < NC; ++ic )
    {
      dPotDiff_dC[i][ic] += multiplier[i] * trans * dPres_dC[i][ic];
    }
  }


  // Step 4: upwinding based on the flow direction

  real64 flux = 0.0;
  if( potDiff >= 0 )  // ** reservoir cell is upstream **
  {

    // loop over phases, compute and upwind phase flux
    // and sum contributions to each component's perforation rate
    for( integer ip = 0; ip < NP; ++ip )
    {

      // skip the rest of the calculation if the phase is absent
      // or if crossflow is disabled for injectors
      bool const phaseExists = (resPhaseVolFrac[ip] > 0);
      if( !phaseExists || disableReservoirToWellFlow )
      {
        continue;
      }

      // here, we have to recompute the reservoir phase mobility (not including density)

      // density
      real64 const resDens = resPhaseDens[ip];
      real64 const dResDens_dP  = dResPhaseDens[ip][Deriv::dP];
      applyChainRule( NC, dResCompFrac_dCompDens,
                      dResPhaseDens[ip],
                      dDens_dC,
                      Deriv::dC );

      // viscosity
      real64 const resVisc = resPhaseVisc[ip];
      real64 const dResVisc_dP  = dResPhaseVisc[ip][Deriv::dP];
      applyChainRule( NC, dResCompFrac_dCompDens,
                      dResPhaseVisc[ip],
                      dVisc_dC,
                      Deriv::dC );

      // relative permeability
      real64 const resRelPerm = resPhaseRelPerm[ip];
      real64 dResRelPerm_dP = 0.0;
      for( integer jc = 0; jc < NC; ++jc )
      {
        dRelPerm_dC[jc] = 0;
      }

      for( integer jp = 0; jp < NP; ++jp )
      {
        real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[ip][jp];
        dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac[jp][Deriv::dP];

        for( integer jc = 0; jc < NC; ++jc )
        {
          dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac[jp][Deriv::dC+jc];
        }
      }

      // compute the reservoir phase mobility, including phase density
      real64 const resPhaseMob = resDens * resRelPerm / resVisc;
      real64 const dResPhaseMob_dPres = dResRelPerm_dP * resDens / resVisc
                                        + resPhaseMob * (dResDens_dP / resDens - dResVisc_dP / resVisc);
      for( integer jc = 0; jc < NC; ++jc )
      {
        dMob_dC[jc] = dRelPerm_dC[jc] * resDens / resVisc
                      + resPhaseMob * (dDens_dC[jc] / resDens - dVisc_dC[jc] / resVisc);
      }

      // compute the phase flux and derivatives using upstream cell mobility
      flux = resPhaseMob * potDiff;
      dFlux_dP[TAG::RES]  = dResPhaseMob_dPres * potDiff + resPhaseMob * dPotDiff_dP[TAG::RES];
      dFlux_dP[TAG::WELL] = resPhaseMob *  dPotDiff_dP[TAG::WELL];

      for( integer ic = 0; ic < NC; ++ic )
      {
        dFlux_dC[TAG::RES][ic] = dMob_dC[ic] * potDiff + resPhaseMob * dPotDiff_dC[TAG::RES][ic];
        dFlux_dC[TAG::WELL][ic] = resPhaseMob * dPotDiff_dC[TAG::WELL][ic];
      }

      // increment component fluxes
      for( integer ic = 0; ic < NC; ++ic )
      {
        compPerfRate[ic] += flux * resPhaseCompFrac[ip][ic];

        dCompPerfRate_dPres[TAG::RES][ic]  += resPhaseCompFrac[ip][ic] * dFlux_dP[TAG::RES];
        dCompPerfRate_dPres[TAG::RES][ic]  += dResPhaseCompFrac[ip][ic][Deriv::dP] * flux;
        dCompPerfRate_dPres[TAG::WELL][ic] += resPhaseCompFrac[ip][ic] * dFlux_dP[TAG::WELL];

        applyChainRule( NC,
                        dResCompFrac_dCompDens,
                        dResPhaseCompFrac[ip][ic],
                        dCompFrac_dCompDens,
                        Deriv::dC );

        for( integer jc = 0; jc < NC; ++jc )
        {
          dCompPerfRate_dComp[TAG::RES][ic][jc]  += dFlux_dC[TAG::RES][jc] * resPhaseCompFrac[ip][ic];
          dCompPerfRate_dComp[TAG::RES][ic][jc]  += flux * dCompFrac_dCompDens[jc];
          dCompPerfRate_dComp[TAG::WELL][ic][jc] += dFlux_dC[TAG::WELL][jc] * resPhaseCompFrac[ip][ic];
        }
      }
    }
  }
  else // ** well is upstream **
  {

    real64 resTotalMob     = 0.0;
    real64 dResTotalMob_dP = 0.0;

    // we re-compute here the total mass (when useMass == 1) or molar (when useMass == 0) density
    real64 wellElemTotalDens = 0;
    for( integer ic = 0; ic < NC; ++ic )
    {
      wellElemTotalDens += wellElemCompDens[ic];
    }

    // first, compute the reservoir total mobility (excluding phase density)
    for( integer ip = 0; ip < NP; ++ip )
    {

      // skip the rest of the calculation if the phase is absent
      bool const phaseExists = (resPhaseVolFrac[ip] > 0);
      if( !phaseExists )
      {
        continue;
      }

      // viscosity
      real64 const resVisc = resPhaseVisc[ip];
      real64 const dResVisc_dP  = dResPhaseVisc[ip][Deriv::dP];
      applyChainRule( NC, dResCompFrac_dCompDens,
                      dResPhaseVisc[ip],
                      dVisc_dC,
                      Deriv::dC );

      // relative permeability
      real64 const resRelPerm = resPhaseRelPerm[ip];
      real64 dResRelPerm_dP = 0.0;
      for( integer jc = 0; jc < NC; ++jc )
      {
        dRelPerm_dC[jc] = 0;
      }

      for( integer jp = 0; jp < NP; ++jp )
      {
        real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[ip][jp];
        dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac[jp][Deriv::dP];

        for( integer jc = 0; jc < NC; ++jc )
        {
          dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac[jp][Deriv::dC+jc];
        }
      }

      // increment total mobility
      resTotalMob     += resRelPerm / resVisc;
      dResTotalMob_dP += ( dResRelPerm_dP * resVisc - resRelPerm * dResVisc_dP )
                         / ( resVisc * resVisc );
      for( integer ic = 0; ic < NC; ++ic )
      {
        dResTotalMob_dC[ic] += ( dRelPerm_dC[ic] * resVisc - resRelPerm * dVisc_dC[ic] )
                               / ( resVisc * resVisc );
      }
    }

    // compute a potdiff multiplier = wellElemTotalDens * resTotalMob
    // wellElemTotalDens is a mass density if useMass == 1 and a molar density otherwise
    real64 const mult   = wellElemTotalDens * resTotalMob;
    dMult_dP[TAG::RES]  = wellElemTotalDens * dResTotalMob_dP;
    dMult_dP[TAG::WELL] = 0.0; // because totalDens does not depend on pressure

    for( integer ic = 0; ic < NC; ++ic )
    {
      dMult_dC[TAG::RES][ic]  = wellElemTotalDens * dResTotalMob_dC[ic];
      dMult_dC[TAG::WELL][ic] = resTotalMob;
    }

    // compute the volumetric flux and derivatives using upstream cell mobility
    flux = mult * potDiff;
    dFlux_dP[TAG::RES]  = dMult_dP[TAG::RES] * potDiff + mult * dPotDiff_dP[TAG::RES];
    dFlux_dP[TAG::WELL] = dMult_dP[TAG::WELL] * potDiff + mult * dPotDiff_dP[TAG::WELL];

    for( integer ic = 0; ic < NC; ++ic )
    {
      dFlux_dC[TAG::RES][ic]  = dMult_dC[TAG::RES][ic] * potDiff + mult * dPotDiff_dC[TAG::RES][ic];
      dFlux_dC[TAG::WELL][ic] = dMult_dC[TAG::WELL][ic] * potDiff + mult * dPotDiff_dC[TAG::WELL][ic];
    }

    // compute component fluxes
    for( integer ic = 0; ic < NC; ++ic )
    {
      compPerfRate[ic] += wellElemCompFrac[ic] * flux;
      dCompPerfRate_dPres[TAG::RES][ic]  = wellElemCompFrac[ic] * dFlux_dP[TAG::RES];
      dCompPerfRate_dPres[TAG::WELL][ic] = wellElemCompFrac[ic] * dFlux_dP[TAG::WELL];

      for( integer jc = 0; jc < NC; ++jc )
      {
        dCompPerfRate_dComp[TAG::RES][ic][jc]  += wellElemCompFrac[ic] * dFlux_dC[TAG::RES][jc];
        dCompPerfRate_dComp[TAG::WELL][ic][jc] += wellElemCompFrac[ic] * dFlux_dC[TAG::WELL][jc];
        dCompPerfRate_dComp[TAG::WELL][ic][jc] += dWellElemCompFrac_dCompDens[ic][jc] * flux;
      }
    }
  }
}

template< integer NC, integer NP >
void
PerforationKernel::
  launch( localIndex const size,
          bool const disableReservoirToWellFlow,
          ElementViewConst< arrayView1d< real64 const > > const & resPres,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dResPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dResCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dResPhaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseVisc,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dResPhaseVisc,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & resPhaseCompFrac,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dResPhaseCompFrac,
          ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & resPhaseRelPerm,
          ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const & dResPhaseRelPerm_dPhaseVolFrac,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 const > const & wellElemPres,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens,
          arrayView1d< real64 const > const & wellElemTotalMassDens,
          arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres,
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< localIndex const > const & perfWellElemIndex,
          arrayView1d< real64 const > const & perfTrans,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView2d< real64 > const & compPerfRate,
          arrayView3d< real64 > const & dCompPerfRate_dPres,
          arrayView4d< real64 > const & dCompPerfRate_dComp )
{

  // loop over the perforations to compute the perforation rates
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iperf )
  {

    // get the index of the reservoir elem
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];

    // get the index of the well elem
    localIndex const iwelem = perfWellElemIndex[iperf];

    compute< NC, NP >( disableReservoirToWellFlow,
                       resPres[er][esr][ei],
                       resPhaseVolFrac[er][esr][ei],
                       dResPhaseVolFrac[er][esr][ei],
                       dResCompFrac_dCompDens[er][esr][ei],
                       resPhaseDens[er][esr][ei][0],
                       dResPhaseDens[er][esr][ei][0],
                       resPhaseVisc[er][esr][ei][0],
                       dResPhaseVisc[er][esr][ei][0],
                       resPhaseCompFrac[er][esr][ei][0],
                       dResPhaseCompFrac[er][esr][ei][0],
                       resPhaseRelPerm[er][esr][ei][0],
                       dResPhaseRelPerm_dPhaseVolFrac[er][esr][ei][0],
                       wellElemGravCoef[iwelem],
                       wellElemPres[iwelem],
                       wellElemCompDens[iwelem],
                       wellElemTotalMassDens[iwelem],
                       dWellElemTotalMassDens_dPres[iwelem],
                       dWellElemTotalMassDens_dCompDens[iwelem],
                       wellElemCompFrac[iwelem],
                       dWellElemCompFrac_dCompDens[iwelem],
                       perfGravCoef[iperf],
                       perfTrans[iperf],
                       compPerfRate[iperf],
                       dCompPerfRate_dPres[iperf],
                       dCompPerfRate_dComp[iperf] );

  } );
}

#define INST_PerforationKernel( NC, NP ) \
  template \
  void PerforationKernel:: \
    launch< NC, NP >( localIndex const size, \
                      bool const disableReservoirToWellFlow, \
                      ElementViewConst< arrayView1d< real64 const > > const & resPres, \
                      ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac, \
                      ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dResPhaseVolFrac, \
                      ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dResCompFrac_dCompDens, \
                      ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseDens, \
                      ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dResPhaseDens, \
                      ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & resPhaseVisc, \
                      ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dResPhaseVisc, \
                      ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & resPhaseCompFrac, \
                      ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dResPhaseCompFrac, \
                      ElementViewConst< arrayView3d< real64 const, relperm::USD_RELPERM > > const & resPhaseRelPerm, \
                      ElementViewConst< arrayView4d< real64 const, relperm::USD_RELPERM_DS > > const & dResPhaseRelPerm_dPhaseVolFrac, \
                      arrayView1d< real64 const > const & wellElemGravCoef, \
                      arrayView1d< real64 const > const & wellElemPres, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompDens, \
                      arrayView1d< real64 const > const & wellElemTotalMassDens, \
                      arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres, \
                      arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac, \
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens, \
                      arrayView1d< real64 const > const & perfGravCoef, \
                      arrayView1d< localIndex const > const & perfWellElemIndex, \
                      arrayView1d< real64 const > const & perfTrans, \
                      arrayView1d< localIndex const > const & resElementRegion, \
                      arrayView1d< localIndex const > const & resElementSubRegion, \
                      arrayView1d< localIndex const > const & resElementIndex, \
                      arrayView2d< real64 > const & compPerfRate, \
                      arrayView3d< real64 > const & dCompPerfRate_dPres, \
                      arrayView4d< real64 > const & dCompPerfRate_dComp )

INST_PerforationKernel( 1, 2 );
INST_PerforationKernel( 2, 2 );
INST_PerforationKernel( 3, 2 );
INST_PerforationKernel( 4, 2 );
INST_PerforationKernel( 5, 2 );
INST_PerforationKernel( 1, 3 );
INST_PerforationKernel( 2, 3 );
INST_PerforationKernel( 3, 3 );
INST_PerforationKernel( 4, 3 );
INST_PerforationKernel( 5, 3 );

/******************************** AccumulationKernel ********************************/

template< integer NC >
GEOS_HOST_DEVICE
void
AccumulationKernel::
  compute( integer const numPhases,
           real64 const & volume,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens,
           arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac,
           arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > const & dPhaseCompFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac_n,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens_n,
           arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac_n,
           real64 ( & localAccum )[NC],
           real64 ( & localAccumJacobian )[NC][NC + 1] )
{
  using Deriv = multifluid::DerivativeOffset;

  // temporary work arrays
  real64 dPhaseAmount_dC[NC]{};
  real64 dPhaseCompFrac_dC[NC]{};

  // reset the local values
  for( integer i = 0; i < NC; ++i )
  {
    localAccum[i] = 0.0;
    for( integer j = 0; j < NC+1; ++j )
    {
      localAccumJacobian[i][j] = 0.0;
    }
  }

  // sum contributions to component accumulation from each phase
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    real64 const phaseAmountNew = volume * phaseVolFrac[ip] * phaseDens[ip];
    real64 const phaseAmount_n = volume * phaseVolFrac_n[ip] * phaseDens_n[ip];

    real64 const dPhaseAmount_dP = volume * ( dPhaseVolFrac[ip][Deriv::dP] * phaseDens[ip]
                                              + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dP] );

    // assemble density dependence
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens[ip], dPhaseAmount_dC, Deriv::dC );
    for( integer jc = 0; jc < NC; ++jc )
    {
      dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                            + phaseDens[ip] * dPhaseVolFrac[ip][Deriv::dC+jc];
      dPhaseAmount_dC[jc] *= volume;
    }

    // ic - index of component whose conservation equation is assembled
    // (i.e. row number in local matrix)
    for( integer ic = 0; ic < NC; ++ic )
    {
      real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
      real64 const phaseCompAmount_n = phaseAmount_n * phaseCompFrac_n[ip][ic];

      real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                         + phaseAmountNew * dPhaseCompFrac[ip][ic][Deriv::dP];

      localAccum[ic] += phaseCompAmountNew - phaseCompAmount_n;
      localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

      // jc - index of component w.r.t. whose compositional var the derivative is being taken
      // (i.e. col number in local matrix)

      // assemble phase composition dependence
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac[ip][ic], dPhaseCompFrac_dC, Deriv::dC );
      for( integer jc = 0; jc < NC; ++jc )
      {
        real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                           + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
        localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
      }
    }
  }
}

template< integer NC >
void
AccumulationKernel::
  launch( localIndex const size,
          integer const numPhases,
          globalIndex const rankOffset,
          integer const useTotalMassEquation,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dWellElemPhaseDens,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac,
          arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dWellElemPhaseCompFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac_n,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens_n,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac_n,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{

  using namespace compositionalMultiphaseUtilities;

  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {

    if( wellElemGhostRank[iwelem] >= 0 )
    {
      return;
    }

    real64 localAccum[NC]{};
    real64 localAccumJacobian[NC][NC+1]{};

    compute< NC >( numPhases,
                   wellElemVolume[iwelem],
                   wellElemPhaseVolFrac[iwelem],
                   dWellElemPhaseVolFrac[iwelem],
                   dWellElemCompFrac_dCompDens[iwelem],
                   wellElemPhaseDens[iwelem][0],
                   dWellElemPhaseDens[iwelem][0],
                   wellElemPhaseCompFrac[iwelem][0],
                   dWellElemPhaseCompFrac[iwelem][0],
                   wellElemPhaseVolFrac_n[iwelem],
                   wellElemPhaseDens_n[iwelem][0],
                   wellElemPhaseCompFrac_n[iwelem][0],
                   localAccum,
                   localAccumJacobian );

    // set the equation row indices to be the mass balance equations for all components
    localIndex eqnRowIndices[NC]{};
    for( integer ic = 0; ic < NC; ++ic )
    {
      eqnRowIndices[ic] = wellElemDofNumber[iwelem] + ROFFSET::MASSBAL + ic - rankOffset;
    }

    // set DOF col indices for this block
    globalIndex dofColIndices[NC+1]{};
    for( integer idof = 0; idof < NC+1; ++idof )
    {
      dofColIndices[idof] = wellElemDofNumber[iwelem] + COFFSET::DPRES + idof;
    }

    if( useTotalMassEquation > 0 )
    {
      // Apply equation/variable change transformation(s)
      real64 work[NC + 1];
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC + 1, localAccumJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, localAccum );
    }

    // add contribution to residual and jacobian
    for( integer ic = 0; ic < NC; ++ic )
    {
      localRhs[eqnRowIndices[ic]] += localAccum[ic];
      localMatrix.addToRow< serialAtomic >( eqnRowIndices[ic],
                                            dofColIndices,
                                            localAccumJacobian[ic],
                                            NC+1 );
    }
  } );
}

#define INST_AccumulationKernel( NC ) \
  template \
  void AccumulationKernel:: \
    launch< NC >( localIndex const size, \
                  integer const numPhases, \
                  globalIndex const rankOffset, \
                  integer const useTotalMassEquation, \
                  arrayView1d< globalIndex const > const & wellElemDofNumber, \
                  arrayView1d< integer const > const & wellElemGhostRank, \
                  arrayView1d< real64 const > const & wellElemVolume, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac, \
                  arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac, \
                  arrayView3d< real64 const, compflow::USD_COMP_DC > const & dWellElemCompFrac_dCompDens, \
                  arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens, \
                  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dWellElemPhaseDens, \
                  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac, \
                  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dWellElemPhaseCompFrac, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac_n, \
                  arrayView3d< real64 const, multifluid::USD_PHASE > const & wellElemPhaseDens_n, \
                  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac_n, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_AccumulationKernel( 1 );
INST_AccumulationKernel( 2 );
INST_AccumulationKernel( 3 );
INST_AccumulationKernel( 4 );
INST_AccumulationKernel( 5 );

/******************************** VolumeBalanceKernel ********************************/

template< integer NC >
GEOS_HOST_DEVICE
void
VolumeBalanceKernel::
  compute( integer const numPhases,
           real64 const & volume,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac,
           real64 & localVolBalance,
           real64 ( & localVolBalanceJacobian )[NC+1] )
{
  using Deriv = multifluid::DerivativeOffset;

  localVolBalance = 1.0;
  for( integer ic = 0; ic < NC+1; ++ic )
  {
    localVolBalanceJacobian[ic] = 0.0;
  }

  // sum contributions to component accumulation from each phase
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    localVolBalance -= phaseVolFrac[ip];
    localVolBalanceJacobian[0] -= dPhaseVolFrac[ip][Deriv::dP];

    for( integer jc = 0; jc < NC; ++jc )
    {
      localVolBalanceJacobian[jc + 1] -= dPhaseVolFrac[ip][Deriv::dC+jc];
    }
  }

  // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
  for( integer idof = 0; idof < NC+1; ++idof )
  {
    localVolBalanceJacobian[idof] *= volume;
  }
  localVolBalance *= volume;
}

template< integer NC >
void
VolumeBalanceKernel::
  launch( localIndex const size,
          integer const numPhases,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac,
          arrayView1d< real64 const > const & wellElemVolume,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
  {

    if( wellElemGhostRank[iwelem] >= 0 )
    {
      return;
    }

    real64 localVolBalance = 1.0;
    real64 localVolBalanceJacobian[NC+1]{};

    compute< NC >( numPhases,
                   wellElemVolume[iwelem],
                   wellElemPhaseVolFrac[iwelem],
                   dWellElemPhaseVolFrac[iwelem],
                   localVolBalance,
                   localVolBalanceJacobian );

    // get equation/dof indices
    localIndex const localVolBalanceEqnIndex = wellElemDofNumber[iwelem] - rankOffset + ROFFSET::MASSBAL + NC;
    globalIndex localVolBalanceDOF[NC+1]{};
    for( integer jdof = 0; jdof < NC+1; ++jdof )
    {
      localVolBalanceDOF[jdof] = wellElemDofNumber[iwelem] + COFFSET::DPRES + jdof;
    }

    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localVolBalanceEqnIndex,
                                                              localVolBalanceDOF,
                                                              localVolBalanceJacobian,
                                                              NC+1 );
    localRhs[localVolBalanceEqnIndex] += localVolBalance;
  } );
}

#define INST_VolumeBalanceKernel( NC ) \
  template \
  void VolumeBalanceKernel:: \
    launch< NC >( localIndex const size, \
                  integer const numPhases, \
                  globalIndex const rankOffset, \
                  arrayView1d< globalIndex const > const & wellElemDofNumber, \
                  arrayView1d< integer const > const & wellElemGhostRank, \
                  arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac, \
                  arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac, \
                  arrayView1d< real64 const > const & wellElemVolume, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_VolumeBalanceKernel( 1 );
INST_VolumeBalanceKernel( 2 );
INST_VolumeBalanceKernel( 3 );
INST_VolumeBalanceKernel( 4 );
INST_VolumeBalanceKernel( 5 );

/******************************** PresTempCompFracInitializationKernel ********************************/

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

  stackArray1d< real64, MAX_NUM_COMP > avgCompFrac( numComps );
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
