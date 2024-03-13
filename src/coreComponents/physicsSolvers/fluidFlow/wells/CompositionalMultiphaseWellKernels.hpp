/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"

namespace geos
{

namespace compositionalMultiphaseWellKernels
{

using namespace constitutive;

static constexpr real64 minDensForDivision = 1e-10;

// tag to access well and reservoir elements in perforation rates computation
struct SubRegionTag
{
  static constexpr integer RES  = 0;
  static constexpr integer WELL = 1;
};

// tag to access the next and current well elements of a connection
struct ElemTag
{
  static constexpr integer CURRENT = 0;
  static constexpr integer NEXT    = 1;
};

// define the column offset of the derivatives
struct ColOffset
{
  static constexpr integer DPRES = 0;
  static constexpr integer DCOMP = 1;
};

template < integer NC >
struct ColOffset_WellJac
{
  static constexpr integer dP = 0;
  static constexpr integer dC = 1;
  static constexpr integer dQ = dC + NC;
  static constexpr integer dT = dQ+1;
};

// define the row offset of the residual equations
struct RowOffset
{
  static constexpr integer CONTROL = 0;
  static constexpr integer MASSBAL = 1;
};

template < integer NC >
struct RowOffset_WellJac
{
  static constexpr integer CONTROL   = 0;
  static constexpr integer MASSBAL   = 1;
  static constexpr integer VOLBAL    = MASSBAL + NC;
  static constexpr integer ENERGYBAL = VOLBAL+1;
};
/******************************** ControlEquationHelper ********************************/
struct ControlEquationHelper
{
  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

  GEOS_HOST_DEVICE
  inline
  static
  void
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
                 WellControls::Control & newControl );

  template< integer NC , integer IS_THERMAL >
  GEOS_HOST_DEVICE
  inline
  static void
  compute( globalIndex const rankOffset,
           WellControls::Control const currentControl,
           integer const targetPhaseIndex,
           real64 const & targetBHP,
           real64 const & targetPhaseRate,
           real64 const & targetTotalRate,
           real64 const & targetMassRate,
           real64 const & currentBHP,
           arrayView1d< real64 const > const & dCurrentBHP,
           real64 const & dCurrentBHP_dPres,
           arrayView1d< real64 const > const & dCurrentBHP_dCompDens,
           arrayView1d< real64 const > const & currentPhaseVolRate,
           arrayView2d< real64 const > const & dCurrentPhaseVolRate,
           arrayView1d< real64 const > const & dCurrentPhaseVolRate_dPres,
           arrayView2d< real64 const > const & dCurrentPhaseVolRate_dCompDens,
           arrayView1d< real64 const > const & dCurrentPhaseVolRate_dRate,
           real64 const & currentTotalVolRate,
           arrayView1d< real64 const > const & dCurrentTotalVolRate,
           real64 const & dCurrentTotalVolRate_dPres,
           arrayView1d< real64 const > const & dCurrentTotalVolRate_dCompDens,
           real64 const & dCurrentTotalVolRate_dRate,
           real64 const & massDensity,
           globalIndex const dofNumber,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs );

};

/******************************** FluxKernel ********************************/

struct FluxKernel
{

  using TAG = compositionalMultiphaseWellKernels::ElemTag;
  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

  template< integer NC >
  GEOS_HOST_DEVICE
  inline
  static void
    computeExit( real64 const & dt,
                 real64 const ( &compFlux )[NC],
                 real64 const ( &dCompFlux_dRate )[NC],
                 real64 const ( &dCompFlux_dPresUp )[NC],
                 real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
                 real64 ( &oneSidedFlux )[NC],
                 real64 ( &oneSidedFluxJacobian_dRate )[NC][1],
                 real64 ( &oneSidedFluxJacobian_dPresCompUp )[NC][NC + 1] );

  template< integer NC >
  GEOS_HOST_DEVICE
  inline
  static void
    compute( real64 const & dt,
             real64 const ( &compFlux )[NC],
             real64 const ( &dCompFlux_dRate )[NC],
             real64 const ( &dCompFlux_dPresUp )[NC],
             real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
             real64 ( &localFlux )[2*NC],
             real64 ( &localFluxJacobian_dRate )[2*NC][1],
             real64 ( &localFluxJacobian_dPresCompUp )[2*NC][NC + 1] );

  template< integer NC >
  static void
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
          arrayView1d< real64 > const & localRhs );

};

/******************************** PressureRelationKernel ********************************/

struct PressureRelationKernel
{
  using Deriv = multifluid::DerivativeOffset;
  using TAG = compositionalMultiphaseWellKernels::ElemTag;
  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

  template< integer NC ,integer IS_THERMAL >
  GEOS_HOST_DEVICE
  inline
  static void
    compute( real64 const & gravCoef,
             real64 const & gravCoefNext,
             real64 const & pres,
             real64 const & presNext,
             real64 const & totalMassDens,
             real64 const & totalMassDensNext,
             arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens,
             arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDensNext,
             real64 const & dTotalMassDens_dPres,
             real64 const & dTotalMassDens_dPresNext,
             arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDens,
             arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dTotalMassDens_dCompDensNext,
             real64 & localPresRel,
             real64 ( &localPresRelJacobian )[2*(NC+1+IS_THERMAL)] );

  template< integer NC , integer IS_THERMAL >
  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          bool const isLocallyOwned,
          localIndex const iwelemControl,
          bool const isThermal,
          integer const targetPhaseIndex,
          WellControls const & wellControls,
          real64 const & timeAtEndOfStep,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & wellElemTotalMassDens,
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens,
          arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres,
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens_dCompDens,
          bool & controlHasSwitched,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

/******************************** PerforationKernel ********************************/

struct PerforationKernel
{

  using TAG = compositionalMultiphaseWellKernels::SubRegionTag;

  using CompFlowAccessors =
    StencilAccessors< fields::flow::pressure,
                      fields::flow::phaseVolumeFraction,
                      fields::flow::dPhaseVolumeFraction,
                      fields::flow::dGlobalCompFraction_dGlobalCompDensity >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseViscosity,
                              fields::multifluid::dPhaseViscosity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  using RelPermAccessors =
    StencilMaterialAccessors< RelativePermeabilityBase,
                              fields::relperm::phaseRelPerm,
                              fields::relperm::dPhaseRelPerm_dPhaseVolFraction >;


  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;


  template< integer NC, integer NP , integer IS_THERMAL>
  GEOS_HOST_DEVICE
  inline
  static void
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
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dWellElemTotalMassDens,
           real64 const & dWellElemTotalMassDens_dPres,
           arraySlice1d< real64 const, compflow::USD_FLUID_DC - 1 > const & dWellElemTotalMassDens_dCompDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & wellElemCompFrac,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dWellElemCompFrac_dCompDens,
           real64 const & perfGravCoef,
           real64 const & trans,
           arraySlice1d< real64 > const & compPerfRate,
           arraySlice3d< real64 > const & dCompPerfRate,
           arraySlice2d< real64 > const & dCompPerfRate_dPres,
           arraySlice3d< real64 > const & dCompPerfRate_dComp );

  template< integer NC, integer NP, integer IS_THERMAL >
  static void
  launch( localIndex const size,
          bool const disableReservoirToWellFlow,
          ElementViewConst< arrayView1d< real64 const > > const & resPres,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dResPhaseVolFrac_dComp,
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
          arrayView2d< real64 const, compflow::USD_FLUID_DC > const & dWellElemTotalMassDens,
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
          arrayView4d< real64 > const & dCompPerfRate,
          arrayView3d< real64 > const & dCompPerfRate_dPres,
          arrayView4d< real64 > const & dCompPerfRate_dComp );

};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

  template< integer NC >
  GEOS_HOST_DEVICE
  inline
  static void
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
             real64 ( &localAccum )[NC],
             real64 ( &localAccumJacobian )[NC][NC + 1] );

  template< integer NC >
  static void
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
          arrayView1d< real64 > const & localRhs );

};

/******************************** VolumeBalanceKernel ********************************/

struct VolumeBalanceKernel
{

  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

  template< integer NC >
  GEOS_HOST_DEVICE
  inline
  static void
    compute( integer const numPhases,
             real64 const & volume,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac,
             real64 & localVolBalance,
             real64 ( &localVolBalanceJacobian )[NC+1] );

  template< integer NC >
  static void
  launch( localIndex const size,
          integer const numPhases,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dWellElemPhaseVolFrac,
          arrayView1d< real64 const > const & wellElemVolume,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

/******************************** PresTempCompFracInitializationKernel ********************************/

struct PresTempCompFracInitializationKernel
{

  using CompFlowAccessors =
    StencilAccessors< fields::flow::pressure,
                      fields::flow::temperature,
                      fields::flow::globalCompDensity,
                      fields::flow::phaseVolumeFraction >;

  using MultiFluidAccessors =
    StencilMaterialAccessors< MultiFluidBase,
                              fields::multifluid::phaseMassDensity >;


  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  static void
  launch( localIndex const perforationSize,
          localIndex const subRegionSize,
          integer const numComponents,
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
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompFrac );

};

/******************************** CompDensInitializationKernel ********************************/

struct CompDensInitializationKernel
{

  static void
  launch( localIndex const subRegionSize,
          integer const numComponents,
          arrayView2d< real64 const, compflow::USD_COMP > const & wellElemCompFrac,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & wellElemTotalDens,
          arrayView2d< real64, compflow::USD_COMP > const & wellElemCompDens );

};

/******************************** RateInitializationKernel ********************************/

struct RateInitializationKernel
{

  static void
  launch( localIndex const subRegionSize,
          integer const targetPhaseIndex,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView2d< real64 const, multifluid::USD_FLUID > const & totalDens,
          arrayView1d< real64 > const & connRate );

};


/******************************** TotalMassDensityKernel ****************************/

/**
 * @class TotalMassDensityKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_PHASE number of fluid phases
 * @brief Define the interface for the property kernel in charge of computing the total mass density
 */
template< integer NUM_COMP, integer NUM_PHASE >
class TotalMassDensityKernel : public isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >
{
public:

  using Base = isothermalCompositionalMultiphaseBaseKernels::PropertyKernelBase< NUM_COMP >;
  using Base::numComp;

  /// Compile time value for the number of phases
  static constexpr integer numPhase = NUM_PHASE;

  /**
   * @brief Constructor
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  TotalMassDensityKernel( ObjectManagerBase & subRegion,
                          MultiFluidBase const & fluid )
    : Base(),
    m_phaseVolFrac( subRegion.getField< fields::well::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::well::dPhaseVolumeFraction >() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseMassDens( fluid.phaseMassDensity() ),
    m_dPhaseMassDens( fluid.dPhaseMassDensity() ),
    m_totalMassDens( subRegion.getField< fields::well::totalMassDensity >() ),
    m_dTotalMassDens( subRegion.getField< fields::well::dTotalMassDensity >() ),
    m_dTotalMassDens_dPres( subRegion.getField< fields::well::dTotalMassDensity_dPressure >() ),
    m_dTotalMassDens_dCompDens( subRegion.getField< fields::well::dTotalMassDensity_dGlobalCompDensity >() )
  {}

  /**
   * @brief Compute the total mass density in an element
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[in] totalMassDensityKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void compute( localIndex const ei,
                FUNC && totalMassDensityKernelOp = NoOpFunc{} ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseMassDens = m_phaseMassDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseMassDens = m_dPhaseMassDens[ei][0];

    real64 & totalMassDens = m_totalMassDens[ei];
    arraySlice1d< real64, compflow::USD_FLUID_DC - 1 > dTotalMassDens = m_dTotalMassDens[ei];
    real64 & dTotalMassDens_dPres = m_dTotalMassDens_dPres[ei];
    arraySlice1d< real64, compflow::USD_FLUID_DC - 1 > dTotalMassDens_dCompDens = m_dTotalMassDens_dCompDens[ei];
    real64 dMassDens_dC[numComp]{};

    totalMassDens = 0.0;
    dTotalMassDens[Deriv::dP]=0.0;
    dTotalMassDens_dPres = 0.0;
    dTotalMassDens[Deriv::dP]=0.0;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      dTotalMassDens_dCompDens[ic] = 0.0;
      dTotalMassDens[Deriv::dC+ic]=0.0;
    }

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      totalMassDens += phaseVolFrac[ip] * phaseMassDens[ip];
      dTotalMassDens_dPres += dPhaseVolFrac[ip][Deriv::dP] * phaseMassDens[ip] + phaseVolFrac[ip] * dPhaseMassDens[ip][Deriv::dP];
      dTotalMassDens[Deriv::dP] += dPhaseVolFrac[ip][Deriv::dP] * phaseMassDens[ip] + phaseVolFrac[ip] * dPhaseMassDens[ip][Deriv::dP];

      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseMassDens[ip], dMassDens_dC, Deriv::dC );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dTotalMassDens_dCompDens[ic] += dPhaseVolFrac[ip][Deriv::dC+ic] * phaseMassDens[ip]
                                        + phaseVolFrac[ip] * dMassDens_dC[ic];
        dTotalMassDens[Deriv::dC+ic] += dPhaseVolFrac[ip][Deriv::dC+ic] * phaseMassDens[ip]
                                        + phaseVolFrac[ip] * dMassDens_dC[ic];
      }

      totalMassDensityKernelOp( ip ); //, phaseVolFrac, dTotalMassDens_dPres, dTotalMassDens_dCompDens );
    }

  }

protected:

  // inputs

  /// Views on phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dPhaseVolFrac;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on phase mass densities
  arrayView3d< real64 const, multifluid::USD_PHASE > m_phaseMassDens;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > m_dPhaseMassDens;

  // outputs

  /// Views on total mass densities
  arrayView1d< real64 > m_totalMassDens;
  arrayView2d< real64, compflow::USD_FLUID_DC > m_dTotalMassDens;
  arrayView1d< real64 > m_dTotalMassDens_dPres;
  arrayView2d< real64, compflow::USD_FLUID_DC > m_dTotalMassDens_dCompDens;

};

/**
 * @class TotalMassDensityKernelFactory
 */
class TotalMassDensityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp the number of fluid components
   * @param[in] numPhase the number of fluid phases
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numPhase,
                   ObjectManagerBase & subRegion,
                   MultiFluidBase const & fluid )
  {
    if( numPhase == 2 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 2 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 2 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
    else if( numPhase == 3 )
    {
      isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        integer constexpr NUM_COMP = NC();
        TotalMassDensityKernel< NUM_COMP, 3 > kernel( subRegion, fluid );
        TotalMassDensityKernel< NUM_COMP, 3 >::template launch< POLICY >( subRegion.size(), kernel );
      } );
    }
  }
};


/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = solverBaseKernels::ResidualNormKernelBase< 1 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      integer const numComp,
                      integer const numDof,
                      integer const targetPhaseIndex,
                      WellElementSubRegion const & subRegion,
                      MultiFluidBase const & fluid,
                      WellControls const & wellControls,
                      real64 const timeAtEndOfStep,
                      real64 const dt,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_numComp( numComp ),
    m_numDof( numDof ),
    m_targetPhaseIndex( targetPhaseIndex ),
    m_dt( dt ),
    m_isLocallyOwned( subRegion.isLocallyOwned() ),
    m_iwelemControl( subRegion.getTopWellElementIndex() ),
    m_isProducer( wellControls.isProducer() ),
    m_currentControl( wellControls.getControl() ),
    m_targetBHP( wellControls.getTargetBHP( timeAtEndOfStep ) ),
    m_targetTotalRate( wellControls.getTargetTotalRate( timeAtEndOfStep ) ),
    m_targetPhaseRate( wellControls.getTargetPhaseRate( timeAtEndOfStep ) ),
    m_targetMassRate( wellControls.getTargetMassRate( timeAtEndOfStep ) ),
    m_volume( subRegion.getElementVolume() ),
    m_phaseDens_n( fluid.phaseDensity_n() ),
    m_totalDens_n( fluid.totalDensity_n() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const iwelem,
                            LinfStackVariables & stack ) const override
  {
    using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;

    real64 normalizer = 0.0;
    for( integer idof = 0; idof < m_numDof; ++idof )
    {

      // Step 1: compute a normalizer for the control or pressure equation

      // for the control equation, we distinguish two cases
      if( idof == ROFFSET::CONTROL )
      {

        // for the top well element, normalize using the current control
        if( m_isLocallyOwned && iwelem == m_iwelemControl )
        {
          if( m_currentControl == WellControls::Control::BHP )
          {
            // the residual entry is in pressure units
            normalizer = m_targetBHP;
          }
          else if( m_currentControl == WellControls::Control::TOTALVOLRATE )
          {
            // the residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetTotalRate ), m_minNormalizer );
          }
          else if( m_currentControl == WellControls::Control::PHASEVOLRATE )
          {
            // the residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetPhaseRate ), m_minNormalizer );
          }
          else if( m_currentControl == WellControls::Control::MASSRATE )
          {
            // the residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetMassRate ), m_minNormalizer );
          }
        }
        // for the pressure difference equation, always normalize by the BHP
        else
        {
          normalizer = m_targetBHP;
        }
      }
      // Step 2: compute a normalizer for the mass balance equations
      else if( idof >= ROFFSET::MASSBAL && idof < ROFFSET::MASSBAL + m_numComp )
      {
        if( m_isProducer ) // only PHASEVOLRATE is supported for now
        {
          // the residual is in mass units
          normalizer = m_dt * LvArray::math::abs( m_targetPhaseRate ) * m_phaseDens_n[iwelem][0][m_targetPhaseIndex];
        }
        else // Type::INJECTOR, only TOTALVOLRATE is supported for now
        {
          if( m_currentControl == WellControls::Control::MASSRATE )
          {
            normalizer = m_dt * LvArray::math::abs( m_targetMassRate );
          }
          else
          {
            // the residual is in mass units
            normalizer = m_dt * LvArray::math::abs( m_targetTotalRate ) * m_totalDens_n[iwelem][0];
          }

        }

        // to make sure that everything still works well if the rate is zero, we add this check
        normalizer = LvArray::math::max( normalizer, m_volume[iwelem] * m_totalDens_n[iwelem][0] );
      }
      // Step 3: compute a normalizer for the volume balance equations
      else
      {
        if( m_isProducer ) // only PHASEVOLRATE is supported for now
        {
          // the residual is in volume units
          normalizer = m_dt * LvArray::math::abs( m_targetPhaseRate );
        }
        else // Type::INJECTOR, only TOTALVOLRATE is supported for now
        {
          if( m_currentControl == WellControls::Control::MASSRATE )
          {
            normalizer = m_dt * LvArray::math::abs( m_targetMassRate/  m_totalDens_n[iwelem][0] );
          }
          else
          {
            normalizer = m_dt * LvArray::math::abs( m_targetTotalRate );
          }

        }

        // to make sure that everything still works well if the rate is zero, we add this check
        normalizer = LvArray::math::max( normalizer, m_volume[iwelem] );
      }
      normalizer = LvArray::math::max( m_minNormalizer, normalizer );

      // Step 4: compute the contribution to the residual
      real64 const val = LvArray::math::abs( m_localResidual[stack.localRow + idof] ) / normalizer;
      if( val > stack.localValue[0] )
      {
        stack.localValue[0] = val;
      }
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const iwelem,
                          L2StackVariables & stack ) const override
  {
    GEOS_UNUSED_VAR( iwelem, stack );
    GEOS_ERROR( "The L2 norm is not implemented for CompositionalMultiphaseWell" );
  }


protected:

  /// Number of fluid components
  integer const m_numComp;

  /// Number of dof per well element
  integer const m_numDof;

  /// Index of the target phase
  integer const m_targetPhaseIndex;

  /// Time step size
  real64 const m_dt;

  /// Flag indicating whether the well is locally owned or not
  bool const m_isLocallyOwned;

  /// Index of the element where the control is enforced
  localIndex const m_iwelemControl;

  /// Flag indicating whether the well is a producer or an injector
  bool const m_isProducer;

  /// Controls
  WellControls::Control const m_currentControl;
  real64 const m_targetBHP;
  real64 const m_targetTotalRate;
  real64 const m_targetPhaseRate;
  real64 const m_targetMassRate;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on phase/total density at the previous converged time step
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseDens_n;
  arrayView2d< real64 const, multifluid::USD_FLUID > const m_totalDens_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComp number of fluid components
   * @param[in] numDof number of dofs per well element
   * @param[in] targetPhaseIndex the index of the target phase (for phase volume control)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the well element subregion
   * @param[in] fluid the fluid model
   * @param[in] wellControls the controls
   * @param[in] timeAtEndOfStep the time at the end of the step (time_n + dt)
   * @param[in] dt the time step size
   * @param[out] residualNorm the residual norm on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComp,
                   integer const numDof,
                   integer const targetPhaseIndex,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   WellElementSubRegion const & subRegion,
                   MultiFluidBase const & fluid,
                   WellControls const & wellControls,
                   real64 const timeAtEndOfStep,
                   real64 const dt,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               numComp, numDof, targetPhaseIndex, subRegion, fluid, wellControls, timeAtEndOfStep, dt, minNormalizer );
    ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
  }

};

/******************************** ScalingForSystemSolutionKernel ********************************/

/**
 * @class ScalingForSystemSolutionKernelFactory
 */
class ScalingForSystemSolutionKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] maxRelativePresChange the max allowed relative pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   */
  template< typename POLICY >
  static isothermalCompositionalMultiphaseBaseKernels::ScalingForSystemSolutionKernel::StackVariables
  createAndLaunch( real64 const maxRelativePresChange,
                   real64 const maxAbsolutePresChange,
                   real64 const maxCompFracChange,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    arrayView1d< real64 const > const pressure =
      subRegion.getField< fields::well::pressure >();
    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getField< fields::well::globalCompDensity >();
    arrayView1d< real64 > pressureScalingFactor =
      subRegion.getField< fields::well::pressureScalingFactor >();
    arrayView1d< real64 > compDensScalingFactor =
      subRegion.getField< fields::well::globalCompDensityScalingFactor >();
    isothermalCompositionalMultiphaseBaseKernels::
      ScalingForSystemSolutionKernel kernel( maxRelativePresChange, maxAbsolutePresChange, maxCompFracChange, rankOffset,
                                             numComp, dofKey, subRegion, localSolution, pressure, compDens, pressureScalingFactor, compDensScalingFactor );
    return isothermalCompositionalMultiphaseBaseKernels::
             ScalingForSystemSolutionKernel::
             launch< POLICY >( subRegion.size(), kernel );
  }

};

/******************************** SolutionCheckKernel ********************************/

/**
 * @class SolutionCheckKernelFactory
 */
class SolutionCheckKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] allowCompDensChopping flag to allow the component density chopping
   * @param[in] scalingFactor the scaling factor
   * @param[in] rankOffset the rank offset
   * @param[in] numComp the number of components
   * @param[in] dofKey the dof key to get dof numbers
   * @param[in] subRegion the subRegion
   * @param[in] localSolution the Newton update
   */
  template< typename POLICY >
  static isothermalCompositionalMultiphaseBaseKernels::SolutionCheckKernel::StackVariables
  createAndLaunch( integer const allowCompDensChopping,
                   CompositionalMultiphaseFVM::ScalingType const scalingType,
                   real64 const scalingFactor,
                   globalIndex const rankOffset,
                   integer const numComp,
                   string const dofKey,
                   ElementSubRegionBase & subRegion,
                   arrayView1d< real64 const > const localSolution )
  {
    arrayView1d< real64 const > const pressure = subRegion.getField< fields::well::pressure >();
    arrayView2d< real64 const, compflow::USD_COMP > const compDens = subRegion.getField< fields::well::globalCompDensity >();
    arrayView1d< real64 > pressureScalingFactor = subRegion.getField< fields::well::pressureScalingFactor >();
    arrayView1d< real64 > compDensScalingFactor = subRegion.getField< fields::well::globalCompDensityScalingFactor >();
    isothermalCompositionalMultiphaseBaseKernels::
      SolutionCheckKernel kernel( allowCompDensChopping, 0, scalingType, scalingFactor, rankOffset, // no negative pressure
                                  numComp, dofKey, subRegion, localSolution, pressure, compDens, pressureScalingFactor, compDensScalingFactor );
    return isothermalCompositionalMultiphaseBaseKernels::
             SolutionCheckKernel::
             launch< POLICY >( subRegion.size(), kernel );
  }

};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of accumulation and volume balance
 */
template< integer NUM_COMP, integer NUM_DOF >
class ElementBasedAssemblyKernel
{
public:
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;
  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_DOF;

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( localIndex const numPhases,
                              globalIndex const rankOffset,
                              string const dofKey,
                              ElementSubRegionBase const & subRegion,
                              MultiFluidBase const & fluid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > const kernelFlags )
    : m_numPhases( numPhases ),
    m_rankOffset( rankOffset ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_volume( subRegion.getElementVolume() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::flow::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseVolFrac_n( subRegion.getField< fields::flow::phaseVolumeFraction_n >() ),
    m_phaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::flow::dPhaseVolumeFraction >() ),
    m_phaseDens_n( fluid.phaseDensity_n() ),
    m_phaseDens( fluid.phaseDensity() ),
    m_dPhaseDens( fluid.dPhaseDensity() ),
    m_phaseCompFrac_n( fluid.phaseCompFraction_n() ),
    m_phaseCompFrac( fluid.phaseCompFraction() ),
    m_dPhaseCompFrac( fluid.dPhaseCompFraction() ),
    m_compDens( subRegion.getField< fields::flow::globalCompDensity >() ),
    m_compDens_n( subRegion.getField< fields::flow::globalCompDensity_n >() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs ),
    m_kernelFlags( kernelFlags )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    //  volume information (used by both accumulation and volume balance)
    real64 volume = 0.0;

    // Residual information

    /// Index of the local row corresponding to this element
    localIndex localRow = -1;

    /// Indices of the matrix rows/columns corresponding to the dofs in this element
     globalIndex dofIndices[numDof]{};
    globalIndex eqnRowIndices[numComp+1]{};
    globalIndex dofColIndices[numComp+1]{};

    /// C-array storage for the element local residual vector (all equations except volume balance)
    real64 localResidual[numEqn]{};

    /// C-array storage for the element local Jacobian matrix (all equations except volume balance, all dofs)
    real64 localJacobian[numEqn][numDof]{};

  };

  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] ei the element index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const ei,
              StackVariables & stack ) const
  {
    // initialize the volume
    stack.volume = m_volume[ei];
    // set row index and degrees of freedom indices for this element (mass + vol bal)
    for( integer ic = 0; ic < numComp+1; ++ic )
    {
      stack.eqnRowIndices[ic] = m_dofNumber[ei] + ROFFSET::MASSBAL + ic - m_rankOffset;
    }

    // set DOF col indices for this block ( mass + vol bal)
    for( integer idof = 0; idof < numComp+1; ++idof )
    {
      stack.dofColIndices[idof] = m_dofNumber[ei] + COFFSET::DPRES + idof;
    }
    if ( 0)
    for ( integer jc = 0; jc < numEqn; ++jc )
    {
       stack.localResidual[jc] = 0.0;
       for ( integer ic = 0; ic < numDof; ++ic )
       {
        stack.localJacobian[jc][ic] = 0.0;
       }

    }

  }

  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] phaseAmountKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const ei,
                            StackVariables & stack,
                            FUNC && phaseAmountKernelOp = NoOpFunc{} ) const
  {

    using Deriv = multifluid::DerivativeOffset;

    // construct the slices for variables accessed multiple times
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac_n = m_phaseVolFrac_n[ei];
    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];

    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens_n = m_phaseDens_n[ei][0];
    arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > phaseDens = m_phaseDens[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > dPhaseDens = m_dPhaseDens[ei][0];

    arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac_n = m_phaseCompFrac_n[ei][0];
    arraySlice2d< real64 const, multifluid::USD_PHASE_COMP - 2 > phaseCompFrac = m_phaseCompFrac[ei][0];
    arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC - 2 > dPhaseCompFrac = m_dPhaseCompFrac[ei][0];

    // temporary work arrays
    real64 dPhaseAmount_dC[numComp]{};
    real64 dPhaseCompFrac_dC[numComp]{};

    // sum contributions to component accumulation from each phase
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      real64 const phaseAmount = stack.volume * phaseVolFrac[ip] * phaseDens[ip];
      real64 const phaseAmount_n = stack.volume * phaseVolFrac_n[ip] * phaseDens_n[ip];

      real64 const dPhaseAmount_dP = stack.volume * ( dPhaseVolFrac[ip][Deriv::dP] * phaseDens[ip]
                                                      + phaseVolFrac[ip] * dPhaseDens[ip][Deriv::dP] );

      // assemble density dependence
      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseDens[ip], dPhaseAmount_dC, Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                              + phaseDens[ip] * dPhaseVolFrac[ip][Deriv::dC+jc];
        dPhaseAmount_dC[jc] *= stack.volume;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( integer ic = 0; ic < numComp; ++ic )
      {
        real64 const phaseCompAmount = phaseAmount * phaseCompFrac[ip][ic];
        real64 const phaseCompAmount_n = phaseAmount_n * phaseCompFrac_n[ip][ic];

        real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                           + phaseAmount * dPhaseCompFrac[ip][ic][Deriv::dP];

        stack.localResidual[ic] += phaseCompAmount - phaseCompAmount_n;
        stack.localJacobian[ic][0] += dPhaseCompAmount_dP;

        // jc - index of component w.r.t. whose compositional var the derivative is being taken
        // (i.e. col number in local matrix)

        // assemble phase composition dependence
        applyChainRule( numComp, dCompFrac_dCompDens, dPhaseCompFrac[ip][ic], dPhaseCompFrac_dC, Deriv::dC );
        for( integer jc = 0; jc < numComp; ++jc )
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmount
                                             + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];

          stack.localJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
        }
      }

      // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
      // possible use: assemble the derivatives wrt temperature, and the accumulation term of the energy equation for this phase
      phaseAmountKernelOp( ip, phaseAmount, phaseAmount_n, dPhaseAmount_dP, dPhaseAmount_dC );

    }

    // check zero diagonal (works only in debug)
    /*
       for( integer ic = 0; ic < numComp; ++ic )
       {
       GEOS_ASSERT_MSG ( LvArray::math::abs( stack.localJacobian[ic][ic] ) > minDensForDivision,
                        GEOS_FMT( "Zero diagonal in Jacobian: equation {}, value = {}", ic, stack.localJacobian[ic][ic] ) );
       }
     */
  }


  /**
   * @brief Compute the local volume balance contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] phaseVolFractionSumKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeVolumeBalance( localIndex const ei,
                             StackVariables & stack,
                             FUNC && phaseVolFractionSumKernelOp = NoOpFunc{} ) const
  {
    using Deriv = multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];

    real64 oneMinusPhaseVolFracSum = 1.0;

    // sum contributions to component accumulation from each phase
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      oneMinusPhaseVolFracSum -= phaseVolFrac[ip];
      stack.localJacobian[numComp][0] -= dPhaseVolFrac[ip][Deriv::dP];

      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localJacobian[numComp][jc+1] -= dPhaseVolFrac[ip][Deriv::dC+jc];
      }
    }

    // call the lambda in the phase loop to allow the reuse of the phase amounts and their derivatives
    // possible use: assemble the derivatives wrt temperature, and use oneMinusPhaseVolFracSum if poreVolume depends on temperature
    // tjb revisit
    phaseVolFractionSumKernelOp( oneMinusPhaseVolFracSum );

    // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
    stack.localResidual[numComp] = stack.volume * oneMinusPhaseVolFracSum;
    for( integer idof = 0; idof < numComp+1; ++idof )
    {
      stack.localJacobian[numComp][idof] *= stack.volume;
    }

  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void complete( localIndex const GEOS_UNUSED_PARAM( ei ),
                 StackVariables & stack ) const
  {
    using namespace compositionalMultiphaseUtilities;


    if( m_kernelFlags.isSet( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation ) )
    {
      // apply equation/variable change transformation to the component mass balance equations
      real64 work[numComp + 1]{};
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numComp+1, stack.localJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.localResidual );
    }

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // - the volume balance equations (i = numComp)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    integer const numRows = numComp+1;
    for( integer i = 0; i < numRows; ++i )
    {
      m_localRhs[stack.eqnRowIndices[i]]  += stack.localResidual[i];
      m_localMatrix.addToRow< serialAtomic >( stack.eqnRowIndices[i],
                                              stack.dofColIndices,
                                              stack.localJacobian[i],
                                              numComp+1 );
    }

  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( kernelComponent.elemGhostRank( ei ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( ei, stack );
      kernelComponent.computeAccumulation( ei, stack );
      kernelComponent.computeVolumeBalance( ei, stack );
      kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Number of fluid phases
  integer const m_numPhases;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_volume;

  /// Views on the porosity
  arrayView2d< real64 const > const m_porosity_n;
  arrayView2d< real64 const > const m_porosity;
  arrayView2d< real64 const > const m_dPoro_dPres;

  /// Views on the derivatives of comp fractions wrt component density
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dCompFrac_dCompDens;

  /// Views on the phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFrac_n;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_phaseVolFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > const m_dPhaseVolFrac;

  /// Views on the phase densities
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseDens_n;
  arrayView3d< real64 const, multifluid::USD_PHASE > const m_phaseDens;
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const m_dPhaseDens;

  /// Views on the phase component fraction
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const m_phaseCompFrac_n;
  arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const m_phaseCompFrac;
  arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const m_dPhaseCompFrac;

  // Views on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens_n;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > const m_kernelFlags;
};


/**
 * @class ElementBasedAssemblyKernelFactory
 */
class ElementBasedAssemblyKernelFactory
{
public:
  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( localIndex const numComps,
                   localIndex const numPhases,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const dofKey,
                   ElementSubRegionBase const & subRegion,
                   MultiFluidBase const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      localIndex constexpr NUM_COMP = NC();
      localIndex constexpr NUM_DOF = NC()+2;

      BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation );

      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >
      kernel( numPhases, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs, kernelFlags );
      ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF >::template
      launch< POLICY, ElementBasedAssemblyKernel< NUM_COMP, NUM_DOF > >( subRegion.size(), kernel );
    } );
  }
};
/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NC, integer NUM_DOF >
class FaceBasedAssemblyKernel
{
public:

  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;
  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using TAG = compositionalMultiphaseWellKernels::ElemTag;


  /// Compile time value for the number of components
  static constexpr integer numComp = NC;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;


  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] capPressureAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  FaceBasedAssemblyKernel( real64 const dt,
                           globalIndex const rankOffset,
                           string const wellDofKey,
                           WellControls const & wellControls,
                           ElementSubRegionBase const & subRegion,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > kernelFlags )
    :
    m_dt( dt ),
    m_rankOffset( rankOffset ),
    m_wellElemDofNumber ( subRegion.getReference< array1d< globalIndex > >( wellDofKey ) ),
    m_nextWellElemIndex ( subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString()) ),
    m_connRate ( subRegion.getField< fields::well::mixtureConnectionRate >() ),
    m_wellElemCompFrac ( subRegion.getField< fields::well::globalCompFraction >() ),
    m_dWellElemCompFrac_dCompDens ( subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_localMatrix( localMatrix ),
    m_localRhs ( localRhs ),
    m_useTotalMassEquation ( kernelFlags.isSet( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation ) ),
    m_isProducer ( wellControls.isProducer() ),
    m_injection ( wellControls.getInjectionStream() )
  { }


  GEOS_HOST_DEVICE
  inline
  void
  computeExit( real64 const & dt,
               real64 const ( &compFlux )[NC],
               real64 const ( &dCompFlux_dRate )[NC],
               real64 const ( &dCompFlux_dPresUp )[NC],
               real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
               real64 ( & oneSidedFlux )[NC],
               real64 ( & oneSidedFluxJacobian_dRate )[NC][1],
               real64 ( & oneSidedFluxJacobian_dPresCompUp )[NC][NC + 1] ) const
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

  GEOS_HOST_DEVICE
  inline
  void
  compute( real64 const & dt,
           real64 const ( &compFlux )[NC],
           real64 const ( &dCompFlux_dRate )[NC],
           real64 const ( &dCompFlux_dPresUp )[NC],
           real64 const ( &dCompFlux_dCompDensUp )[NC][NC],
           real64 ( & localFlux )[2*NC],
           real64 ( & localFluxJacobian_dRate )[2*NC][1],
           real64 ( & localFluxJacobian_dPresCompUp )[2*NC][NC + 1] ) const
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

    using namespace compositionalMultiphaseUtilities;

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

    localIndex const iwelemNext = m_nextWellElemIndex[iwelem];
    real64 const currentConnRate = m_connRate[iwelem];
    localIndex iwelemUp = -1;

    if( iwelemNext < 0 && !m_isProducer )  // exit connection, injector
    {
      // we still need to define iwelemUp for Jacobian assembly
      iwelemUp = iwelem;

      // just copy the injection stream into compFrac
      for( integer ic = 0; ic < NC; ++ic )
      {
        compFracUp[ic] = m_injection[ic];
        for( integer jc = 0; jc < NC; ++jc )
        {
          dCompFrac_dCompDensUp[ic][jc] = 0.0;
        }
      }
    }
    else
    {
      // first set iwelemUp to the upstream cell
      if( ( iwelemNext < 0 && m_isProducer )  // exit connection, producer
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
        compFracUp[ic] = m_wellElemCompFrac[iwelemUp][ic];
        for( integer jc = 0; jc < NC; ++jc )
        {
          dCompFrac_dCompDensUp[ic][jc] = m_dWellElemCompFrac_dCompDens[iwelemUp][ic][jc];
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

    globalIndex const offsetUp = m_wellElemDofNumber[iwelemUp];
    globalIndex const offsetCurrent = m_wellElemDofNumber[iwelem];

    if( iwelemNext < 0 )  // exit connection
    {
      // for this case, we only need NC mass conservation equations
      // so we do not use the arrays initialized before the loop
      real64 oneSidedFlux[NC]{};
      real64 oneSidedFluxJacobian_dRate[NC][1]{};
      real64 oneSidedFluxJacobian_dPresCompUp[NC][NC+1]{};

      computeExit ( m_dt,
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
        oneSidedEqnRowIndices[ic] = offsetUp + ROFFSET::MASSBAL + ic - m_rankOffset;
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

      if( m_useTotalMassEquation > 0 )
      {
        // Apply equation/variable change transformation(s)
        real64 work[NC + 1]{};
        shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, 1, oneSidedFluxJacobian_dRate, work );
        shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC + 1, oneSidedFluxJacobian_dPresCompUp, work );
        shiftElementsAheadByOneAndReplaceFirstElementWithSum( NC, oneSidedFlux );
      }

      for( integer i = 0; i < NC; ++i )
      {
        if( oneSidedEqnRowIndices[i] >= 0 && oneSidedEqnRowIndices[i] < m_localMatrix.numRows() )
        {
          m_localMatrix.addToRow< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                          &oneSidedDofColIndices_dRate,
                                                          oneSidedFluxJacobian_dRate[i],
                                                          1 );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( oneSidedEqnRowIndices[i],
                                                                              oneSidedDofColIndices_dPresCompUp,
                                                                              oneSidedFluxJacobian_dPresCompUp[i],
                                                                              NC+1 );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[oneSidedEqnRowIndices[i]], oneSidedFlux[i] );
        }
      }
    }
    else // not an exit connection
    {
      real64 localFlux[2*NC]{};
      real64 localFluxJacobian_dRate[2*NC][1]{};
      real64 localFluxJacobian_dPresCompUp[2*NC][NC+1]{};

      compute( m_dt,
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

      globalIndex const offsetNext = m_wellElemDofNumber[iwelemNext];

      // jacobian indices
      for( integer ic = 0; ic < NC; ++ic )
      {
        // mass balance equations for all components
        eqnRowIndices[TAG::NEXT *NC+ic]    = offsetNext + ROFFSET::MASSBAL + ic - m_rankOffset;
        eqnRowIndices[TAG::CURRENT *NC+ic] = offsetCurrent + ROFFSET::MASSBAL + ic - m_rankOffset;
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

      if( m_useTotalMassEquation > 0 )
      {
        // Apply equation/variable change transformation(s)
        real64 work[NC + 1]{};
        shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC, 1, 2, localFluxJacobian_dRate, work );
        shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, NC, NC + 1, 2, localFluxJacobian_dPresCompUp, work );
        shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( NC, NC, 2, localFlux );
      }

      for( integer i = 0; i < 2*NC; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < m_localMatrix.numRows() )
        {
          m_localMatrix.addToRow< parallelDeviceAtomic >( eqnRowIndices[i],
                                                          &dofColIndices_dRate,
                                                          localFluxJacobian_dRate[i],
                                                          1 );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                              dofColIndices_dPresCompUp,
                                                                              localFluxJacobian_dPresCompUp[i],
                                                                              NC+1 );
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
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElements,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numElements, [=] GEOS_HOST_DEVICE ( localIndex const ie )
    {
      //typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
      //                                            kernelComponent.numPointsInFlux( iconn ) );

      //kernelComponent.setup( iconn, stack );
      kernelComponent.computeFlux( ie );
      //kernelComponent.complete( iconn, stack );
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
  /// Element component fraction derivatives
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dWellElemCompFrac_dCompDens;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

  /// Kernel option flag
  integer const m_useTotalMassEquation;

  /// Well type
  bool const m_isProducer;

  /// Injection stream composition
  arrayView1d< real64 const > const m_injection;


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
   * @param[in] numComps the number of fluid components
   * @param[in] dt time step size
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] useTotalMassEquation flag specifying whether to replace one component bal eqn with total mass eqn
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
                   integer const useTotalMassEquation,
                   string const dofKey,
                   WellControls const & wellControls,
                   ElementSubRegionBase const & subRegion,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 2;

      BitFlags< isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::ElementBasedAssemblyKernelFlags::TotalMassEquation );


      using kernelType = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF >;


      kernelType kernel( dt, rankOffset, dofKey, wellControls, subRegion, localMatrix, localRhs, kernelFlags );
      kernelType::template launch< POLICY >( subRegion.size(), kernel );
    } );
  }
};
} // end namespace compositionalMultiphaseWellKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP
