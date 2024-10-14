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

// define the row offset of the residual equations
struct RowOffset
{
  static constexpr integer CONTROL = 0;
  static constexpr integer MASSBAL = 1;
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

  template< integer NC >
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

  using TAG = compositionalMultiphaseWellKernels::ElemTag;
  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

  template< integer NC >
  GEOS_HOST_DEVICE
  inline
  static void
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
             real64 ( &localPresRelJacobian )[2*(NC+1)] );

  template< integer NC >
  static void
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
    StencilMaterialAccessors< constitutive::MultiFluidBase,
                              fields::multifluid::phaseDensity,
                              fields::multifluid::dPhaseDensity,
                              fields::multifluid::phaseViscosity,
                              fields::multifluid::dPhaseViscosity,
                              fields::multifluid::phaseCompFraction,
                              fields::multifluid::dPhaseCompFraction >;

  using RelPermAccessors =
    StencilMaterialAccessors< constitutive::RelativePermeabilityBase,
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


  template< integer NC, integer NP >
  GEOS_HOST_DEVICE
  inline
  static void
  compute( bool const & disableReservoirToWellFlow,
           real64 const & resPres,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & resPhaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dResPhaseVolFrac,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dResCompFrac_dCompDens,
           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & resPhaseDens,
           arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dResPhaseDens,
           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & resPhaseVisc,
           arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dResPhaseVisc,
           arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const & resPhaseCompFrac,
           arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 2 > const & dResPhaseCompFrac,
           arraySlice1d< real64 const, constitutive::relperm::USD_RELPERM - 2 > const & resPhaseRelPerm,
           arraySlice2d< real64 const, constitutive::relperm::USD_RELPERM_DS - 2 > const & dResPhaseRelPerm_dPhaseVolFrac,
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
           arraySlice3d< real64 > const & dCompPerfRate_dComp );

  template< integer NC, integer NP >
  static void
  launch( localIndex const size,
          bool const disableReservoirToWellFlow,
          ElementViewConst< arrayView1d< real64 const > > const & resPres,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & resPhaseVolFrac,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dResPhaseVolFrac_dComp,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dResCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & resPhaseDens,
          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dResPhaseDens,
          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & resPhaseVisc,
          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dResPhaseVisc,
          ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > > const & resPhaseCompFrac,
          ElementViewConst< arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > > const & dResPhaseCompFrac,
          ElementViewConst< arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > > const & resPhaseRelPerm,
          ElementViewConst< arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > > const & dResPhaseRelPerm_dPhaseVolFrac,
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
             arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & phaseDens,
             arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dPhaseDens,
             arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac,
             arraySlice3d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC - 2 > const & dPhaseCompFrac,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac_n,
             arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & phaseDens_n,
             arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP - 2 > const & phaseCompFrac_n,
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
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & wellElemPhaseDens,
          arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dWellElemPhaseDens,
          arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac,
          arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > const & dWellElemPhaseCompFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & wellElemPhaseVolFrac_n,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & wellElemPhaseDens_n,
          arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > const & wellElemPhaseCompFrac_n,
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
    StencilMaterialAccessors< constitutive::MultiFluidBase,
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
          ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & resPhaseMassDens,
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
          arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & wellElemTotalDens,
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
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens,
          arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const & totalDens,
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
                          constitutive::MultiFluidBase const & fluid )
    : Base(),
    m_phaseVolFrac( subRegion.getField< fields::well::phaseVolumeFraction >() ),
    m_dPhaseVolFrac( subRegion.getField< fields::well::dPhaseVolumeFraction >() ),
    m_dCompFrac_dCompDens( subRegion.getField< fields::well::dGlobalCompFraction_dGlobalCompDensity >() ),
    m_phaseMassDens( fluid.phaseMassDensity() ),
    m_dPhaseMassDens( fluid.dPhaseMassDensity() ),
    m_totalMassDens( subRegion.getField< fields::well::totalMassDensity >() ),
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
    using Deriv = constitutive::multifluid::DerivativeOffset;

    arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac = m_phaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > dPhaseVolFrac = m_dPhaseVolFrac[ei];
    arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens = m_dCompFrac_dCompDens[ei];
    arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseMassDens = m_phaseMassDens[ei][0];
    arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > dPhaseMassDens = m_dPhaseMassDens[ei][0];
    real64 & totalMassDens = m_totalMassDens[ei];
    real64 & dTotalMassDens_dPres = m_dTotalMassDens_dPres[ei];
    arraySlice1d< real64, compflow::USD_FLUID_DC - 1 > dTotalMassDens_dCompDens = m_dTotalMassDens_dCompDens[ei];

    real64 dMassDens_dC[numComp]{};

    totalMassDens = 0.0;
    dTotalMassDens_dPres = 0.0;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      dTotalMassDens_dCompDens[ic] = 0.0;
    }

    for( integer ip = 0; ip < numPhase; ++ip )
    {
      totalMassDens += phaseVolFrac[ip] * phaseMassDens[ip];
      dTotalMassDens_dPres += dPhaseVolFrac[ip][Deriv::dP] * phaseMassDens[ip] + phaseVolFrac[ip] * dPhaseMassDens[ip][Deriv::dP];

      applyChainRule( numComp, dCompFrac_dCompDens, dPhaseMassDens[ip], dMassDens_dC, Deriv::dC );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        dTotalMassDens_dCompDens[ic] += dPhaseVolFrac[ip][Deriv::dC+ic] * phaseMassDens[ip]
                                        + phaseVolFrac[ip] * dMassDens_dC[ic];
      }

      totalMassDensityKernelOp( ip, totalMassDens, dTotalMassDens_dPres, dTotalMassDens_dCompDens );
    }

  }

protected:

  // inputs

  /// Views on phase volume fractions
  arrayView2d< real64 const, compflow::USD_PHASE > m_phaseVolFrac;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > m_dPhaseVolFrac;
  arrayView3d< real64 const, compflow::USD_COMP_DC > m_dCompFrac_dCompDens;

  /// Views on phase mass densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_phaseMassDens;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dPhaseMassDens;

  // outputs

  /// Views on total mass densities
  arrayView1d< real64 > m_totalMassDens;
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
                   constitutive::MultiFluidBase const & fluid )
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
                      constitutive::MultiFluidBase const & fluid,
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
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_phaseDens_n;
  arrayView2d< real64 const, constitutive::multifluid::USD_FLUID > const m_totalDens_n;

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
                   constitutive::MultiFluidBase const & fluid,
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
   * @param[in] maxAbsolutePresChange the max allowed absolute pressure change
   * @param[in] maxCompFracChange the max allowed comp fraction change
   * @param[in] maxRelativeCompDensChange the max allowed relative comp density change
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
                   real64 const maxRelativeCompDensChange,
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
      ScalingForSystemSolutionKernel kernel( maxRelativePresChange, maxAbsolutePresChange, maxCompFracChange, maxRelativeCompDensChange, rankOffset,
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


} // end namespace compositionalMultiphaseWellKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELLKERNELS_HPP
