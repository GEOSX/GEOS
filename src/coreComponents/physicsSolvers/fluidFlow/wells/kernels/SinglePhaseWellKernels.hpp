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
 * @file SinglePhaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "physicsSolvers/KernelLaunchSelectors.hpp"

namespace geos
{

namespace singlePhaseWellKernels
{

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
  static constexpr integer DRATE = 1;
};

template< integer IS_THERMAL >
struct ColOffset_WellJac;

template<>
struct ColOffset_WellJac< 0 >
{
  static constexpr integer dP = 0;
  static constexpr integer dQ = dP + 1;
  static integer constexpr nDer =  dQ + 1;

};

template<>
struct ColOffset_WellJac< 1 >
{
  static constexpr integer dP = 0;
  static constexpr integer dQ = dP + 1;
  static constexpr integer dT = dQ+1;
/// number of derivatives
  static integer constexpr nDer =  dT + 1;
};

// define the row offset of the residual equations
struct RowOffset
{
  static constexpr integer CONTROL = 0;
  static constexpr integer MASSBAL = 1;
};

template< integer IS_THERMAL >
struct RowOffset_WellJac;

template<>
struct RowOffset_WellJac< 0 >
{
  static constexpr integer CONTROL   = 0;
  static constexpr integer MASSBAL   = 1;
  static constexpr integer nEqn      = MASSBAL+1;
};
template<>
struct RowOffset_WellJac< 1 >
{
  static constexpr integer CONTROL   = 0;
  static constexpr integer MASSBAL   = 1;
  static constexpr integer ENERGYBAL = MASSBAL+1;
  static constexpr integer nEqn      = ENERGYBAL+1;

};
/******************************** ControlEquationHelper ********************************/

struct ControlEquationHelper
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  // add an epsilon to the checks to avoid control changes due to tiny pressure/rate updates
  static constexpr real64 EPS = 1e-15;

  GEOS_HOST_DEVICE
  inline
  static
  void
  switchControl( bool const isProducer,
                 WellControls::Control const & currentControl,
                 real64 const & targetBHP,
                 real64 const & targetRate,
                 real64 const & currentBHP,
                 real64 const & currentVolRate,
                 WellControls::Control & newControl );

  template< integer IS_THERMAL >
  GEOS_HOST_DEVICE
  inline
  static
  void
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
           arrayView1d< real64 > const & localRhs );

};


/******************************** FluxKernel ********************************/

struct FluxKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;
  using TAG = singlePhaseWellKernels::ElemTag;

  template< integer IS_THERMAL >
  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};


/******************************** PressureRelationKernel ********************************/

struct PressureRelationKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;
  using TAG = singlePhaseWellKernels::ElemTag;

  template< integer IS_THERMAL >
  static localIndex
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
          arrayView1d< real64 > const & localRhs );

};


/******************************** PerforationKernel ********************************/

struct PerforationKernel
{

  using TAG = singlePhaseWellKernels::SubRegionTag;

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::pressure >;

  using SingleFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity,
                              fields::singlefluid::viscosity,
                              fields::singlefluid::dViscosity >;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< integer IS_THERMAL >
  GEOS_HOST_DEVICE
  inline
  static
  void
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
           arraySlice1d< real64 > const & dPerfRate_dPres );

  template< integer IS_THERMAL >
  static void
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
          arrayView2d< real64 > const & dPerfRate_dPres );

};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  using ROFFSET = singlePhaseWellKernels::RowOffset;
  using COFFSET = singlePhaseWellKernels::ColOffset;

  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< integer const > const & wellElemGhostRank,
          arrayView1d< real64 const > const & wellElemVolume,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity,
          arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dWellElemDensity,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDensity_n,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};

/******************************** ElementBasedAssemblyKernel ********************************/

/**
 * @class ElementBasedAssemblyKernel
 * @tparam IS_THERMAL thermal switch
 * @brief Define the interface for the assembly kernel in charge of accumulation
 */
template< integer IS_THERMAL >
class ElementBasedAssemblyKernel
{
public:
  using COFFSET = singlePhaseWellKernels::ColOffset;
  using ROFFSET = singlePhaseWellKernels::RowOffset;

  // Well jacobian column and row indicies
  using FLUID_PROP_COFFSET = constitutive::singlefluid::DerivativeOffsetC< IS_THERMAL >;
  using WJ_COFFSET = singlePhaseWellKernels::ColOffset_WellJac< IS_THERMAL >;
  using WJ_ROFFSET = singlePhaseWellKernels::RowOffset_WellJac< IS_THERMAL >;

  /// Number of Dof's set in this kernal   - no dQ in accum
  static constexpr integer numDof = 1 + IS_THERMAL;  // tjb review

  /// Compute time value for the number of equations  mass bal + momentum + energy bal
  static constexpr integer numEqn =  1 + IS_THERMAL;  // tjb review


  /**
   * @brief Constructor
   * @param[in] isProducer well type
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  ElementBasedAssemblyKernel( integer const isProducer,
                              globalIndex const rankOffset,
                              string const dofKey,
                              WellElementSubRegion const & subRegion,
                              constitutive::SingleFluidBase const & fluid,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs )
    : m_isProducer( isProducer ),
    m_rankOffset( rankOffset ),
    m_iwelemControl( subRegion.getTopWellElementIndex() ),
    m_dofNumber( subRegion.getReference< array1d< globalIndex > >( dofKey ) ),
    m_elemGhostRank( subRegion.ghostRank() ),
    m_wellElemVolume( subRegion.getElementVolume() ),
    m_wellElemDensity( fluid.density() ),
    m_dWellElemDensity( fluid.dDensity() ),
    m_wellElemDensity_n( fluid.density_n() ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}


  /**
   * @brief Getter for the ghost rank of an element
   * @param[in] ei the element index
   * @return the ghost rank of the element
   */
  GEOS_HOST_DEVICE
  integer elemGhostRank( localIndex const ei ) const
  { return m_elemGhostRank( ei ); }



  /**
   * @brief Compute the local accumulation contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] ei the element index
   * @param[inout] stack the stack variables
   * @param[in] phaseAmountKernelOp the function used to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeAccumulation( localIndex const iwelem,
                            FUNC && kernelOp = NoOpFunc{} ) const
  {

    localIndex const eqnRowIndex = m_dofNumber[iwelem] + ROFFSET::MASSBAL - m_rankOffset;
    globalIndex const presDofColIndex = m_dofNumber[iwelem] + COFFSET::DPRES;

    real64 const localAccum = m_wellElemVolume[iwelem] * ( m_wellElemDensity[iwelem][0] - m_wellElemDensity_n[iwelem][0] );
    real64 const localAccumJacobian = m_wellElemVolume[iwelem] * m_dWellElemDensity[iwelem][0][FLUID_PROP_COFFSET::dP];

    // add contribution to global residual and jacobian (no need for atomics here)
    m_localMatrix.addToRow< serialAtomic >( eqnRowIndex, &presDofColIndex, &localAccumJacobian, 1 );
    m_localRhs[eqnRowIndex] += localAccum;

    //kernelOp

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

    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const iwelem )
    {
      if( kernelComponent.elemGhostRank( iwelem ) >= 0 )
      {
        return;
      }

      //typename KERNEL_TYPE::StackVariables stack;

      //kernelComponent.setup( ei, stack );
      kernelComponent.computeAccumulation( iwelem );
      //kernelComponent.computeVolumeBalance( ei, stack );
      //kernelComponent.complete( ei, stack );
    } );
  }

protected:

  /// Well type
  integer const m_isProducer;

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Index of the element where the control is enforced
  localIndex const m_iwelemControl;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_elemGhostRank;

  /// View on the element volumes
  arrayView1d< real64 const > const m_wellElemVolume;

  /// Views on the density
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_wellElemDensity;
  arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const m_dWellElemDensity;
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_wellElemDensity_n;

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;

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
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const isProducer,
                   globalIndex const rankOffset,
                   string const dofKey,
                   WellElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer isThermal=0;
    geos::internal::kernelLaunchSelectorThermalSwitch( isThermal, [&]( auto IS_THERMAL )
    {

      integer constexpr istherm = IS_THERMAL();

      ElementBasedAssemblyKernel< istherm >
      kernel( isProducer, rankOffset, dofKey, subRegion, fluid, localMatrix, localRhs );
      ElementBasedAssemblyKernel< istherm >::template
      launch< POLICY, ElementBasedAssemblyKernel< istherm > >( subRegion.size(), kernel );
    } );
  }
};

/******************************** PressureTemperatyrInitializationKernel ********************************/

// tjb make this templated on thermal
struct PresTempInitializationKernel
{

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::flow::pressure,
                      fields::flow::temperature >;

  using SingleFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density >;

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
  launch( integer const isThermal,
          localIndex const perforationSize,
          localIndex const subRegionSize,
          localIndex const numPerforations,
          WellControls const & wellControls,
          real64 const & currentTime,
          ElementViewConst< arrayView1d< real64 const > > const & resPressure,
          ElementViewConst< arrayView1d< real64 const > > const & resTemperature,
          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & resDensity,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 > const & wellElemPressure,
          arrayView1d< real64 > const & wellElemTemperature );

};

/******************************** RateInitializationKernel ********************************/

struct RateInitializationKernel
{

  static void
  launch( localIndex const subRegionSize,
          WellControls const & wellControls,
          real64 const & currentTime,
          arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & wellElemDens,
          arrayView1d< real64 > const & connRate );

};


/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public physicsSolverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = physicsSolverBaseKernels::ResidualNormKernelBase< 1 >;
  using Base::m_minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      WellElementSubRegion const & subRegion,
                      constitutive::SingleFluidBase const & fluid,
                      WellControls const & wellControls,
                      real64 const timeAtEndOfStep,
                      real64 const dt,
                      real64 const minNormalizer )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank,
            minNormalizer ),
    m_dt( dt ),
    m_isLocallyOwned( subRegion.isLocallyOwned() ),
    m_iwelemControl( subRegion.getTopWellElementIndex() ),
    m_currentControl( wellControls.getControl() ),
    m_targetBHP( wellControls.getTargetBHP( timeAtEndOfStep ) ),
    m_targetRate( wellControls.getTargetTotalRate( timeAtEndOfStep ) ),
    m_volume( subRegion.getElementVolume() ),
    m_density_n( fluid.density_n() )
  {}

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const iwelem,
                            LinfStackVariables & stack ) const override
  {
    for( localIndex idof = 0; idof < 2; ++idof )
    {
      real64 normalizer = 0.0;
      if( idof == singlePhaseWellKernels::RowOffset::CONTROL )
      {
        // for the top well element, normalize using the current control
        if( m_isLocallyOwned && iwelem == m_iwelemControl )
        {
          if( m_currentControl == WellControls::Control::BHP )
          {
            // this residual entry is in pressure units
            normalizer = m_targetBHP;
          }
          else if( m_currentControl == WellControls::Control::TOTALVOLRATE )
          {
            // this residual entry is in volume / time units
            normalizer = LvArray::math::max( LvArray::math::abs( m_targetRate ), m_minNormalizer );
          }
        }
        // for the pressure difference equation, always normalize by the BHP
        else
        {
          // this residual is in pressure units
          normalizer = m_targetBHP;
        }
      }
      else // SinglePhaseWell::RowOffset::MASSBAL
      {
        // this residual entry is in mass units
        normalizer = m_dt * LvArray::math::abs( m_targetRate ) * m_density_n[iwelem][0];

        // to make sure that everything still works well if the rate is zero, we add this check
        normalizer = LvArray::math::max( normalizer, m_volume[iwelem] * m_density_n[iwelem][0] );
      }

      // we have the normalizer now, we can compute a dimensionless Linfty norm contribution
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
    GEOS_ERROR( "The L2 norm is not implemented for SinglePhaseWell" );
  }


protected:

  /// Time step size
  real64 const m_dt;

  /// Flag indicating whether the well is locally owned or not
  bool const m_isLocallyOwned;

  /// Index of the element where the control is enforced
  localIndex const m_iwelemControl;

  /// Controls
  WellControls::Control const m_currentControl;
  real64 const m_targetBHP;
  real64 const m_targetRate;

  /// View on the volume
  arrayView1d< real64 const > const m_volume;

  /// View on total density at the previous converged time step
  arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const m_density_n;

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
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   WellElementSubRegion const & subRegion,
                   constitutive::SingleFluidBase const & fluid,
                   WellControls const & wellControls,
                   real64 const timeAtEndOfStep,
                   real64 const dt,
                   real64 const minNormalizer,
                   real64 (& residualNorm)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, subRegion, fluid, wellControls, timeAtEndOfStep, dt, minNormalizer );
    ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
  }

};

} // end namespace singlePhaseWellKernels

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELLKERNELS_HPP
