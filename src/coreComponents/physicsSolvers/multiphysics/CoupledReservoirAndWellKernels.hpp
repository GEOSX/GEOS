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
 * @file ThermalCompositionalMultiphaseWellKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLS_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLS_HPP


#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "physicsSolvers/fluidFlow/wells/kernels/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellFields.hpp"
#include "physicsSolvers/fluidFlow/wells/WellTags.hpp"
#include "physicsSolvers/fluidFlow/wells/WellFields.hpp"
namespace geos
{

namespace coupledReservoirAndWellKernels
{

using namespace constitutive;

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NC, integer IS_THERMAL >
class IsothermalCompositionalMultiPhaseFluxKernel
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NC;
  static constexpr integer resNumDOF  = NC+1+IS_THERMAL;

  // Well jacobian column and row indicies
  using WJ_COFFSET = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
  using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NC, IS_THERMAL >;

  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

  using CP_Deriv = multifluid::DerivativeOffsetC< NC, IS_THERMAL >;

  using TAG = compositionalMultiphaseWellKernels::SubRegionTag;



  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = WJ_COFFSET::nDer;

  /// Compile time value for the number of equations except volume and momentum
  static constexpr integer numEqn = WJ_ROFFSET::nEqn - 2;

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
  IsothermalCompositionalMultiPhaseFluxKernel( real64 const dt,
                                               globalIndex const rankOffset,
                                               string const wellDofKey,
                                               WellElementSubRegion const & subRegion,
                                               ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber,
                                               PerforationData const * const perforationData,
                                               MultiFluidBase const & fluid,

                                               arrayView1d< real64 > const & localRhs,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               bool const & detectCrossflow,
                                               integer & numCrossFlowPerforations,
                                               BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags )
    :
    m_dt( dt ),
    m_numPhases ( fluid.numFluidPhases()),
    m_rankOffset( rankOffset ),
    m_compPerfRate( perforationData->getField< fields::well::compPerforationRate >() ),
    m_dCompPerfRate( perforationData->getField< fields::well::dCompPerforationRate >() ),
    m_perfWellElemIndex( perforationData->getField< fields::perforation::wellElementIndex >() ),
    m_wellElemDofNumber( subRegion.getReference< array1d< globalIndex > >( wellDofKey ) ),
    m_resElemDofNumber( resDofNumber ),
    m_resElementRegion( perforationData->getField< fields::perforation::reservoirElementRegion >() ),
    m_resElementSubRegion( perforationData->getField< fields::perforation::reservoirElementSubRegion >() ),
    m_resElementIndex( perforationData->getField< fields::perforation::reservoirElementIndex >() ),
    m_localRhs( localRhs ),
    m_localMatrix( localMatrix ),
    m_detectCrossflow( detectCrossflow ),
    m_numCrossFlowPerforations( numCrossFlowPerforations ),
    m_useTotalMassEquation ( kernelFlags.isSet( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation ) )
  { }


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
  void computeFlux( localIndex const iperf,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {

    using namespace compositionalMultiphaseUtilities;
    // local working variables and arrays
    stackArray1d< localIndex, 2* numComp > eqnRowIndices( 2 * numComp );
    stackArray1d< globalIndex, 2*resNumDOF > dofColIndices( 2 * resNumDOF );

    stackArray1d< real64, 2 * numComp > localPerf( 2 * numComp );
    stackArray2d< real64, 2 * resNumDOF * 2 * numComp > localPerfJacobian( 2 * numComp, 2 * resNumDOF );

    // get the reservoir (sub)region and element indices
    localIndex const er  = m_resElementRegion[iperf];
    localIndex const esr = m_resElementSubRegion[iperf];
    localIndex const ei  = m_resElementIndex[iperf];

    // get the well element index for this perforation
    localIndex const iwelem = m_perfWellElemIndex[iperf];
    globalIndex const resOffset = m_resElemDofNumber[er][esr][ei];
    globalIndex const wellElemOffset = m_wellElemDofNumber[iwelem];

    for( integer ic = 0; ic < numComp; ++ic )
    {
      eqnRowIndices[TAG::RES * numComp + ic] = LvArray::integerConversion< localIndex >( resOffset - m_rankOffset ) + ic;
      eqnRowIndices[TAG::WELL * numComp + ic] = LvArray::integerConversion< localIndex >( wellElemOffset - m_rankOffset ) + WJ_ROFFSET::MASSBAL + ic;
    }
    // Note res and well have same col lineup for P and compdens
    for( integer jdof = 0; jdof < NC+1; ++jdof )
    {
      dofColIndices[TAG::RES * resNumDOF + jdof] = resOffset + jdof;
      dofColIndices[TAG::WELL * resNumDOF + jdof] = wellElemOffset + WJ_COFFSET::dP + jdof;
    }
    // For temp its different
    if constexpr ( IS_THERMAL )
    {
      dofColIndices[TAG::RES * resNumDOF + NC+1 ] = resOffset + NC+1;
      dofColIndices[TAG::WELL * resNumDOF + NC+1 ] = wellElemOffset + WJ_COFFSET::dT;
    }
    // populate local flux vector and derivatives
    for( integer ic = 0; ic < numComp; ++ic )
    {
      localPerf[TAG::RES * numComp + ic] = m_dt * m_compPerfRate[iperf][ic];
      localPerf[TAG::WELL * numComp + ic] = -m_dt * m_compPerfRate[iperf][ic];

      if( m_detectCrossflow )
      {
        if( m_compPerfRate[iperf][ic] > LvArray::NumericLimits< real64 >::epsilon )
        {
          m_numCrossFlowPerforations += 1;
        }
      }
      for( integer ke = 0; ke < 2; ++ke )
      {
        localIndex localDofIndexPres = ke * resNumDOF;

        localPerfJacobian[TAG::RES * numComp + ic][localDofIndexPres] = m_dt *  m_dCompPerfRate[iperf][ke][ic][CP_Deriv::dP];
        localPerfJacobian[TAG::WELL * numComp + ic][localDofIndexPres] = -m_dt *  m_dCompPerfRate[iperf][ke][ic][CP_Deriv::dP];
        for( integer jc = 0; jc < numComp; ++jc )
        {
          localIndex const localDofIndexComp = localDofIndexPres + jc + 1;

          localPerfJacobian[TAG::RES * numComp + ic][localDofIndexComp] = m_dt * m_dCompPerfRate[iperf][ke][ic][CP_Deriv::dC+jc];
          localPerfJacobian[TAG::WELL * numComp + ic][localDofIndexComp] = -m_dt * m_dCompPerfRate[iperf][ke][ic][CP_Deriv::dC+jc];
        }
        if constexpr ( IS_THERMAL )
        {
          localIndex localDofIndexTemp  = localDofIndexPres + NC + 1;
          localPerfJacobian[TAG::RES * numComp + ic][localDofIndexTemp] = m_dt *  m_dCompPerfRate[iperf][ke][ic][CP_Deriv::dT];
          localPerfJacobian[TAG::WELL * numComp + ic][localDofIndexTemp] = -m_dt *  m_dCompPerfRate[iperf][ke][ic][CP_Deriv::dT];
        }
      }
    }

    if( m_useTotalMassEquation )
    {
      // Apply equation/variable change transformation(s)
      stackArray1d< real64, 2 * resNumDOF > work( 2 * resNumDOF );
      shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numComp, resNumDOF * 2, 2, localPerfJacobian, work );
      shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numComp, numComp, 2, localPerf );
    }

    for( localIndex i = 0; i < localPerf.size(); ++i )
    {
      if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < m_localMatrix.numRows() )
      {
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                            dofColIndices.data(),
                                                                            localPerfJacobian[i].dataIfContiguous(),
                                                                            2 * resNumDOF );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[eqnRowIndices[i]], localPerf[i] );
      }
    }
    compFluxKernelOp( resOffset, wellElemOffset, dofColIndices, iwelem );

  }


  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElements the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack
   * variables
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

  /// Number of phases
  integer const m_numPhases;

  globalIndex const m_rankOffset;
  // Perfoation variables
  arrayView2d< real64 const > const m_compPerfRate;
  arrayView4d< real64 const > const m_dCompPerfRate;
  arrayView1d< localIndex const > const m_perfWellElemIndex;

  // Element region, subregion, index
  arrayView1d< globalIndex const > const m_wellElemDofNumber;
  ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const m_resElemDofNumber;
  arrayView1d< localIndex const > const m_resElementRegion;
  arrayView1d< localIndex const > const m_resElementSubRegion;
  arrayView1d< localIndex const > const m_resElementIndex;

  // RHS and Jacobian
  arrayView1d< real64 > const m_localRhs;
  CRSMatrixView< real64, globalIndex const >  m_localMatrix;

  bool const m_detectCrossflow;
  integer & m_numCrossFlowPerforations;
  integer const m_useTotalMassEquation;
};

/**
 * @class FaceBasedAssemblyKernelFactory
 */
class IsothermalCompositionalMultiPhaseFluxKernelFactory
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
                   string const wellDofKey,
                   WellElementSubRegion const & subRegion,
                   ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber,
                   PerforationData const * const perforationData,
                   MultiFluidBase const & fluid,
                   integer const & useTotalMassEquation,
                   bool const & detectCrossflow,
                   integer & numCrossFlowPerforations,
                   arrayView1d< real64 > const & localRhs,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix
                   )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();


      BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );


      using kernelType = IsothermalCompositionalMultiPhaseFluxKernel< NUM_COMP, 0 >;


      kernelType kernel( dt, rankOffset, wellDofKey, subRegion, resDofNumber, perforationData, fluid, localRhs, localMatrix, detectCrossflow, numCrossFlowPerforations, kernelFlags );
      kernelType::template launch< POLICY >( perforationData->size(), kernel );
    } );

  }
};


/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NC, integer IS_THERMAL >
class ThermalCompositionalMultiPhaseFluxKernel : public IsothermalCompositionalMultiPhaseFluxKernel< NC, IS_THERMAL >
{
public:
  using Base = IsothermalCompositionalMultiPhaseFluxKernel< NC, IS_THERMAL >;
  /// Compile time value for the number of components
  static constexpr integer numComp = NC;
  static constexpr integer resNumDOF  = NC+1+IS_THERMAL;

  // Well jacobian column and row indicies
  using WJ_COFFSET = compositionalMultiphaseWellKernels::ColOffset_WellJac< NC, IS_THERMAL >;
  using WJ_ROFFSET = compositionalMultiphaseWellKernels::RowOffset_WellJac< NC, IS_THERMAL >;

  using ROFFSET = compositionalMultiphaseWellKernels::RowOffset;
  using COFFSET = compositionalMultiphaseWellKernels::ColOffset;

  using CP_Deriv = multifluid::DerivativeOffsetC< NC, IS_THERMAL >;

  using TAG = compositionalMultiphaseWellKernels::SubRegionTag;

  using Base::m_dt;
  using Base::m_localRhs;
  using Base::m_localMatrix;
  using Base::m_rankOffset;



  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = WJ_COFFSET::nDer;

  /// Compile time value for the number of equations except volume and momentum
  static constexpr integer numEqn = WJ_ROFFSET::nEqn - 2;

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
  ThermalCompositionalMultiPhaseFluxKernel( real64 const dt,
                                            integer const isProducer,
                                            globalIndex const rankOffset,
                                            string const wellDofKey,
                                            WellElementSubRegion const & subRegion,
                                            ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber,
                                            PerforationData const * const perforationData,
                                            MultiFluidBase const & fluid,
                                            arrayView1d< real64 > const & localRhs,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            bool const & detectCrossflow,
                                            integer & numCrossFlowPerforations,
                                            BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags )
    : Base( dt,
            rankOffset,
            wellDofKey,
            subRegion,
            resDofNumber,
            perforationData,
            fluid,
            localRhs,
            localMatrix,
            detectCrossflow,
            numCrossFlowPerforations,
            kernelFlags ),
    m_isProducer( isProducer ),
    m_globalWellElementIndex( subRegion.getGlobalWellElementIndex() ),
    m_energyPerfFlux( perforationData->getField< fields::well::energyPerforationFlux >()),
    m_dEnergyPerfFlux( perforationData->getField< fields::well::dEnergyPerforationFlux >())

  { }


  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] ie the element index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */

  GEOS_HOST_DEVICE
  inline
  void computeFlux( localIndex const iperf ) const
  {
    Base::computeFlux( iperf, [&] ( globalIndex const & resOffset,
                                    globalIndex const & wellElemOffset,
                                    stackArray1d< globalIndex, 2*resNumDOF > & dofColIndices,
                                    localIndex const iwelem )
    {
      // No energy equation if top element and Injector
      // Top element defined by global index == 0
      // Assumption is global index == 0 is top segment with fixed temp BC
      if( !m_isProducer )
      {
        if( m_globalWellElementIndex[iwelem] == 0 )
          return;
      }
      // local working variables and arrays
      stackArray1d< localIndex, 2* numComp > eqnRowIndices( 2 );

      stackArray1d< real64, 2 * numComp > localPerf( 2 );
      stackArray2d< real64, 2 * resNumDOF * 2 * numComp > localPerfJacobian( 2, 2 * resNumDOF );


      // equantion offsets - note res and well have different equation lineups
      eqnRowIndices[TAG::RES  ] = LvArray::integerConversion< localIndex >( resOffset - m_rankOffset )      + NC + 1;
      eqnRowIndices[TAG::WELL ] = LvArray::integerConversion< localIndex >( wellElemOffset - m_rankOffset ) + WJ_ROFFSET::ENERGYBAL;

      // populate local flux vector and derivatives
      localPerf[TAG::RES  ]   = m_dt * m_energyPerfFlux[iperf];
      localPerf[TAG::WELL ]   = -m_dt * m_energyPerfFlux[iperf];

      for( integer ke = 0; ke < 2; ++ke )
      {
        localIndex localDofIndexPres = ke * resNumDOF;
        localPerfJacobian[TAG::RES  ][localDofIndexPres] = m_dt *  m_dEnergyPerfFlux[iperf][ke][CP_Deriv::dP];
        localPerfJacobian[TAG::WELL ][localDofIndexPres] = -m_dt *  m_dEnergyPerfFlux[iperf][ke][CP_Deriv::dP];

        // populate local flux vector and derivatives
        for( integer ic = 0; ic < numComp; ++ic )
        {
          localIndex const localDofIndexComp = localDofIndexPres + ic + 1;
          localPerfJacobian[TAG::RES ][localDofIndexComp] = m_dt * m_dEnergyPerfFlux[iperf][ke][CP_Deriv::dC+ic];
          localPerfJacobian[TAG::WELL][localDofIndexComp] = -m_dt * m_dEnergyPerfFlux[iperf][ke][CP_Deriv::dC+ic];
        }
        localPerfJacobian[TAG::RES ][localDofIndexPres+NC+1] = m_dt * m_dEnergyPerfFlux[iperf][ke][CP_Deriv::dT];
        localPerfJacobian[TAG::WELL][localDofIndexPres+NC+1] = -m_dt * m_dEnergyPerfFlux[iperf][ke][CP_Deriv::dT];
      }


      for( localIndex i = 0; i < localPerf.size(); ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < m_localMatrix.numRows() )
        {
          m_localMatrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                                       dofColIndices.data(),
                                                                                       localPerfJacobian[i].dataIfContiguous(),
                                                                                       2 * resNumDOF );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[eqnRowIndices[i]], localPerf[i] );
        }
      }
    } );


  }


  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElements the number of elements
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack
   * variables
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

  /// Well type
  integer const m_isProducer;

  /// Global index of local element
  arrayView1d< globalIndex const >  m_globalWellElementIndex;

  /// Views on energy flux
  arrayView1d< real64 const > const m_energyPerfFlux;
  arrayView3d< real64 const > const m_dEnergyPerfFlux;
};

/**
 * @class ThermalCompositionalMultiPhaseFluxKernelFactory
 */
class ThermalCompositionalMultiPhaseFluxKernelFactory
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
                   integer const isProducer,
                   real64 const dt,
                   globalIndex const rankOffset,
                   string const wellDofKey,
                   WellElementSubRegion const & subRegion,
                   ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber,
                   PerforationData const * const perforationData,
                   MultiFluidBase const & fluid,
                   integer const & useTotalMassEquation,
                   bool const & detectCrossflow,
                   integer & numCrossFlowPerforations,
                   arrayView1d< real64 > const & localRhs,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix
                   )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();


      BitFlags< isothermalCompositionalMultiphaseBaseKernels::KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseBaseKernels::KernelFlags::TotalMassEquation );


      using kernelType = ThermalCompositionalMultiPhaseFluxKernel< NUM_COMP, 1 >;


      kernelType kernel( dt, isProducer, rankOffset, wellDofKey, subRegion, resDofNumber, perforationData, fluid, localRhs, localMatrix, detectCrossflow, numCrossFlowPerforations, kernelFlags );
      kernelType::template launch< POLICY >( perforationData->size(), kernel );
    } );

  }
};

} // end namespace coupledReservoirAndWellKernels

} // end namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_COUPLEDRESERVOIRANDWELLS_HPP
