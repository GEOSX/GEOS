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
 * @file SinglePhaseFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SlurryFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"

namespace geos
{

namespace singlePhaseFVMKernels
{

/******************************** FaceBasedAssemblyKernelBase ********************************/

/**
 * @brief Base class for FaceBasedAssemblyKernel that holds all data not dependent
 *        on template parameters (like stencil type and number of dofs).
 */
class FaceBasedAssemblyKernelBase
{
public:

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using DofNumberAccessor = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::ghostRank,
                      fields::flow::pressure,
                      fields::flow::pressure_n,
                      fields::flow::gravityCoefficient,
                      fields::flow::mobility,
                      fields::flow::dMobility_dPressure >;

  using SinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity_dPressure >;

  using SlurryFluidAccessors =
    StencilMaterialAccessors< constitutive::SlurryFluidBase,
                              fields::singlefluid::density,
                              fields::singlefluid::dDensity_dPressure >;

  using PermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure >;

  using ProppantPermeabilityAccessors =
    StencilMaterialAccessors< constitutive::PermeabilityBase,
                              fields::permeability::permeability,
                              fields::permeability::dPerm_dPressure,
                              fields::permeability::dPerm_dDispJump,
                              fields::permeability::permeabilityMultiplier >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofNumberAccessor accessor for the dof numbers
   * @param[in] singleFlowAccessors accessor for wrappers registered by the solver
   * @param[in] singlePhaseFluidAccessors accessor for wrappers registered by the singlefluid model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernelBase( globalIndex const rankOffset,
                               DofNumberAccessor const & dofNumberAccessor,
                               SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                               SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                               PermeabilityAccessors const & permeabilityAccessors,
                               real64 const & dt,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs )
    : m_rankOffset( rankOffset ),
    m_dt( dt ),
    m_dofNumber( dofNumberAccessor.toNestedViewConst() ),
    m_permeability( permeabilityAccessors.get( fields::permeability::permeability {} ) ),
    m_dPerm_dPres( permeabilityAccessors.get( fields::permeability::dPerm_dPressure {} ) ),
    m_ghostRank( singlePhaseFlowAccessors.get( fields::ghostRank {} ) ),
    m_gravCoef( singlePhaseFlowAccessors.get( fields::flow::gravityCoefficient {} ) ),
    m_pres( singlePhaseFlowAccessors.get( fields::flow::pressure {} ) ),
    m_mob( singlePhaseFlowAccessors.get( fields::flow::mobility {} ) ),
    m_dMob_dPres( singlePhaseFlowAccessors.get( fields::flow::dMobility_dPressure {} ) ),
    m_dens( singlePhaseFluidAccessors.get( fields::singlefluid::density {} ) ),
    m_dDens_dPres( singlePhaseFluidAccessors.get( fields::singlefluid::dDensity_dPressure {} ) ),
    m_localMatrix( localMatrix ),
    m_localRhs( localRhs )
  {}

protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// Time step size
  real64 const m_dt;

  /// Views on dof numbers
  ElementViewConst< arrayView1d< globalIndex const > > const m_dofNumber;

  /// Views on permeability
  ElementViewConst< arrayView3d< real64 const > > m_permeability;
  ElementViewConst< arrayView3d< real64 const > > m_dPerm_dPres;

  /// Views on ghost rank numbers and gravity coefficients
  ElementViewConst< arrayView1d< integer const > > const m_ghostRank;
  ElementViewConst< arrayView1d< real64 const > > const m_gravCoef;

  // Primary and secondary variables
  /// Views on pressure
  ElementViewConst< arrayView1d< real64 const > > const m_pres;

  /// Views on fluid mobility
  ElementViewConst< arrayView1d< real64 const > > const m_mob;
  ElementViewConst< arrayView1d< real64 const > > const m_dMob_dPres;

  /// Views on fluid density
  ElementViewConst< arrayView2d< real64 const > > const m_dens;
  ElementViewConst< arrayView2d< real64 const > > const m_dDens_dPres;

  // Residual and jacobian

  /// View on the local CRS matrix
  CRSMatrixView< real64, globalIndex const > const m_localMatrix;
  /// View on the local RHS
  arrayView1d< real64 > const m_localRhs;
};

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public FaceBasedAssemblyKernelBase
{
public:

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations
  static constexpr integer numEqn = NUM_EQN;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] singlePhaseFlowAccessors
   * @param[in] singlePhaseFluidAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  FaceBasedAssemblyKernel( globalIndex const rankOffset,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                           SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
    : FaceBasedAssemblyKernelBase( rankOffset,
                                   dofNumberAccessor,
                                   singlePhaseFlowAccessors,
                                   singlePhaseFluidAccessors,
                                   permeabilityAccessors,
                                   dt,
                                   localMatrix,
                                   localRhs ),
    m_stencilWrapper( stencilWrapper ),
    m_seri( stencilWrapper.getElementRegionIndices() ),
    m_sesri( stencilWrapper.getElementSubRegionIndices() ),
    m_sei( stencilWrapper.getElementIndices() )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : stencilSize( size ),
      numFluxElems( numElems ),
      dofColIndices( size * numDof ),
      localFlux( numElems * numEqn ),
      localFluxJacobian( numElems * numEqn, size * numDof )
    {}

    // Stencil information

    /// Stencil size for a given connection
    localIndex const stencilSize;

    /// Number of elements for a given connection
    localIndex const numFluxElems;

    // Transmissibility and derivatives

    /// Transmissibility
    real64 transmissibility[maxNumConns][2]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dPres[maxNumConns][2]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    stackArray1d< globalIndex, maxNumElems * numDof > dofColIndices;

    /// Storage for the face local residual vector (all equations except volume balance)
    stackArray1d< real64, maxNumElems * numEqn > localFlux;
    /// Storage for the face local Jacobian matrix
    stackArray2d< real64, maxNumElems * numEqn * maxStencilSize * numDof > localFluxJacobian;

  };

  /**
   * @brief Getter for the stencil size at this connection
   * @param[in] iconn the connection index
   * @return the size of the stencil at this connection
   */
  GEOS_HOST_DEVICE
  localIndex stencilSize( localIndex const iconn ) const
  { return m_sei[iconn].size(); }

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOS_HOST_DEVICE
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }

  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    // set degrees of freedom indices for this face
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      globalIndex const offset = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];

      for( integer jdof = 0; jdof < numDof; ++jdof )
      {
        stack.dofColIndices[i * numDof + jdof] = offset + jdof;
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const
  {
    // first, compute the transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     stack.transmissibility,
                                     stack.dTrans_dPres );

    localIndex k[2];
    localIndex connectionIndex = 0;

    for( k[0] = 0; k[0] < stack.numFluxElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numFluxElems; ++k[1] )
      {
        real64 fluxVal = 0.0;
        real64 dFlux_dTrans = 0.0;
        real64 alpha = 0.0;
        real64 mobility = 0.0;
        real64 potGrad = 0.0;
        real64 trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };
        real64 dTrans[2] = { stack.dTrans_dPres[connectionIndex][0], stack.dTrans_dPres[connectionIndex][1] };
        real64 dFlux_dP[2] = {0.0, 0.0};
        localIndex const regionIndex[2]    = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const subRegionIndex[2] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const elementIndex[2]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        fluxKernelsHelper::computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                                                   trans,
                                                   dTrans,
                                                   m_pres,
                                                   m_gravCoef,
                                                   m_dens,
                                                   m_dDens_dPres,
                                                   m_mob,
                                                   m_dMob_dPres,
                                                   alpha,
                                                   mobility,
                                                   potGrad,
                                                   fluxVal,
                                                   dFlux_dP,
                                                   dFlux_dTrans );

        // populate local flux vector and derivatives
        stack.localFlux[k[0]*numEqn] += m_dt * fluxVal;
        stack.localFlux[k[1]*numEqn] -= m_dt * fluxVal;

        for( integer ke = 0; ke < 2; ++ke )
        {
          localIndex const localDofIndexPres = k[ke] * numDof;
          stack.localFluxJacobian[k[0]*numEqn][localDofIndexPres] += m_dt * dFlux_dP[ke];
          stack.localFluxJacobian[k[1]*numEqn][localDofIndexPres] -= m_dt * dFlux_dP[ke];
        }

        // Customize the kernel with this lambda
        kernelOp( k, regionIndex, subRegionIndex, elementIndex, connectionIndex, alpha, mobility, potGrad, fluxVal, dFlux_dP );

        connectionIndex++;
      }
    }
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && kernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const
  {
    // add contribution to residual and jacobian into:
    // - the mass balance equation
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    for( integer i = 0; i < stack.numFluxElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GT( m_localMatrix.numRows(), localRow );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow], stack.localFlux[i * numEqn] );
        m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            stack.dofColIndices.data(),
                                                                            stack.localFluxJacobian[i * numEqn].dataIfContiguous(),
                                                                            stack.stencilSize * numDof );

        // call the lambda to assemble additional terms, such as thermal terms
        kernelOp( i, localRow );
      }
    }
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numConnections the number of connections
   * @param[inout] kernelComponent the kernel component providing access to setup/compute/complete functions and stack variables
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numConnections,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    forAll< POLICY >( numConnections, [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
                                                  kernelComponent.numPointsInFlux( iconn ) );

      kernelComponent.setup( iconn, stack );
      kernelComponent.computeFlux( iconn, stack );
      kernelComponent.complete( iconn, stack );
    } );
  }


protected:

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;
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
   * @tparam STENCILWRAPPER the type of the stencil wrapper
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_EQN = 1;
    integer constexpr NUM_DOF = 1;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, STENCILWRAPPER >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );

    kernelType kernel( rankOffset, stencilWrapper, dofNumberAccessor,
                       flowAccessors, fluidAccessors, permAccessors,
                       dt, localMatrix, localRhs );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

/******************************** DirichletFaceBasedAssemblyKernel ********************************/

/**
 * @class DirichFaceBasedAssemblyKernel
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename FLUIDWRAPPER >
class DirichletFaceBasedAssemblyKernel : public FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF,
                                                                         BoundaryStencilWrapper >
{
public:

  using AbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;
  using AbstractBase::m_mob;
  using AbstractBase::m_dMob_dPres;
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens_dPres;
  using AbstractBase::m_permeability;
  using AbstractBase::m_dPerm_dPres;

  using AbstractBase::m_localMatrix;
  using AbstractBase::m_localRhs;

  using Base = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF,
                                                               BoundaryStencilWrapper >;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
   * @param[in] dofNumberAccessor
   * @param[in] singlePhaseFlowAccessors
   * @param[in] singlePhaseFluidAccessors
   * @param[in] permeabilityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  DirichletFaceBasedAssemblyKernel( globalIndex const rankOffset,
                                    FaceManager const & faceManager,
                                    BoundaryStencilWrapper const & stencilWrapper,
                                    FLUIDWRAPPER const & fluidWrapper,
                                    DofNumberAccessor const & dofNumberAccessor,
                                    SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                    SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                    PermeabilityAccessors const & permeabilityAccessors,
                                    real64 const & dt,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
    : Base( rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_facePres( faceManager.getField< fields::flow::facePressure >() ),
    m_faceGravCoef( faceManager.getField< fields::flow::gravityCoefficient >() ),
    m_fluidWrapper( fluidWrapper )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const GEOS_UNUSED_PARAM( size ),
                    localIndex GEOS_UNUSED_PARAM( numElems ) )
    {}

    /// Transmissibility
    real64 transmissibility = 0.0;

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    globalIndex dofColIndices[numDof]{};

    /// Storage for the face local residual
    real64 localFlux[numEqn]{};

    /// Storage for the face local Jacobian
    real64 localFluxJacobian[numEqn][numDof]{};

  };


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const iconn,
              StackVariables & stack ) const
  {
    globalIndex const offset =
      m_dofNumber[m_seri( iconn, BoundaryStencil::Order::ELEM )][m_sesri( iconn, BoundaryStencil::Order::ELEM )][m_sei( iconn, BoundaryStencil::Order::ELEM )];

    for( integer jdof = 0; jdof < numDof; ++jdof )
    {
      stack.dofColIndices[jdof] = offset + jdof;
    }
  }

  /**
   * @brief Compute the local Dirichlet face flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] compFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const
  {
    using Order = BoundaryStencil::Order;
    localIndex constexpr numElems = BoundaryStencil::maxNumPointsInFlux;

    stackArray1d< real64, numElems > mobility( numElems );
    stackArray1d< real64, numElems > dMobility_dP( numElems );

    localIndex const er  = m_seri( iconn, Order::ELEM );
    localIndex const esr = m_sesri( iconn, Order::ELEM );
    localIndex const ei  = m_sei( iconn, Order::ELEM );
    localIndex const kf  = m_sei( iconn, Order::FACE );

    // Get flow quantities on the elem/face
    real64 faceDens, faceVisc;
    constitutive::SingleFluidBaseUpdate::computeValues( m_fluidWrapper, m_facePres[kf], faceDens, faceVisc );

    mobility[Order::ELEM] = m_mob[er][esr][ei];
    singlePhaseBaseKernels::MobilityKernel::compute( faceDens, faceVisc, mobility[Order::FACE] );

    dMobility_dP[Order::ELEM] = m_dMob_dPres[er][esr][ei];
    dMobility_dP[Order::FACE] = 0.0;

    // Compute average density
    real64 const densMean = 0.5 * ( m_dens[er][esr][ei][0] + faceDens );
    real64 const dDens_dP = 0.5 * m_dDens_dPres[er][esr][ei][0];

    // Evaluate potential difference
    real64 const potDif = (m_pres[er][esr][ei] - m_facePres[kf])
                          - densMean * (m_gravCoef[er][esr][ei] - m_faceGravCoef[kf]);
    real64 const dPotDif_dP = 1.0 - dDens_dP * m_gravCoef[er][esr][ei];

    real64 dTrans_dPerm[3];
    m_stencilWrapper.computeWeights( iconn, m_permeability, stack.transmissibility, dTrans_dPerm );
    real64 const dTrans_dPres = LvArray::tensorOps::AiBi< 3 >( dTrans_dPerm, m_dPerm_dPres[er][esr][ei][0] );

    real64 const f = stack.transmissibility * potDif;
    real64 const dF_dP = stack.transmissibility * dPotDif_dP + dTrans_dPres * potDif;

    // Upwind mobility
    localIndex const k_up = ( potDif >= 0 ) ? Order::ELEM : Order::FACE;
    real64 const mobility_up = mobility[k_up];
    real64 const dMobility_dP_up = dMobility_dP[k_up];

    // call the lambda in the phase loop to allow the reuse of the fluxes and their derivatives
    // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase

    compFluxKernelOp( er, esr, ei, kf, f, dF_dP, mobility_up, dMobility_dP_up );

    // *** end of upwinding

    // Populate local flux vector and derivatives

    stack.localFlux[0] = m_dt * mobility[k_up] * f;
    stack.localFluxJacobian[0][0] = m_dt * ( mobility_up * dF_dP + dMobility_dP_up * f );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  template< typename FUNC = singlePhaseBaseKernels::NoOpFunc >
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = singlePhaseBaseKernels::NoOpFunc{} ) const
  {
    using Order = BoundaryStencil::Order;

    localIndex const er = m_seri( iconn, Order::ELEM );
    localIndex const esr = m_sesri( iconn, Order::ELEM );
    localIndex const ei = m_sei( iconn, Order::ELEM );

    if( m_ghostRank[er][esr][ei] < 0 )
    {
      // Add to global residual/jacobian
      globalIndex const dofIndex = m_dofNumber[er][esr][ei];
      localIndex const localRow = LvArray::integerConversion< localIndex >( dofIndex - m_rankOffset );

      RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow], stack.localFlux[0] );

      AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
        ( localRow,
        stack.dofColIndices,
        stack.localFluxJacobian[0],
        numDof );

      assemblyKernelOp( localRow );
    }
  }

protected:

  /// Views on face pressure, temperature, and composition
  arrayView1d< real64 const > const m_facePres;

  /// View on the face gravity coefficient
  arrayView1d< real64 const > const m_faceGravCoef;

  /// Reference to the fluid wrapper
  FLUIDWRAPPER const m_fluidWrapper;

};


/**
 * @class DirichletFaceBasedAssemblyKernelFactory
 */
class DirichletFaceBasedAssemblyKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] faceManager reference to the face manager
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the boundary stencil wrapper
   * @param[in] fluidBase the single phase fluid constitutive model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & dofKey,
                   string const & solverName,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   BoundaryStencilWrapper const & stencilWrapper,
                   constitutive::SingleFluidBase & fluidBase,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
    {
      using FluidType = TYPEOFREF( fluid );
      typename FluidType::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

      integer constexpr NUM_DOF = 1;
      integer constexpr NUM_EQN = 1;
      using kernelType = DirichletFaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, typename FluidType::KernelWrapper >;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );

      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      typename kernelType::SinglePhaseFlowAccessors singlePhaseFlowAccessors( elemManager, solverName );
      typename kernelType::SinglePhaseFluidAccessors singlePhaseFluidAccessors( elemManager, solverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

      kernelType kernel( rankOffset,
                         faceManager,
                         stencilWrapper,
                         fluidWrapper,
                         dofNumberAccessor,
                         singlePhaseFlowAccessors,
                         singlePhaseFluidAccessors,
                         permeabilityAccessors,
                         dt,
                         localMatrix,
                         localRhs );

      kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }

};


/******************************** AquiferBCKernel ********************************/

/**
 * @brief Functions to assemble aquifer boundary condition contributions to residual and Jacobian
 */
struct AquiferBCKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  GEOS_HOST_DEVICE
  static void
  compute( real64 const & aquiferVolFlux,
           real64 const & dAquiferVolFlux_dPres,
           real64 const & aquiferDens,
           real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & dt,
           real64 & localFlux,
           real64 & localFluxJacobian )
  {
    if( aquiferVolFlux > 0 ) // aquifer is upstream
    {
      localFlux -= dt * aquiferVolFlux * aquiferDens;
      localFluxJacobian -= dt * dAquiferVolFlux_dPres * aquiferDens;
    }
    else // reservoir is upstream
    {
      localFlux -= dt * aquiferVolFlux * dens;
      localFluxJacobian -= dt * (dAquiferVolFlux_dPres * dens + aquiferVolFlux * dDens_dPres);
    }
  }

  static void
  launch( BoundaryStencil const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          real64 const & aquiferDens,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & pres_n,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          real64 const & timeAtBeginningOfStep,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    using Order = BoundaryStencil::Order;

    BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
    BoundaryStencil::WeightContainerViewConstType const & weight = stencil.getWeights();

    forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      // working variables
      real64 localFlux = 0.0;
      real64 localFluxJacobian = 0.0;

      localIndex const er  = seri( iconn, Order::ELEM );
      localIndex const esr = sesri( iconn, Order::ELEM );
      localIndex const ei  = sefi( iconn, Order::ELEM );
      real64 const areaFraction = weight( iconn, Order::ELEM );

      // compute the aquifer influx rate using the pressure influence function and the aquifer props
      real64 dAquiferVolFlux_dPres = 0.0;
      real64 const aquiferVolFlux = aquiferBCWrapper.compute( timeAtBeginningOfStep,
                                                              dt,
                                                              pres[er][esr][ei],
                                                              pres_n[er][esr][ei],
                                                              gravCoef[er][esr][ei],
                                                              areaFraction,
                                                              dAquiferVolFlux_dPres );

      // compute the phase/component aquifer flux
      AquiferBCKernel::compute( aquiferVolFlux,
                                dAquiferVolFlux_dPres,
                                aquiferDens,
                                dens[er][esr][ei][0],
                                dDens_dPres[er][esr][ei][0],
                                dt,
                                localFlux,
                                localFluxJacobian );

      // Add to residual/jacobian
      if( ghostRank[er][esr][ei] < 0 )
      {
        globalIndex const globalRow = dofNumber[er][esr][ei];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GT( localMatrix.numRows(), localRow );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux );
        localMatrix.addToRow< parallelDeviceAtomic >( localRow,
                                                      &dofNumber[er][esr][ei],
                                                      &localFluxJacobian,
                                                      1 );
      }
    } );
  }

};


} // namespace singlePhaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
