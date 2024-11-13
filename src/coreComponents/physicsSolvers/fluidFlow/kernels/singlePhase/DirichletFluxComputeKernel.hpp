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
 * @file DirichletFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_DIRICHLETFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_DIRICHLETFLUXCOMPUTEKERNEL_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/MobilityKernel.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace singlePhaseFVMKernels
{

/******************************** DirichletFluxComputeKernel ********************************/

/**
 * @class DirichletFluxComputeKernel
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_EQN, integer NUM_DOF, typename FLUIDWRAPPER >
class DirichletFluxComputeKernel : public FluxComputeKernel< NUM_EQN, NUM_DOF,
                                                             BoundaryStencilWrapper >
{
public:

  using AbstractBase = singlePhaseFVMKernels::FluxComputeKernelBase;
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

  using Base = singlePhaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF,
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
  DirichletFluxComputeKernel( globalIndex const rankOffset,
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
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
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
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = NoOpFunc{} ) const
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
 * @class DirichletFluxComputeKernelFactory
 */
class DirichletFluxComputeKernelFactory
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
      using kernelType = DirichletFluxComputeKernel< NUM_EQN, NUM_DOF, typename FluidType::KernelWrapper >;

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

} // namespace singlePhaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_DIRICHLETFLUXCOMPUTEKERNEL_HPP
