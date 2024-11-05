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
 * @file SinglePhasePoromechanicsEmbeddedFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP

#include "physicsSolvers/fluidFlow/kernels/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/FluxKernelsHelper.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace singlePhasePoromechanicsEmbeddedFracturesKernels
{

template< integer NUM_EQN, integer NUM_DOF >
class ConnectorBasedAssemblyKernel : public singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper >
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

  using AbstractBase = singlePhaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using SinglePhaseFlowAccessors = AbstractBase::SinglePhaseFlowAccessors;
  using SinglePhaseFluidAccessors = AbstractBase::SinglePhaseFluidAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;
  using FracturePermeabilityAccessors = StencilMaterialAccessors< constitutive::PermeabilityBase,
                                                                  fields::permeability::dPerm_dDispJump >;
  using AbstractBase::m_dt;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_permeability;
  using AbstractBase::m_dPerm_dPres;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;
  using AbstractBase::m_mob;
  using AbstractBase::m_dMob_dPres;
  using AbstractBase::m_dens;
  using AbstractBase::m_dDens_dPres;

  using Base = singlePhaseFVMKernels::FaceBasedAssemblyKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;
  using Base::m_ghostRank;



  ConnectorBasedAssemblyKernel( globalIndex const rankOffset,
                                SurfaceElementStencilWrapper const & stencilWrapper,
                                DofNumberAccessor const & flowDofNumberAccessor,
                                DofNumberAccessor const & dispJumpDofNumberAccessor,
                                SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                PermeabilityAccessors const & permeabilityAccessors,
                                FracturePermeabilityAccessors const & edfmPermeabilityAccessors,
                                real64 const & dt,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs )
    : Base( rankOffset,
            stencilWrapper,
            flowDofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_dispJumpDofNumber( dispJumpDofNumberAccessor.toNestedViewConst() ),
    m_dPerm_dDispJump( edfmPermeabilityAccessors.get( fields::permeability::dPerm_dDispJump {} ) )
  {}


  /**
   * @struct StackVariables
   * @brief Kernel variables (dof numbers, jacobian and residual) located on the stack
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /**
     * @brief Constructor for the stack variables
     * @param[in] size size of the stencil for this connection
     * @param[in] numElems number of elements for this connection
     */
    GEOS_HOST_DEVICE
    StackVariables( localIndex const size, localIndex numElems )
      : Base::StackVariables( size, numElems )
    {}

    /// Derivatives of transmissibility with respect to the dispJump
    real64 dTrans_dDispJump[maxNumConns][2][3]{};
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
    // set degrees of freedom indices for this face
    for( integer i = 0; i < stack.stencilSize; ++i )
    {
      localIndex localDofIndex = numDof * i;
      stack.dofColIndices[ localDofIndex ]     = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
      stack.dofColIndices[ localDofIndex + 1 ] = m_dispJumpDofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
      stack.dofColIndices[ localDofIndex + 2 ] = m_dispJumpDofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] + 1;
      stack.dofColIndices[ localDofIndex + 3 ] = m_dispJumpDofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] + 2;
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the flux
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] NoOpFunc the function used to customize the computation of the flux
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    StackVariables & stack,
                    FUNC && kernelOp = NoOpFunc{} ) const

  {

    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     m_dPerm_dDispJump,
                                     stack.transmissibility,
                                     stack.dTrans_dPres,
                                     stack.dTrans_dDispJump );


    real64 fluxVal = 0.0;
    real64 dFlux_dTrans = 0.0;
    /// EDFM connections are always only between 2 elements. There are no star connections.
    real64 trans[2] = {stack.transmissibility[0][0], stack.transmissibility[0][1]};
    real64 dTrans[2] = { stack.dTrans_dPres[0][0], stack.dTrans_dPres[0][1] };
    real64 dFlux_dP[2] = {0.0, 0.0};
    localIndex const regionIndex[2]    = {m_seri[iconn][0], m_seri[iconn][1]};
    localIndex const subRegionIndex[2] = {m_sesri[iconn][0], m_sesri[iconn][1]};
    localIndex const elementIndex[2]   = {m_sei[iconn][0], m_sei[iconn][1]};
    real64 alpha = 0.0;
    real64 mobility = 0.0;
    real64 potGrad = 0.0;

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
    stack.localFlux[0] =  m_dt * fluxVal;
    stack.localFlux[1] = -m_dt * fluxVal;

    real64 dFlux_dDispJump[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    for( localIndex i=0; i < 3; i++ )
    {
      dFlux_dDispJump[0][i] =   dFlux_dTrans * stack.dTrans_dDispJump[0][0][i];
      dFlux_dDispJump[1][i] = -dFlux_dTrans * stack.dTrans_dDispJump[0][1][i];
    }
    for( localIndex ke = 0; ke < 2; ++ke )
    {
      localIndex const dofIndex = numDof*ke;

      stack.localFluxJacobian[0][dofIndex]   =  m_dt * dFlux_dP[ke];
      stack.localFluxJacobian[0][dofIndex+1] =  m_dt * dFlux_dDispJump[ke][0];
      stack.localFluxJacobian[0][dofIndex+2] =  m_dt * dFlux_dDispJump[ke][1];
      stack.localFluxJacobian[0][dofIndex+3] =  m_dt * dFlux_dDispJump[ke][2];

      stack.localFluxJacobian[1][dofIndex]   = -m_dt * dFlux_dP[ke];
      stack.localFluxJacobian[1][dofIndex+1] = -m_dt * dFlux_dDispJump[ke][0];
      stack.localFluxJacobian[1][dofIndex+2] = -m_dt * dFlux_dDispJump[ke][1];
      stack.localFluxJacobian[1][dofIndex+3] = -m_dt * dFlux_dDispJump[ke][2];
    }

    kernelOp( regionIndex, subRegionIndex, elementIndex, iconn, alpha, mobility, potGrad, fluxVal, dFlux_dTrans, dFlux_dP );
  }


private:

  ElementViewConst< arrayView1d< globalIndex const > > const m_dispJumpDofNumber;

  ElementViewConst< arrayView4d< real64 const > > const m_dPerm_dDispJump;

};



/**
 * @class FaceBasedAssemblyKernelFactory
 */
class ConnectorBasedAssemblyKernelFactory
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
  template< typename POLICY >
  static void
  createAndLaunch( globalIndex const rankOffset,
                   string const & pressureDofKey,
                   string const & dispJumpDofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   SurfaceElementStencilWrapper const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    integer constexpr NUM_DOF = 4; // pressure + jumps
    integer constexpr NUM_EQN = 1; // pressure + jumps


    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > pressureDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( pressureDofKey );
    pressureDofNumberAccessor.setName( solverName + "/accessors/" + pressureDofKey );

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dispJumpDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dispJumpDofKey );
    dispJumpDofNumberAccessor.setName( solverName + "/accessors/" + dispJumpDofKey );

    using kernelType = ConnectorBasedAssemblyKernel< NUM_EQN, NUM_DOF >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
    typename kernelType::FracturePermeabilityAccessors edfmPermAccessors( elemManager, solverName );


    kernelType kernel( rankOffset, stencilWrapper,
                       pressureDofNumberAccessor, dispJumpDofNumberAccessor,
                       flowAccessors, fluidAccessors, permAccessors, edfmPermAccessors,
                       dt, localMatrix, localRhs );

    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};



} // namespace SinglePhaseProppantFluxKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEMBEDDEDFRACTURES_HPP
