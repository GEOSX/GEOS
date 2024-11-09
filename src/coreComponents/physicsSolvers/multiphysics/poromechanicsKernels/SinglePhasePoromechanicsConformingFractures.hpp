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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP

#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxComputeKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/FluxKernelsHelper.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace singlePhasePoromechanicsConformingFracturesKernels
{

template< integer NUM_EQN, integer NUM_DOF >
class ConnectorBasedAssemblyKernel : public singlePhaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper >
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

  using AbstractBase = singlePhaseFVMKernels::FluxComputeKernelBase;
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

  using Base = singlePhaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  ConnectorBasedAssemblyKernel( globalIndex const rankOffset,
                                SurfaceElementStencilWrapper const & stencilWrapper,
                                DofNumberAccessor const & flowDofNumberAccessor,
                                SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                                SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                                PermeabilityAccessors const & permeabilityAccessors,
                                FracturePermeabilityAccessors const & fracturePermeabilityAccessors,
                                real64 const & dt,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs,
                                CRSMatrixView< real64, localIndex const > const & dR_dAper )
    : Base( rankOffset,
            stencilWrapper,
            flowDofNumberAccessor,
            singlePhaseFlowAccessors,
            singlePhaseFluidAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs ),
    m_dR_dAper( dR_dAper ),
    m_dPerm_dDispJump( fracturePermeabilityAccessors.get( fields::permeability::dPerm_dDispJump {} ) )
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
      : Base::StackVariables( size, numElems ),
      localColIndices( numElems ),
      dFlux_dAperture( numElems, size )
    {}

    stackArray1d< localIndex, maxNumElems > localColIndices;

    stackArray2d< real64, maxNumElems * maxStencilSize > dFlux_dAperture;

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
      globalIndex const offset = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
      for( integer jdof = 0; jdof < numDof; ++jdof )
      {
        stack.dofColIndices[i * numDof + jdof] = offset + jdof;
      }
      stack.localColIndices[ i ] = m_sei( iconn, i );
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


    localIndex k[2];
    localIndex connectionIndex = 0;
    for( k[0]=0; k[0]<stack.numFluxElems; ++k[0] )
    {
      for( k[1]=k[0]+1; k[1]<stack.numFluxElems; ++k[1] )
      {
        real64 fluxVal = 0.0;
        real64 dFlux_dTrans = 0.0;
        real64 alpha = 0.0;
        real64 mobility = 0.0;
        real64 potGrad = 0.0;
        real64 const trans[2] = { stack.transmissibility[connectionIndex][0], stack.transmissibility[connectionIndex][1] };
        real64 const dTrans[2] = { stack.dTrans_dPres[connectionIndex][0], stack.dTrans_dPres[connectionIndex][1] };
        real64 dFlux_dP[2] = { 0.0, 0.0 };
        localIndex const regionIndex[2]    = {m_seri[iconn][k[0]], m_seri[iconn][k[1]]};
        localIndex const subRegionIndex[2] = {m_sesri[iconn][k[0]], m_sesri[iconn][k[1]]};
        localIndex const elementIndex[2]   = {m_sei[iconn][k[0]], m_sei[iconn][k[1]]};

        singlePhaseFluxKernelsHelper::computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
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
        stack.localFlux[k[0]* numDof] += m_dt * fluxVal;
        stack.localFlux[k[1]* numDof] -= m_dt * fluxVal;

        real64 dFlux_dAper[2] = {0.0, 0.0};
        dFlux_dAper[0] =  m_dt * dFlux_dTrans * stack.dTrans_dDispJump[connectionIndex][0][0];
        dFlux_dAper[1] = -m_dt * dFlux_dTrans * stack.dTrans_dDispJump[connectionIndex][1][0];

        stack.localFluxJacobian[k[0]*numEqn][k[0]* numDof] += dFlux_dP[0] * m_dt;
        stack.localFluxJacobian[k[0]*numEqn][k[1]* numDof] += dFlux_dP[1] * m_dt;
        stack.localFluxJacobian[k[1]*numEqn][k[0]* numDof] -= dFlux_dP[0] * m_dt;
        stack.localFluxJacobian[k[1]*numEqn][k[1]* numDof] -= dFlux_dP[1] * m_dt;

        stack.dFlux_dAperture[k[0]][k[0]] += dFlux_dAper[0];
        stack.dFlux_dAperture[k[0]][k[1]] += dFlux_dAper[1];
        stack.dFlux_dAperture[k[1]][k[0]] -= dFlux_dAper[0];
        stack.dFlux_dAperture[k[1]][k[1]] -= dFlux_dAper[1];

        kernelOp( k, regionIndex, subRegionIndex, elementIndex, iconn, alpha, mobility, potGrad, fluxVal, dFlux_dTrans, dFlux_dP );
        connectionIndex++;
      }
    }
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
                 FUNC && kernelOp = NoOpFunc{} ) const
  {
    // Call Base::complete to assemble the mass balance equations
    // In the lambda, fill the dR_dAper matrix
    Base::complete( iconn, stack, [&] ( integer const i,
                                        localIndex const localRow )
    {

      localIndex const row = LvArray::integerConversion< localIndex >( m_sei( iconn, i ) );

      m_dR_dAper.addToRowBinarySearch< parallelDeviceAtomic >( row,
                                                               stack.localColIndices.data(),
                                                               stack.dFlux_dAperture[i].dataIfContiguous(),
                                                               stack.stencilSize );
      // call the lambda to assemble additional terms, such as thermal terms
      kernelOp( i, localRow );
    } );
  }

private:

  CRSMatrixView< real64, localIndex const > m_dR_dAper;

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
                   string const & dofKey,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   SurfaceElementStencilWrapper const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   CRSMatrixView< real64, localIndex const > const & dR_dAper )
  {
    integer constexpr NUM_DOF = 1; // pressure
    integer constexpr NUM_EQN = 1;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > flowDofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    flowDofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = ConnectorBasedAssemblyKernel< NUM_EQN, NUM_DOF >;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permAccessors( elemManager, solverName );
    typename kernelType::FracturePermeabilityAccessors fracPermAccessors( elemManager, solverName );

    kernelType kernel( rankOffset, stencilWrapper, flowDofNumberAccessor,
                       flowAccessors, fluidAccessors, permAccessors, fracPermAccessors,
                       dt, localMatrix, localRhs, dR_dAper );

    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  }
};

} // namespace SinglePhasePoromechanicsConformingFracturesKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
