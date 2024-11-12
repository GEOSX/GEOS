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
 * @file MultiphasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP

#include "physicsSolvers/fluidFlow/kernels/compositional/FluxComputeKernel.hpp"
//#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"

namespace geos
{

namespace multiphasePoromechanicsConformingFracturesKernels
{

template< integer NUM_EQN, integer NUM_DOF >
class FluxComputeKernel : public isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper >
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

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;
  using FracturePermeabilityAccessors = StencilMaterialAccessors< constitutive::PermeabilityBase,
                                                                  fields::permeability::dPerm_dDispJump >;

  using AbstractBase::m_dt;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernel< NUM_EQN, NUM_DOF, SurfaceElementStencilWrapper >;
  using Base::numDof;
  using Base::numEqn;
  using Base::maxNumElems;
  using Base::maxNumConns;
  using Base::maxStencilSize;
  using Base::m_stencilWrapper;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;
  using Base::numFluxSupportPoints;
  using Base::numComp;
  using Base::m_numPhases;
  using Base::m_kernelFlags;

  using Base::m_permeability;
  using Base::m_dPerm_dPres;
  using Base::m_phaseMob;
  using Base::m_dPhaseMob;
  using Base::m_phaseVolFrac;
  using Base::m_dPhaseVolFrac;
  using Base::m_phaseCompFrac;
  using Base::m_dPhaseCompFrac;
  using Base::m_dCompFrac_dCompDens;
  using Base::m_phaseMassDens;
  using Base::m_dPhaseMassDens;
  using Base::m_phaseCapPressure;
  using Base::m_dPhaseCapPressure_dPhaseVolFrac;

  FluxComputeKernel( integer const numPhases,
                     globalIndex const rankOffset,
                     SurfaceElementStencilWrapper const & stencilWrapper,
                     DofNumberAccessor const & dofNumberAccessor,
                     CompFlowAccessors const & compFlowAccessors,
                     MultiFluidAccessors const & multiFluidAccessors,
                     CapPressureAccessors const & capPressureAccessors,
                     PermeabilityAccessors const & permeabilityAccessors,
                     FracturePermeabilityAccessors const & fracturePermeabilityAccessors,
                     real64 const dt,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs,
                     BitFlags< isothermalCompositionalMultiphaseFVMKernels::KernelFlags > kernelFlags,
                     CRSMatrixView< real64, localIndex const > const & dR_dAper )
    : Base( numPhases,
            rankOffset,
            stencilWrapper,
            dofNumberAccessor,
            compFlowAccessors,
            multiFluidAccessors,
            capPressureAccessors,
            permeabilityAccessors,
            dt,
            localMatrix,
            localRhs,
            kernelFlags ),
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
      // TODO
      dFlux_dAperture( numComp, numElems, size )
    {}

    stackArray1d< localIndex, maxNumElems > localColIndices;

    // TODO
    stackArray3d< real64, numEqn * maxNumElems * maxStencilSize > dFlux_dAperture;

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
                    FUNC && compFluxKernelOp = NoOpFunc{} ) const
  {
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     m_dPerm_dPres,
                                     m_dPerm_dDispJump,
                                     stack.transmissibility,
                                     stack.dTrans_dPres,
                                     stack.dTrans_dDispJump );

    localIndex k[numFluxSupportPoints];
    localIndex connectionIndex = 0;
    for( k[0] = 0; k[0] < stack.numConnectedElems; ++k[0] )
    {
      for( k[1] = k[0] + 1; k[1] < stack.numConnectedElems; ++k[1] )
      {
        /// cell indices
        localIndex const seri[numFluxSupportPoints]  = {m_seri( iconn, k[0] ), m_seri( iconn, k[1] )};
        localIndex const sesri[numFluxSupportPoints] = {m_sesri( iconn, k[0] ), m_sesri( iconn, k[1] )};
        localIndex const sei[numFluxSupportPoints]   = {m_sei( iconn, k[0] ), m_sei( iconn, k[1] )};

        // clear working arrays
        real64 compFlux[numComp]{};
        real64 dCompFlux_dP[numFluxSupportPoints][numComp]{};
        real64 dCompFlux_dC[numFluxSupportPoints][numComp][numComp]{};
        real64 dCompFlux_dTrans[numComp]{};

        real64 const trans[numFluxSupportPoints] = { stack.transmissibility[connectionIndex][0],
                                                     stack.transmissibility[connectionIndex][1] };

        real64 const dTrans_dPres[numFluxSupportPoints] = { stack.dTrans_dPres[connectionIndex][0],
                                                            stack.dTrans_dPres[connectionIndex][1] };

        //***** calculation of flux *****
        // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {
          // create local work arrays
          real64 potGrad = 0.0;
          real64 phaseFlux = 0.0;
          real64 dPhaseFlux_dP[numFluxSupportPoints]{};
          real64 dPhaseFlux_dC[numFluxSupportPoints][numComp]{};

          localIndex k_up = -1;


          isothermalCompositionalMultiphaseFVMKernelUtilities::PPUPhaseFlux::compute< numComp, numFluxSupportPoints >
            ( m_numPhases,
            ip,
            m_kernelFlags.isSet( isothermalCompositionalMultiphaseFVMKernels::KernelFlags::CapPressure ),
            seri, sesri, sei,
            trans,
            dTrans_dPres,
            m_pres,
            m_gravCoef,
            m_phaseMob, m_dPhaseMob,
            m_phaseVolFrac,
            m_dPhaseVolFrac,
            m_phaseCompFrac, m_dPhaseCompFrac,
            m_dCompFrac_dCompDens,
            m_phaseMassDens, m_dPhaseMassDens,
            m_phaseCapPressure, m_dPhaseCapPressure_dPhaseVolFrac,
            k_up,
            potGrad,
            phaseFlux,
            dPhaseFlux_dP,
            dPhaseFlux_dC,
            compFlux,
            dCompFlux_dP,
            dCompFlux_dC,
            dCompFlux_dTrans );

          // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
          // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
          compFluxKernelOp( ip, k, seri, sesri, sei, connectionIndex,
                            k_up, seri[k_up], sesri[k_up], sei[k_up], potGrad,
                            phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

        } // loop over phases

        /// populate local flux vector and derivatives
        for( integer ic = 0; ic < numComp; ++ic )
        {
          integer const eqIndex0 = k[0] * numEqn + ic;
          integer const eqIndex1 = k[1] * numEqn + ic;

          stack.localFlux[eqIndex0]  +=  m_dt * compFlux[ic];
          stack.localFlux[eqIndex1]  -=  m_dt * compFlux[ic];

          for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
          {
            localIndex const localDofIndexPres = k[ke] * numDof;
            stack.localFluxJacobian[eqIndex0][localDofIndexPres] += m_dt * dCompFlux_dP[ke][ic];
            stack.localFluxJacobian[eqIndex1][localDofIndexPres] -= m_dt * dCompFlux_dP[ke][ic];

            for( integer jc = 0; jc < numComp; ++jc )
            {
              localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
              stack.localFluxJacobian[eqIndex0][localDofIndexComp] += m_dt * dCompFlux_dC[ke][ic][jc];
              stack.localFluxJacobian[eqIndex1][localDofIndexComp] -= m_dt * dCompFlux_dC[ke][ic][jc];
            }
          }
        }

        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 dFlux_dAper[2] = {0.0, 0.0};
          dFlux_dAper[0] =  m_dt * dCompFlux_dTrans[ic] * stack.dTrans_dDispJump[connectionIndex][0][0];
          dFlux_dAper[1] = -m_dt * dCompFlux_dTrans[ic] * stack.dTrans_dDispJump[connectionIndex][1][0];

          stack.dFlux_dAperture[ic][k[0]][k[0]] += dFlux_dAper[0];
          stack.dFlux_dAperture[ic][k[0]][k[1]] += dFlux_dAper[1];
          stack.dFlux_dAperture[ic][k[1]][k[0]] -= dFlux_dAper[0];
          stack.dFlux_dAperture[ic][k[1]][k[1]] -= dFlux_dAper[1];
        }
//
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
      // TODO
      localIndex const ei = LvArray::integerConversion< localIndex >( m_sei( iconn, i ) );
      for( integer ic = 0; ic < numComp; ++ic )
      {
        localIndex const row = ei * numComp + ic;

        m_dR_dAper.addToRowBinarySearch< parallelDeviceAtomic >( row,
                                                                 stack.localColIndices.data(),
                                                                 stack.dFlux_dAperture[ic][i].dataIfContiguous(),
                                                                 stack.stencilSize );
      }

      // call the lambda to assemble additional terms, such as thermal terms
      kernelOp( i, localRow );
    } );
  }

private:

  CRSMatrixView< real64, localIndex const > m_dR_dAper;

  ElementViewConst< arrayView4d< real64 const > > const m_dPerm_dDispJump;
};



/**
 * @class FluxComputeKernelFactory
 */
class FluxComputeKernelFactory
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
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasCapPressure,
                   integer const useTotalMassEquation,
                   UpwindingParameters GEOS_UNUSED_PARAM( upwindingParams ),
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   SurfaceElementStencilWrapper const & stencilWrapper,
                   real64 const dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   CRSMatrixView< real64, localIndex const > const & dR_dAper )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      BitFlags< isothermalCompositionalMultiphaseFVMKernels::KernelFlags > kernelFlags;
      if( hasCapPressure )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::KernelFlags::CapPressure );
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::KernelFlags::TotalMassEquation );

      using kernelType = FluxComputeKernel< NUM_COMP, NUM_DOF >;
      typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename kernelType::FracturePermeabilityAccessors fracPermAccessors( elemManager, solverName );

      kernelType kernel( numPhases, rankOffset, stencilWrapper, dofNumberAccessor,
                         compFlowAccessors, multiFluidAccessors, capPressureAccessors, permeabilityAccessors, fracPermAccessors,
                         dt, localMatrix, localRhs, kernelFlags, dR_dAper );

      kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace multiphasePoromechanicsConformingFracturesKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
