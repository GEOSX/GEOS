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
 * @file DissipationCompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_DISSIPATIONCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_DISSIPATIONCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "IsothermalCompositionalMultiphaseFVMKernels.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"

namespace geos
{

namespace dissipationCompositionalMultiphaseFVMKernels
{

static constexpr integer newtonContinuationCutoffIteration = 5;
static constexpr real64 initialDirectionalCoef = 100;
static constexpr real64 multiplierDirectionalCoef = 1000;

/******************************** FaceBasedAssemblyKernel ********************************/

/**
 * @class FaceBasedAssemblyKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class FaceBasedAssemblyKernel : public isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >
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

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelBase;
  using DofNumberAccessor = AbstractBase::DofNumberAccessor;
  using CompFlowAccessors = AbstractBase::CompFlowAccessors;
  using MultiFluidAccessors = AbstractBase::MultiFluidAccessors;
  using CapPressureAccessors = AbstractBase::CapPressureAccessors;
  using PermeabilityAccessors = AbstractBase::PermeabilityAccessors;

  using AbstractBase::m_dt;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_dCompFrac_dCompDens;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::numFluxSupportPoints;
  using Base::m_permeability;
  using Base::m_dPerm_dPres;

  using DissCompFlowAccessors =
    StencilAccessors< fields::flow::pressure_n,
                      fields::flow::globalCompDensity,
                      fields::flow::globalCompFraction,
                      fields::elementVolume >;

  using PorosityAccessors =
    StencilMaterialAccessors< constitutive::PorosityBase, fields::porosity::porosity_n >;

  using Deriv = constitutive::multifluid::DerivativeOffset;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor accessor for the dofs numbers
   * @param[in] compFlowAccessor accessor for wrappers registered by the solver
   * @param[in] dissCompFlowAccessor accessor for wrappers registered by the solver needed for dissipation
   * @param[in] multiFluidAccessor accessor for wrappers registered by the multifluid model
   * @param[in] capPressureAccessors accessor for wrappers registered by the cap pressure model
   * @param[in] permeabilityAccessors accessor for wrappers registered by the permeability model
   * @param[in] porosityAccessors accessor for wrappers registered by the porosity model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] omega omega coefficient for DBC
   * @param[in] curNewton current Newton iteration number
   * @param[in] continuation flag indicating if continuation is used or not
   * @param[in] miscible flag to trigger some treatment for miscible cases
   * @param[in] kappamin minimum value for kappa coefficient in DBC
   * @param[in] contMultiplier continuation multiplier factor (should be < 1)
   */
  FaceBasedAssemblyKernel( integer const numPhases,
                           globalIndex const rankOffset,
                           STENCILWRAPPER const & stencilWrapper,
                           DofNumberAccessor const & dofNumberAccessor,
                           CompFlowAccessors const & compFlowAccessors,
                           DissCompFlowAccessors const & dissCompFlowAccessors,
                           MultiFluidAccessors const & multiFluidAccessors,
                           CapPressureAccessors const & capPressureAccessors,
                           PermeabilityAccessors const & permeabilityAccessors,
                           PorosityAccessors const & porosityAccessors,
                           real64 const & dt,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs,
                           BitFlags< isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags > kernelFlags,
                           real64 const omega,
                           integer const curNewton,
                           integer const continuation,
                           integer const miscible,
                           real64 const kappamin,
                           real64 const contMultiplier )
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
    m_pres_n( dissCompFlowAccessors.get( fields::flow::pressure_n {} ) ),
    m_porosity_n( porosityAccessors.get( fields::porosity::porosity_n {} ) ),
    m_volume( dissCompFlowAccessors.get( fields::elementVolume {} ) ),
    m_compFrac( dissCompFlowAccessors.get( fields::flow::globalCompFraction {} ) ),
    m_omegaDBC( omega ),
    m_miscibleDBC( miscible )
  {
    /// Step 1. Calculate the continuation parameter based on the current Newton iteration
    m_kappaDBC = 1.0;   // default value

    if( continuation )   // if continuation is enabled
    {
      if( curNewton >= newtonContinuationCutoffIteration )
      {
        m_kappaDBC = kappamin;
      }
      else
      {
        for( int mp = 0; mp < curNewton; mp++ )
        {
          m_kappaDBC *= contMultiplier;
        }
        m_kappaDBC = std::max( m_kappaDBC, kappamin );
      }
    }
  }

  /**
   * @brief Compute the local flux contributions to the residual and Jacobian, including dissipation
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFlux( localIndex const iconn,
                    typename Base::StackVariables & stack ) const
  {
    // ***********************************************
    // First, we call the base computeFlux to compute:
    //  1) compFlux and its derivatives,
    //
    // We use the lambda below (called **inside** the phase loop of the base computeFlux) to compute dissipation terms
    Base::computeFlux( iconn, stack, [&] ( integer const ip,
                                           localIndex const (&k)[2],
                                           localIndex const (&seri)[2],
                                           localIndex const (&sesri)[2],
                                           localIndex const (&sei)[2],
                                           localIndex const connectionIndex,
                                           localIndex const k_up,
                                           localIndex const er_up,
                                           localIndex const esr_up,
                                           localIndex const ei_up,
                                           real64 const & potGrad,
                                           real64 const & phaseFlux,
                                           real64 const (&dPhaseFlux_dP)[2],
                                           real64 const (&dPhaseFlux_dC)[2][numComp] )
    {
      GEOS_UNUSED_VAR( k_up, potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, er_up, esr_up, ei_up );

      /// Storing dissipation flux and its derivatives locally
      real64 dissFlux[numComp]{};
      real64 dDissFlux_dP[numFluxSupportPoints][numComp]{};
      real64 dDissFlux_dC[numFluxSupportPoints][numComp][numComp]{};
      real64 fluxPointCoef[numFluxSupportPoints] = {1.0, -1.0}; // for gradients

      real64 viscosityMult[3] = {1.0, 1.0, 1.0}; // for viscosity

      /// Step 2. Collect all contributions
      real64 poreVolume_n = 0.0; // Pore volume contribution
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        poreVolume_n += m_volume[seri[ke]][sesri[ke]][sei[ke]] * m_porosity_n[seri[ke]][sesri[ke]][sei[ke]][0];
      }
      poreVolume_n /= numFluxSupportPoints;

      // potential gradient contribution
      // pressure
      real64 pressure_gradient = 0;
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        pressure_gradient += fluxPointCoef[ke] * m_pres_n[seri[ke]][sesri[ke]][sei[ke]];

      real64 const potential_gradient = LvArray::math::abs( pressure_gradient );

      real64 grad_depth = 0;
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        grad_depth += fluxPointCoef[ke] * m_gravCoef[seri[ke]][sesri[ke]][sei[ke]];
      }

      // bias towards x and y direction for miscible
      /* The point of this "directional coefficient" was to add less dissipation flux in the direction of gravity and
         more in all other directions.  These hard-coded values of 1000 and 100 were picked by some tuning and as
         expected are not robust. Further research/testing is needed to generalize this idea.
       */
      real64 directional_coef = 1.0;
      if( m_miscibleDBC )
      {
        directional_coef = initialDirectionalCoef;
        if( LvArray::math::abs( grad_depth ) > 0.0 )
        {
          real64 const d2 = LvArray::math::abs( grad_depth * grad_depth );
          if( multiplierDirectionalCoef / d2 < initialDirectionalCoef )
            directional_coef = multiplierDirectionalCoef / d2;
        }
      }

      // multiplier with all contributions
      real64 const multiplier = m_kappaDBC * m_omegaDBC / poreVolume_n * m_dt * potential_gradient * directional_coef;

      /// Step 3. Compute the dissipation flux and its derivative
      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        real64 const coef = multiplier * viscosityMult[ip] * stack.transmissibility[connectionIndex][ke];
        for( integer ic = 0; ic < numComp; ++ic )
        {
          localIndex const er  = seri[ke];
          localIndex const esr = sesri[ke];
          localIndex const ei  = sei[ke];

          // composition gradient contribution to the dissipation flux
          dissFlux[ic] += coef * m_compFrac[er][esr][ei][ic];

          dDissFlux_dP[ke][ic] = 0; // no dependency on pressure at n+1, compFrad is global component fraction (z_c)

          for( integer jc = 0; jc < numComp; ++jc )
          {
            // composition gradient derivative with respect to component density contribution to the dissipation flux
            dDissFlux_dC[ke][ic][jc] += coef * m_dCompFrac_dCompDens[er][esr][ei][ic][jc];
          }
        }
      }

      /// Step 4: add the dissipation flux and its derivatives to the residual and Jacobian
      for( integer ic = 0; ic < numComp; ++ic )
      {
        integer const eqIndex0 = k[0] * numEqn + ic;
        integer const eqIndex1 = k[1] * numEqn + ic;

        stack.localFlux[eqIndex0] += m_dt * dissFlux[ic];
        stack.localFlux[eqIndex1] -= m_dt * dissFlux[ic];

        for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
        {
          localIndex const localDofIndexPres = k[ke] * numDof;
          stack.localFluxJacobian[eqIndex0][localDofIndexPres] +=  m_dt * dDissFlux_dP[ke][ic];
          stack.localFluxJacobian[eqIndex1][localDofIndexPres] -=  m_dt * dDissFlux_dP[ke][ic];

          for( integer jc = 0; jc < numComp; ++jc )
          {
            localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
            stack.localFluxJacobian[eqIndex0][localDofIndexComp] += m_dt * dDissFlux_dC[ke][ic][jc];
            stack.localFluxJacobian[eqIndex1][localDofIndexComp] -= m_dt * dDissFlux_dC[ke][ic][jc];
          }
        }
      }

    } ); // end call to Base::computeFlux

  }

protected:

  /// Views on flow properties at the previous converged time step
  ElementViewConst< arrayView1d< real64 const > > const m_pres_n;
  ElementViewConst< arrayView2d< real64 const > > const m_porosity_n;

  ElementViewConst< arrayView1d< real64 const > > const m_volume;
  ElementViewConst< arrayView2d< real64 const, compflow::USD_COMP > > const m_compFrac;

  // DBC specific parameters
  // input
  real64 m_omegaDBC;
  integer m_miscibleDBC;
  // computed
  real64 m_kappaDBC;
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
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] hasCapPressure flag specifying whether capillary pressure is used or not
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY, typename STENCILWRAPPER >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   string const & dofKey,
                   integer const hasCapPressure,
                   integer const useTotalMassEquation,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const & dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs,
                   real64 const omega,
                   integer const curNewton,
                   integer const continuation,
                   integer const miscible,
                   real64 const kappamin,
                   real64 const contMultiplier )
  {
    isothermalCompositionalMultiphaseBaseKernels::
      internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      BitFlags< isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags > kernelFlags;
      if( hasCapPressure )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags::CapPressure );
      if( useTotalMassEquation )
        kernelFlags.set( isothermalCompositionalMultiphaseFVMKernels::FaceBasedAssemblyKernelFlags::TotalMassEquation );

      using KERNEL_TYPE = FaceBasedAssemblyKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename KERNEL_TYPE::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename KERNEL_TYPE::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename KERNEL_TYPE::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
      typename KERNEL_TYPE::PorosityAccessors porosityAccessors( elemManager, solverName );
      typename KERNEL_TYPE::DissCompFlowAccessors dissCompFlowAccessors( elemManager, solverName );

      KERNEL_TYPE kernel( numPhases, rankOffset, stencilWrapper, dofNumberAccessor, compFlowAccessors, dissCompFlowAccessors,
                          multiFluidAccessors, capPressureAccessors, permeabilityAccessors, porosityAccessors,
                          dt, localMatrix, localRhs, kernelFlags, omega, curNewton, continuation, miscible, kappamin, contMultiplier );
      KERNEL_TYPE::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace DissipationCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_DISSIPATIONCOMPOSITIONALMULTIPHASEFVMKERNELS_HPP
