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

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIRICHLETFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIRICHLETFLUXCOMPUTEKERNEL_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"


namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

/******************************** DirichletFluxComputeKernel ********************************/

/**
 * @class DirichletFluxComputeKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam FLUIDWRAPPER the type of the fluid wrapper
 * @brief Define the interface for the assembly kernel in charge of Dirichlet face flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename FLUIDWRAPPER >
class DirichletFluxComputeKernel : public FluxComputeKernel< NUM_COMP,
                                                             NUM_DOF,
                                                             BoundaryStencilWrapper >
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

  using AbstractBase::m_dt;
  using AbstractBase::m_numPhases;
  using AbstractBase::m_rankOffset;
  using AbstractBase::m_dofNumber;
  using AbstractBase::m_ghostRank;
  using AbstractBase::m_gravCoef;
  using AbstractBase::m_pres;
  using AbstractBase::m_phaseCompFrac;
  using AbstractBase::m_dPhaseCompFrac;
  using AbstractBase::m_dCompFrac_dCompDens;
  using AbstractBase::m_localMatrix;
  using AbstractBase::m_localRhs;
  using AbstractBase::m_kernelFlags;

  using Base = isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernel< NUM_COMP, NUM_DOF, BoundaryStencilWrapper >;
  using Base::numComp;
  using Base::numDof;
  using Base::numEqn;
  using Base::m_stencilWrapper;
  using Base::m_phaseMob;
  using Base::m_dPhaseMob;
  using Base::m_phaseMassDens;
  using Base::m_dPhaseMassDens;
  using Base::m_permeability;
  using Base::m_dPerm_dPres;
  using Base::m_seri;
  using Base::m_sesri;
  using Base::m_sei;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] faceManager the face manager
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] fluidWrapper reference to the fluid wrapper
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
  DirichletFluxComputeKernel( integer const numPhases,
                              globalIndex const rankOffset,
                              FaceManager const & faceManager,
                              BoundaryStencilWrapper const & stencilWrapper,
                              FLUIDWRAPPER const & fluidWrapper,
                              DofNumberAccessor const & dofNumberAccessor,
                              CompFlowAccessors const & compFlowAccessors,
                              MultiFluidAccessors const & multiFluidAccessors,
                              CapPressureAccessors const & capPressureAccessors,
                              PermeabilityAccessors const & permeabilityAccessors,
                              real64 const dt,
                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                              arrayView1d< real64 > const & localRhs,
                              BitFlags< KernelFlags > kernelFlags )
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
    m_facePres( faceManager.getField< fields::flow::facePressure >() ),
    m_faceTemp( faceManager.getField< fields::flow::faceTemperature >() ),
    m_faceCompFrac( faceManager.getField< fields::flow::faceGlobalCompFraction >() ),
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
                    localIndex GEOS_UNUSED_PARAM( numElems )) {}

    // Transmissibility
    real64 transmissibility = 0.0;

    // Component fluxes and derivatives

    /// Component fluxes
    real64 compFlux[numComp]{};
    /// Derivatives of component fluxes wrt pressure
    real64 dCompFlux_dP[numComp]{};
    /// Derivatives of component fluxes wrt component densities
    real64 dCompFlux_dC[numComp][numComp]{};

    // Local degrees of freedom and local residual/jacobian

    /// Indices of the matrix rows/columns corresponding to the dofs in this face
    globalIndex dofColIndices[numDof]{};

    /// Storage for the face local residual vector
    real64 localFlux[numEqn]{};
    /// Storage for the face local Jacobian matrix
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
    using Deriv = constitutive::multifluid::DerivativeOffset;
    using Order = BoundaryStencil::Order;

    localIndex const er  = m_seri( iconn, Order::ELEM );
    localIndex const esr = m_sesri( iconn, Order::ELEM );
    localIndex const ei  = m_sei( iconn, Order::ELEM );
    localIndex const kf  = m_sei( iconn, Order::FACE );

    // Step 1: compute the transmissibility at the boundary face

    real64 dTrans_dPerm[3]{};
    m_stencilWrapper.computeWeights( iconn,
                                     m_permeability,
                                     stack.transmissibility,
                                     dTrans_dPerm );
    real64 const dTrans_dPres = LvArray::tensorOps::AiBi< 3 >( dTrans_dPerm, m_dPerm_dPres[er][esr][ei][0] );

    // Step 2: compute the fluid properties on the face
    // This is needed to get the phase mass density and the phase comp fraction at the face
    // Because we approximate the face mobility using the total element mobility

    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseFrac( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseDens( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseMassDens( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseVisc( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseEnthalpy( 1, 1, m_numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > facePhaseInternalEnergy( 1, 1, m_numPhases );
    StackArray< real64, 4, constitutive::MultiFluidBase::MAX_NUM_PHASES * NUM_COMP,
                constitutive::multifluid::LAYOUT_PHASE_COMP > facePhaseCompFrac( 1, 1, m_numPhases, NUM_COMP );
    real64 faceTotalDens = 0.0;

    constitutive::MultiFluidBase::KernelWrapper::computeValues( m_fluidWrapper,
                                                                m_facePres[kf],
                                                                m_faceTemp[kf],
                                                                m_faceCompFrac[kf],
                                                                facePhaseFrac[0][0],
                                                                facePhaseDens[0][0],
                                                                facePhaseMassDens[0][0],
                                                                facePhaseVisc[0][0],
                                                                facePhaseEnthalpy[0][0],
                                                                facePhaseInternalEnergy[0][0],
                                                                facePhaseCompFrac[0][0],
                                                                faceTotalDens );

    // Step 3: loop over phases, compute and upwind phase flux and sum contributions to each component's flux

    for( integer ip = 0; ip < m_numPhases; ++ip )
    {

      // working variables
      real64 dDensMean_dC[numComp]{};
      real64 dF_dC[numComp]{};
      real64 dProp_dC[numComp]{};

      real64 phaseFlux = 0.0; // for the lambda
      real64 dPhaseFlux_dP = 0.0;
      real64 dPhaseFlux_dC[numComp]{};


      // Step 3.1: compute the average phase mass density at the face

      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er][esr][ei],
                      m_dPhaseMassDens[er][esr][ei][0][ip],
                      dProp_dC,
                      Deriv::dC );

      // average density and derivatives
      real64 const densMean = 0.5 * ( m_phaseMassDens[er][esr][ei][0][ip] + facePhaseMassDens[0][0][ip] );
      real64 const dDensMean_dP = 0.5 * m_dPhaseMassDens[er][esr][ei][0][ip][Deriv::dP];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dDensMean_dC[jc] = 0.5 * dProp_dC[jc];
      }


      // Step 3.2: compute the (TPFA) potential difference at the face

      real64 const gravTimesDz = m_gravCoef[er][esr][ei] - m_faceGravCoef[kf];
      real64 const potDif = m_pres[er][esr][ei] - m_facePres[kf] - densMean * gravTimesDz;
      real64 const f = stack.transmissibility * potDif;
      real64 const dF_dP = stack.transmissibility * ( 1.0 - dDensMean_dP * gravTimesDz ) + dTrans_dPres * potDif;
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dF_dC[jc] = -stack.transmissibility * dDensMean_dC[jc] * gravTimesDz;
      }

      // Step 3.3: computation of the mobility
      // We do that before the if/else statement to be able to pass it to the compFluxOpKernel

      // recomputing the exact mobility at the face would be quite complex, as it would require:
      //   1) computing the saturation
      //   2) computing the relperm
      //   3) computing the mobility as \lambda_p = \rho_p kr_p( S_p ) / \mu_p
      // the second step in particular would require yet another dispatch to get the relperm model
      // so, for simplicity, we approximate the face mobility as
      //    \lambda^approx_p = \rho_p S_p / \mu_p
      //                     = \rho_p ( (nu_p / rho_p) * rho_t ) / \mu_p (plugging the expression of saturation)
      //                     = \nu_p * rho_t / \mu_p
      // fortunately, we don't need the derivatives
      real64 const facePhaseMob = ( facePhaseFrac[0][0][ip] > 0.0 )
  ? facePhaseFrac[0][0][ip] * faceTotalDens / facePhaseVisc[0][0][ip]
  : 0.0;

      // *** upwinding ***
      // Step 3.4: upwinding based on the sign of the phase potential gradient
      // It is easier to hard-code the if/else because it is difficult to address elem and face variables in a uniform way

      if( potDif >= 0 ) // the element is upstream
      {

        // compute the phase flux and derivatives using the element mobility
        phaseFlux = m_phaseMob[er][esr][ei][ip] * f;
        dPhaseFlux_dP = m_phaseMob[er][esr][ei][ip] * dF_dP + m_dPhaseMob[er][esr][ei][ip][Deriv::dP] * f;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[jc] =
            m_phaseMob[er][esr][ei][ip] * dF_dC[jc] + m_dPhaseMob[er][esr][ei][ip][Deriv::dC+jc] * f;
        }

        // slice some constitutive arrays to avoid too much indexing in component loop
        arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE_COMP-3 > phaseCompFracSub =
          m_phaseCompFrac[er][esr][ei][0][ip];
        arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC-3 > dPhaseCompFracSub =
          m_dPhaseCompFrac[er][esr][ei][0][ip];

        // compute component fluxes and derivatives using element composition
        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const ycp = phaseCompFracSub[ic];
          stack.compFlux[ic] += phaseFlux * ycp;
          stack.dCompFlux_dP[ic] += dPhaseFlux_dP * ycp + phaseFlux * dPhaseCompFracSub[ic][Deriv::dP];

          applyChainRule( numComp,
                          m_dCompFrac_dCompDens[er][esr][ei],
                          dPhaseCompFracSub[ic],
                          dProp_dC,
                          Deriv::dC );
          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dCompFlux_dC[ic][jc] += dPhaseFlux_dC[jc] * ycp + phaseFlux * dProp_dC[jc];
          }
        }

      }
      else // the face is upstream
      {

        // compute the phase flux and derivatives using the approximated face mobility
        // we only have to take derivatives of the potential gradient in this case
        phaseFlux = facePhaseMob * f;
        dPhaseFlux_dP = facePhaseMob * dF_dP;
        for( integer jc = 0; jc < numComp; ++jc )
        {
          dPhaseFlux_dC[jc] = facePhaseMob * dF_dC[jc];
        }

        // compute component fluxes and derivatives using the face composition
        for( integer ic = 0; ic < numComp; ++ic )
        {
          real64 const ycp = facePhaseCompFrac[0][0][ip][ic];
          stack.compFlux[ic] += phaseFlux * ycp;
          stack.dCompFlux_dP[ic] += dPhaseFlux_dP * ycp;
          for( integer jc = 0; jc < numComp; ++jc )
          {
            stack.dCompFlux_dC[ic][jc] += dPhaseFlux_dC[jc] * ycp;
          }
        }
      }

      // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
      // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
      compFluxKernelOp( ip, er, esr, ei, kf, f,
                        facePhaseMob, facePhaseEnthalpy[0][0], facePhaseCompFrac[0][0],
                        phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC );

    }

    // *** end of upwinding

    // Step 4: populate local flux vector and derivatives
    for( integer ic = 0; ic < numComp; ++ic )
    {
      stack.localFlux[ic]            = m_dt * stack.compFlux[ic];
      stack.localFluxJacobian[ic][0] = m_dt * stack.dCompFlux_dP[ic];
      for( integer jc = 0; jc < numComp; ++jc )
      {
        stack.localFluxJacobian[ic][jc+1] = m_dt * stack.dCompFlux_dC[ic][jc];
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
                 FUNC && assemblyKernelOp = NoOpFunc{} ) const
  {
    using namespace compositionalMultiphaseUtilities;
    using Order = BoundaryStencil::Order;

    if( AbstractBase::m_kernelFlags.isSet( KernelFlags::TotalMassEquation ) )
    {
      // Apply equation/variable change transformation(s)
      real64 work[numDof]{};
      shiftRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numDof, stack.localFluxJacobian, work );
      shiftElementsAheadByOneAndReplaceFirstElementWithSum( numComp, stack.localFlux );
    }

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    if( m_ghostRank[m_seri( iconn, Order::ELEM )][m_sesri( iconn, Order::ELEM )][m_sei( iconn, Order::ELEM )] < 0 )
    {
      globalIndex const globalRow = m_dofNumber[m_seri( iconn, Order::ELEM )][m_sesri( iconn, Order::ELEM )][m_sei( iconn, Order::ELEM )];
      localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
      GEOS_ASSERT_GE( localRow, 0 );
      GEOS_ASSERT_GT( AbstractBase::m_localMatrix.numRows(), localRow + numComp );

      for( integer ic = 0; ic < numComp; ++ic )
      {
        RAJA::atomicAdd( parallelDeviceAtomic{}, &AbstractBase::m_localRhs[localRow + ic], stack.localFlux[ic] );
        AbstractBase::m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
          ( localRow + ic,
          stack.dofColIndices,
          stack.localFluxJacobian[ic],
          numDof );
      }

      // call the lambda to assemble additional terms, such as thermal terms
      assemblyKernelOp( localRow );
    }
  }

protected:

  /// Views on face pressure, temperature, and composition
  arrayView1d< real64 const > const m_facePres;
  arrayView1d< real64 const > const m_faceTemp;
  arrayView2d< real64 const, compflow::USD_COMP > const m_faceCompFrac;

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
   * @param[in] numComps the number of fluid components
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey string to get the element degrees of freedom numbers
   * @param[in] solverName name of the solver (to name accessors)
   * @param[in] faceManager reference to the face manager
   * @param[in] elemManager reference to the element region manager
   * @param[in] stencilWrapper reference to the boundary stencil wrapper
   * @param[in] fluidBase the multifluid constitutive model
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numComps,
                   integer const numPhases,
                   globalIndex const rankOffset,
                   integer const useTotalMassEquation,
                   string const & dofKey,
                   string const & solverName,
                   FaceManager const & faceManager,
                   ElementRegionManager const & elemManager,
                   BoundaryStencilWrapper const & stencilWrapper,
                   constitutive::MultiFluidBase & fluidBase,
                   real64 const dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    constitutive::constitutiveComponentUpdatePassThru( fluidBase, numComps, [&]( auto & fluid, auto NC )
    {
      using FluidType = TYPEOFREF( fluid );
      typename FluidType::KernelWrapper const fluidWrapper = fluid.createKernelWrapper();

      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      // for now, we neglect capillary pressure in the kernel
      BitFlags< KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( KernelFlags::TotalMassEquation );

      using kernelType = DirichletFluxComputeKernel< NUM_COMP, NUM_DOF, typename FluidType::KernelWrapper >;
      typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
      typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

      kernelType kernel( numPhases, rankOffset, faceManager, stencilWrapper, fluidWrapper,
                         dofNumberAccessor, compFlowAccessors, multiFluidAccessors, capPressureAccessors, permeabilityAccessors,
                         dt, localMatrix, localRhs, kernelFlags );
      kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
    } );
  }
};

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIRICHLETFLUXCOMPUTEKERNEL_HPP
