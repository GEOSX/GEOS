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
 * @file DiffusionDispersionFluxComputeKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIFFUSIONDISPERSIONFLUXCOMPUTEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIFFUSIONDISPERSIONFLUXCOMPUTEKERNEL_HPP

#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/diffusion/DiffusionFields.hpp"
#include "constitutive/diffusion/DiffusionBase.hpp"
#include "constitutive/dispersion/DispersionFields.hpp"
#include "constitutive/dispersion/DispersionBase.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/KernelLaunchSelectors.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

/******************************** DiffusionDispersionFluxComputeKernel ********************************/

/**
 * @class DiffusionDispersionFluxComputeKernel
 * @tparam NUM_COMP number of fluid components
 * @tparam NUM_DOF number of degrees of freedom
 * @tparam STENCILWRAPPER the type of the stencil wrapper
 * @brief Define the interface for the assembly kernel in charge of diffusion/dispersion flux terms
 */
template< integer NUM_COMP, integer NUM_DOF, typename STENCILWRAPPER >
class DiffusionDispersionFluxComputeKernel : public FluxComputeKernelBase
{
public:

  /// Compile time value for the number of components
  static constexpr integer numComp = NUM_COMP;

  /// Compute time value for the number of degrees of freedom
  static constexpr integer numDof = NUM_DOF;

  /// Compute time value for the number of equations (all of them, except the volume balance equation)
  static constexpr integer numEqn = NUM_DOF-1;

  /// Maximum number of elements at the face
  static constexpr localIndex maxNumElems = STENCILWRAPPER::maxNumPointsInFlux;

  /// Maximum number of connections at the face
  static constexpr localIndex maxNumConns = STENCILWRAPPER::maxNumConnections;

  /// Maximum number of points in the stencil
  static constexpr localIndex maxStencilSize = STENCILWRAPPER::maxStencilSize;

  /// Number of flux support points (hard-coded for TFPA)
  static constexpr integer numFluxSupportPoints = 2;

  using AbstractBase = isothermalCompositionalMultiphaseFVMKernels::FluxComputeKernelBase;
  using AbstractBase::m_dPhaseVolFrac;
  using AbstractBase::m_kernelFlags;

  using DiffusionAccessors =
    StencilMaterialAccessors< constitutive::DiffusionBase,
                              fields::diffusion::diffusivity,
                              fields::diffusion::dDiffusivity_dTemperature,
                              fields::diffusion::phaseDiffusivityMultiplier >;

  using DispersionAccessors =
    StencilMaterialAccessors< constitutive::DispersionBase,
                              fields::dispersion::dispersivity >;

  using PorosityAccessors =
    StencilMaterialAccessors< constitutive::PorosityBase,
                              fields::porosity::referencePorosity >;

  /**
   * @brief Constructor for the kernel interface
   * @param[in] numPhases the number of fluid phases
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] stencilWrapper reference to the stencil wrapper
   * @param[in] dofNumberAccessor
   * @param[in] compFlowAccessors
   * @param[in] multiFluidAccessors
   * @param[in] diffusionAccessors
   * @param[in] dispersionAccessors
   * @param[in] porosityAccessors
   * @param[in] dt time step size
   * @param[inout] localMatrix the local CRS matrix
   * @param[inout] localRhs the local right-hand side vector
   * @param[in] kernelFlags flags packed together
   */
  DiffusionDispersionFluxComputeKernel( integer const numPhases,
                                        globalIndex const rankOffset,
                                        STENCILWRAPPER const & stencilWrapper,
                                        DofNumberAccessor const & dofNumberAccessor,
                                        CompFlowAccessors const & compFlowAccessors,
                                        MultiFluidAccessors const & multiFluidAccessors,
                                        DiffusionAccessors const & diffusionAccessors,
                                        DispersionAccessors const & dispersionAccessors,
                                        PorosityAccessors const & porosityAccessors,
                                        real64 const dt,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs,
                                        BitFlags< KernelFlags > kernelFlags )
    : FluxComputeKernelBase( numPhases,
                             rankOffset,
                             dofNumberAccessor,
                             compFlowAccessors,
                             multiFluidAccessors,
                             dt,
                             localMatrix,
                             localRhs,
                             kernelFlags ),
    m_phaseVolFrac( compFlowAccessors.get( fields::flow::phaseVolumeFraction {} ) ),
    m_phaseDens( multiFluidAccessors.get( fields::multifluid::phaseDensity {} ) ),
    m_dPhaseDens( multiFluidAccessors.get( fields::multifluid::dPhaseDensity {} ) ),
    m_diffusivity( diffusionAccessors.get( fields::diffusion::diffusivity {} ) ),
    m_dDiffusivity_dTemp( diffusionAccessors.get( fields::diffusion::dDiffusivity_dTemperature {} ) ),
    m_phaseDiffusivityMultiplier( diffusionAccessors.get( fields::diffusion::phaseDiffusivityMultiplier {} ) ),
    m_dispersivity( dispersionAccessors.get( fields::dispersion::dispersivity {} ) ),
    m_referencePorosity( porosityAccessors.get( fields::porosity::referencePorosity {} ) ),
    m_stencilWrapper( stencilWrapper ),
    m_seri( stencilWrapper.getElementRegionIndices() ),
    m_sesri( stencilWrapper.getElementSubRegionIndices() ),
    m_sei( stencilWrapper.getElementIndices() )
  { }

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
      numConnectedElems( numElems ),
      dofColIndices( size * numDof ),
      localFlux( numElems * numEqn ),
      localFluxJacobian( numElems * numEqn, size * numDof )
    {}

    // Stencil information

    /// Stencil size for a given connection
    localIndex const stencilSize;
    /// Number of elements connected at a given connection
    localIndex const numConnectedElems;

    /// Transmissibility
    real64 transmissibility[maxNumConns][numFluxSupportPoints]{};
    /// Derivatives of transmissibility with respect to pressure
    real64 dTrans_dTemp[maxNumConns][numFluxSupportPoints]{};

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
  inline
  localIndex stencilSize( localIndex const iconn ) const
  { return m_sei[iconn].size(); }

  /**
   * @brief Getter for the number of elements at this connection
   * @param[in] iconn the connection index
   * @return the number of elements at this connection
   */
  GEOS_HOST_DEVICE
  inline
  localIndex numPointsInFlux( localIndex const iconn ) const
  { return m_stencilWrapper.numPointsInFlux( iconn ); }


  /**
   * @brief Performs the setup phase for the kernel.
   * @param[in] iconn the connection index
   * @param[in] stack the stack variables
   */
  GEOS_HOST_DEVICE
  inline
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
   * @brief Compute the local diffusion flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] diffusionFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void computeDiffusionFlux( localIndex const iconn,
                             StackVariables & stack,
                             FUNC && diffusionFluxKernelOp = NoOpFunc{} ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    // first, compute the transmissibilities at this face
    m_stencilWrapper.computeWeights( iconn,
                                     m_diffusivity,
                                     m_dDiffusivity_dTemp,
                                     stack.transmissibility,
                                     stack.dTrans_dTemp );


    localIndex k[numFluxSupportPoints]{};
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
        real64 diffusionFlux[numComp]{};
        real64 dDiffusionFlux_dP[numFluxSupportPoints][numComp]{};
        real64 dDiffusionFlux_dC[numFluxSupportPoints][numComp][numComp]{};
        real64 dDens_dC[numComp]{};

        real64 const trans[numFluxSupportPoints] = { stack.transmissibility[connectionIndex][0],
                                                     stack.transmissibility[connectionIndex][1] };

        //***** calculation of flux *****
        // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {

          // loop over components
          for( integer ic = 0; ic < numComp; ++ic )
          {

            real64 compFracGrad = 0.0;
            real64 dCompFracGrad_dP[numFluxSupportPoints]{};
            real64 dCompFracGrad_dC[numFluxSupportPoints][numComp]{};

            // compute the component fraction gradient using the diffusion transmissibility
            computeCompFractionGradient( ip, ic,
                                         seri, sesri, sei,
                                         trans,
                                         compFracGrad,
                                         dCompFracGrad_dP,
                                         dCompFracGrad_dC );

            // choose upstream cell for composition upwinding
            localIndex const k_up = (compFracGrad >= 0) ? 0 : 1;

            localIndex const er_up  = seri[k_up];
            localIndex const esr_up = sesri[k_up];
            localIndex const ei_up  = sei[k_up];

            // computation of the upwinded mass flux
            real64 const upwindCoefficient =
              m_referencePorosity[er_up][esr_up][ei_up] *
              m_phaseDiffusivityMultiplier[er_up][esr_up][ei_up][0][ip] *
              m_phaseDens[er_up][esr_up][ei_up][0][ip] * m_phaseVolFrac[er_up][esr_up][ei_up][ip];
            diffusionFlux[ic] += upwindCoefficient * compFracGrad;

            // add contributions of the derivatives of component fractions wrt pressure/component fractions
            for( integer ke = 0; ke < numFluxSupportPoints; ke++ )
            {
              dDiffusionFlux_dP[ke][ic] += upwindCoefficient * dCompFracGrad_dP[ke];
              for( integer jc = 0; jc < numComp; ++jc )
              {
                dDiffusionFlux_dC[ke][ic][jc] += upwindCoefficient * dCompFracGrad_dC[ke][jc];
              }
            }

            // add contributions of the derivatives of upwind coefficient wrt pressure/component fractions
            real64 const dUpwindCoefficient_dP =
              m_referencePorosity[er_up][esr_up][ei_up] *
              m_phaseDiffusivityMultiplier[er_up][esr_up][ei_up][0][ip] *
              ( m_dPhaseDens[er_up][esr_up][ei_up][0][ip][Deriv::dP] * m_phaseVolFrac[er_up][esr_up][ei_up][ip]
                + m_phaseDens[er_up][esr_up][ei_up][0][ip] * m_dPhaseVolFrac[er_up][esr_up][ei_up][ip][Deriv::dP] );
            dDiffusionFlux_dP[k_up][ic] += dUpwindCoefficient_dP * compFracGrad;

            applyChainRule( numComp,
                            m_dCompFrac_dCompDens[er_up][esr_up][ei_up],
                            m_dPhaseDens[er_up][esr_up][ei_up][0][ip],
                            dDens_dC,
                            Deriv::dC );
            for( integer jc = 0; jc < numComp; ++jc )
            {
              real64 const dUpwindCoefficient_dC =
                m_referencePorosity[er_up][esr_up][ei_up] *
                m_phaseDiffusivityMultiplier[er_up][esr_up][ei_up][0][ip] *
                ( dDens_dC[jc] * m_phaseVolFrac[er_up][esr_up][ei_up][ip]
                  + m_phaseDens[er_up][esr_up][ei_up][0][ip] * m_dPhaseVolFrac[er_up][esr_up][ei_up][ip][Deriv::dC+jc] );
              dDiffusionFlux_dC[k_up][ic][jc] += dUpwindCoefficient_dC * compFracGrad;
            }

            // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
            // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
            diffusionFluxKernelOp( ip, ic, k, seri, sesri, sei, connectionIndex,
                                   k_up, seri[k_up], sesri[k_up], sei[k_up],
                                   compFracGrad, upwindCoefficient );

          } // loop over components
        } // loop over phases

        // add diffusion flux to local flux and local flux jacobian
        addToLocalFluxAndJacobian( k,
                                   stack,
                                   diffusionFlux,
                                   dDiffusionFlux_dP,
                                   dDiffusionFlux_dC );

        connectionIndex++;
      }   // loop over k[1]
    }   // loop over k[0]
  }

  /**
   * @brief Compute the local dispersion flux contributions to the residual and Jacobian
   * @tparam FUNC the type of the function that can be used to customize the computation of the phase fluxes
   * @param[in] iconn the connection index
   * @param[inout] stack the stack variables
   * @param[in] dispersionFluxKernelOp the function used to customize the computation of the component fluxes
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  inline
  void computeDispersionFlux( localIndex const iconn,
                              StackVariables & stack,
                              FUNC && dispersionFluxKernelOp = NoOpFunc{} ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    // first, compute the transmissibilities at this face
    // note that the dispersion tensor is lagged in iteration
    m_stencilWrapper.computeWeights( iconn,
                                     m_dispersivity,
                                     m_dispersivity, // this is just to pass something, but the resulting derivative won't be used
                                     stack.transmissibility,
                                     stack.dTrans_dTemp ); // will not be used


    localIndex k[numFluxSupportPoints]{};
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
        real64 dispersionFlux[numComp]{};
        real64 dDispersionFlux_dP[numFluxSupportPoints][numComp]{};
        real64 dDispersionFlux_dC[numFluxSupportPoints][numComp][numComp]{};
        real64 dDens_dC[numComp]{};

        real64 const trans[numFluxSupportPoints] = { stack.transmissibility[connectionIndex][0],
                                                     stack.transmissibility[connectionIndex][1] };

        //***** calculation of flux *****
        // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
        for( integer ip = 0; ip < m_numPhases; ++ip )
        {

          // loop over components
          for( integer ic = 0; ic < numComp; ++ic )
          {

            real64 compFracGrad = 0.0;
            real64 dCompFracGrad_dP[numFluxSupportPoints]{};
            real64 dCompFracGrad_dC[numFluxSupportPoints][numComp]{};

            // compute the component fraction gradient using the dispersion transmissibility
            computeCompFractionGradient( ip, ic,
                                         seri, sesri, sei,
                                         trans,
                                         compFracGrad,
                                         dCompFracGrad_dP,
                                         dCompFracGrad_dC );

            // choose upstream cell for composition upwinding
            localIndex const k_up = (compFracGrad >= 0) ? 0 : 1;

            localIndex const er_up  = seri[k_up];
            localIndex const esr_up = sesri[k_up];
            localIndex const ei_up  = sei[k_up];

            // computation of the upwinded mass flux
            dispersionFlux[ic] += m_phaseDens[er_up][esr_up][ei_up][0][ip] * compFracGrad;

            // add contributions of the derivatives of component fractions wrt pressure/component fractions
            for( integer ke = 0; ke < numFluxSupportPoints; ke++ )
            {
              dDispersionFlux_dP[ke][ic] += m_phaseDens[er_up][esr_up][ei_up][0][ip] * dCompFracGrad_dP[ke];
              for( integer jc = 0; jc < numComp; ++jc )
              {
                dDispersionFlux_dC[ke][ic][jc] += m_phaseDens[er_up][esr_up][ei_up][0][ip] * dCompFracGrad_dC[ke][jc];
              }
            }

            // add contributions of the derivatives of upwind coefficient wrt pressure/component fractions
            dDispersionFlux_dP[k_up][ic] += m_dPhaseDens[er_up][esr_up][ei_up][0][ip][Deriv::dP] * compFracGrad;

            applyChainRule( numComp,
                            m_dCompFrac_dCompDens[er_up][esr_up][ei_up],
                            m_dPhaseDens[er_up][esr_up][ei_up][0][ip],
                            dDens_dC,
                            Deriv::dC );
            for( integer jc = 0; jc < numComp; ++jc )
            {
              dDispersionFlux_dC[k_up][ic][jc] += dDens_dC[jc] * compFracGrad;
            }

            // call the lambda in the phase loop to allow the reuse of the phase fluxes and their derivatives
            // possible use: assemble the derivatives wrt temperature, and the flux term of the energy equation for this phase
            dispersionFluxKernelOp( ip, ic, k, seri, sesri, sei, connectionIndex,
                                    k_up, seri[k_up], sesri[k_up], sei[k_up],
                                    compFracGrad );

          } // loop over components
        } // loop over phases

        // add dispersion flux to local flux and local flux jacobian
        addToLocalFluxAndJacobian( k,
                                   stack,
                                   dispersionFlux,
                                   dDispersionFlux_dP,
                                   dDispersionFlux_dC );

        connectionIndex++;
      }   // loop over k[1]
    }   // loop over k[0]
  }

  /**
   * @brief Compute the component fraction gradient at this interface
   * @param[in] ip the phase index
   * @param[in] ic the component index
   * @param[in] seri the region indices
   * @param[in] sesri the subregion indices
   * @param[in] sei the element indices
   * @param[out] compFracGrad the component fraction gradient
   * @param[out] dCompFracGrad_dP the derivatives of the component fraction gradient wrt pressure
   * @param[out] dCompFracGrad_dC the derivatives of the component fraction gradient wrt component densities
   */
  GEOS_HOST_DEVICE
  inline
  void computeCompFractionGradient( integer const ip,
                                    integer const ic,
                                    localIndex const (&seri)[numFluxSupportPoints],
                                    localIndex const (&sesri)[numFluxSupportPoints],
                                    localIndex const (&sei)[numFluxSupportPoints],
                                    real64 const (&trans)[numFluxSupportPoints],
                                    real64 & compFracGrad,
                                    real64 (& dCompFracGrad_dP)[numFluxSupportPoints],
                                    real64 (& dCompFracGrad_dC)[numFluxSupportPoints][numComp] ) const
  {
    using Deriv = constitutive::multifluid::DerivativeOffset;

    real64 dCompFrac_dC[numComp]{};

    for( integer i = 0; i < numFluxSupportPoints; i++ )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      compFracGrad += trans[i] * m_phaseCompFrac[er][esr][ei][0][ip][ic];
      dCompFracGrad_dP[i] += trans[i] * m_dPhaseCompFrac[er][esr][ei][0][ip][ic][Deriv::dP];

      applyChainRule( numComp,
                      m_dCompFrac_dCompDens[er][esr][ei],
                      m_dPhaseCompFrac[er][esr][ei][0][ip][ic],
                      dCompFrac_dC,
                      Deriv::dC );
      for( integer jc = 0; jc < numComp; ++jc )
      {
        dCompFracGrad_dC[i][jc] += trans[i] * dCompFrac_dC[jc];
      }
    }
  }

  /**
   * @brief Add the local diffusion/dispersion flux contributions to the residual and Jacobian
   * @param[in] k the cell indices
   * @param[in] stack the stack variables
   * @param[in] flux the diffusion/dispersion flux
   * @param[in] dFlux_dP the derivative of the diffusion/dispersion flux wrt pressure
   * @param[in] dFlux_dC the derivative of the diffusion/dispersion flux wrt compositions
   */
  GEOS_HOST_DEVICE
  inline
  void addToLocalFluxAndJacobian( localIndex const (&k)[numFluxSupportPoints],
                                  StackVariables & stack,
                                  real64 const (&flux)[numComp],
                                  real64 const (&dFlux_dP)[numFluxSupportPoints][numComp],
                                  real64 const (&dFlux_dC)[numFluxSupportPoints][numComp][numComp] ) const
  {
    // loop over components
    for( integer ic = 0; ic < numComp; ++ic )
    {
      // finally, increment local flux and local Jacobian
      integer const eqIndex0 = k[0] * numEqn + ic;
      integer const eqIndex1 = k[1] * numEqn + ic;

      stack.localFlux[eqIndex0] += m_dt * flux[ic];
      stack.localFlux[eqIndex1] -= m_dt * flux[ic];

      for( integer ke = 0; ke < numFluxSupportPoints; ++ke )
      {
        localIndex const localDofIndexPres = k[ke] * numDof;
        stack.localFluxJacobian[eqIndex0][localDofIndexPres] += m_dt * dFlux_dP[ke][ic];
        stack.localFluxJacobian[eqIndex1][localDofIndexPres] -= m_dt * dFlux_dP[ke][ic];

        for( integer jc = 0; jc < numComp; ++jc )
        {
          localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
          stack.localFluxJacobian[eqIndex0][localDofIndexComp] += m_dt * dFlux_dC[ke][ic][jc];
          stack.localFluxJacobian[eqIndex1][localDofIndexComp] -= m_dt * dFlux_dC[ke][ic][jc];
        }
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
  inline
  void complete( localIndex const iconn,
                 StackVariables & stack,
                 FUNC && assemblyKernelOp = NoOpFunc{} ) const
  {
    using namespace compositionalMultiphaseUtilities;

    if( m_kernelFlags.isSet( KernelFlags::TotalMassEquation ) )
    {
      // Apply equation/variable change transformation(s)
      stackArray1d< real64, maxStencilSize * numDof > work( stack.stencilSize * numDof );
      shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( numComp, numEqn, numDof*stack.stencilSize, stack.numConnectedElems,
                                                               stack.localFluxJacobian, work );
      shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( numComp, numEqn, stack.numConnectedElems,
                                                                 stack.localFlux );
    }

    // add contribution to residual and jacobian into:
    // - the component mass balance equations (i = 0 to i = numComp-1)
    // note that numDof includes derivatives wrt temperature if this class is derived in ThermalKernels
    for( integer i = 0; i < stack.numConnectedElems; ++i )
    {
      if( m_ghostRank[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = m_dofNumber[m_seri( iconn, i )][m_sesri( iconn, i )][m_sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - m_rankOffset );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GT( m_localMatrix.numRows(), localRow + numComp );

        for( integer ic = 0; ic < numComp; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &m_localRhs[localRow + ic], stack.localFlux[i * numEqn + ic] );
          m_localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >
            ( localRow + ic,
            stack.dofColIndices.data(),
            stack.localFluxJacobian[i * numEqn + ic].dataIfContiguous(),
            stack.stencilSize * numDof );
        }

        // call the lambda to assemble additional terms, such as thermal terms
        assemblyKernelOp( i, localRow );
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
          integer const hasDiffusion,
          integer const hasDispersion,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numConnections, [=] GEOS_HOST_DEVICE ( localIndex const iconn )
    {
      typename KERNEL_TYPE::StackVariables stack( kernelComponent.stencilSize( iconn ),
                                                  kernelComponent.numPointsInFlux( iconn ) );

      kernelComponent.setup( iconn, stack );
      if( hasDiffusion )
      {
        kernelComponent.computeDiffusionFlux( iconn, stack );
      }
      if( hasDispersion )
      {
        kernelComponent.computeDispersionFlux( iconn, stack );
      }
      kernelComponent.complete( iconn, stack );
    } );
  }

protected:

  /// Views on phase volume fraction
  ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const m_phaseVolFrac;

  /// Views on phase densities
  ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const m_phaseDens;
  ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const m_dPhaseDens;

  /// Views on diffusivity
  ElementViewConst< arrayView3d< real64 const > > const m_diffusivity;
  ElementViewConst< arrayView3d< real64 const > > const m_dDiffusivity_dTemp;
  ElementViewConst< arrayView3d< real64 const > > const m_phaseDiffusivityMultiplier;

  /// Views on dispersivity
  ElementViewConst< arrayView3d< real64 const > > const m_dispersivity;

  /// View on the reference porosity
  ElementViewConst< arrayView1d< real64 const > > const m_referencePorosity;

  // Stencil information

  /// Reference to the stencil wrapper
  STENCILWRAPPER const m_stencilWrapper;

  /// Connection to element maps
  typename STENCILWRAPPER::IndexContainerViewConstType const m_seri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sesri;
  typename STENCILWRAPPER::IndexContainerViewConstType const m_sei;

};

/**
 * @class DiffusionDispersionFluxComputeKernelFactory
 */
class DiffusionDispersionFluxComputeKernelFactory
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
   * @param[in] hasDiffusion flag specifying whether diffusion is used or not
   * @param[in] hasDispersion flag specifying whether dispersion is used or not
   * @param[in] solverName the name of the solver
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
                   integer const hasDiffusion,
                   integer const hasDispersion,
                   integer const useTotalMassEquation,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   STENCILWRAPPER const & stencilWrapper,
                   real64 const dt,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs )
  {
    isothermalCompositionalMultiphaseBaseKernels::internal::kernelLaunchSelectorCompSwitch( numComps, [&]( auto NC )
    {
      integer constexpr NUM_COMP = NC();
      integer constexpr NUM_DOF = NC() + 1;

      ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
        elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
      dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

      BitFlags< KernelFlags > kernelFlags;
      if( useTotalMassEquation )
        kernelFlags.set( KernelFlags::TotalMassEquation );

      using kernelType = DiffusionDispersionFluxComputeKernel< NUM_COMP, NUM_DOF, STENCILWRAPPER >;
      typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
      typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
      typename kernelType::DiffusionAccessors diffusionAccessors( elemManager, solverName );
      typename kernelType::DispersionAccessors dispersionAccessors( elemManager, solverName );
      typename kernelType::PorosityAccessors porosityAccessors( elemManager, solverName );

      kernelType kernel( numPhases, rankOffset, stencilWrapper,
                         dofNumberAccessor, compFlowAccessors, multiFluidAccessors,
                         diffusionAccessors, dispersionAccessors, porosityAccessors,
                         dt, localMatrix, localRhs, kernelFlags );
      kernelType::template launch< POLICY >( stencilWrapper.size(),
                                             hasDiffusion, hasDispersion,
                                             kernel );
    } );
  }
};

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_DIFFUSIONDISPERSIONFLUXCOMPUTEKERNEL_HPP
