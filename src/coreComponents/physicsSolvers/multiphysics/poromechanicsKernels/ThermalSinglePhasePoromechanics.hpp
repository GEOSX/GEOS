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
 * @file ThermalSinglePhasePoromechanics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICS_HPP_

#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"

namespace geos
{

namespace thermalPoromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static thermal single-phase poromechanics.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ThermalSinglePhasePoromechanics :
  public poromechanicsKernels::SinglePhasePoromechanics< SUBREGION_TYPE,
                                                         CONSTITUTIVE_TYPE,
                                                         FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = poromechanicsKernels::SinglePhasePoromechanics< SUBREGION_TYPE,
                                                               CONSTITUTIVE_TYPE,
                                                               FE_TYPE >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_gravityAcceleration;
  using Base::m_gravityVector;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_pressure;
  using Base::m_pressure_n;
  using Base::m_fluidDensity;
  using Base::m_fluidDensity_n;
  using Base::m_dFluidDensity_dPressure;
  using Base::m_solidDensity;
  using Base::m_flowDofNumber;
  using Base::m_dt;
  using Base::m_performStressInitialization;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param gravityVector The gravity vector.
   */
  ThermalSinglePhasePoromechanics( NodeManager const & nodeManager,
                                   EdgeManager const & edgeManager,
                                   FaceManager const & faceManager,
                                   localIndex const targetRegionIndex,
                                   SUBREGION_TYPE const & elementSubRegion,
                                   FE_TYPE const & finiteElementSpace,
                                   CONSTITUTIVE_TYPE & inputConstitutiveType,
                                   arrayView1d< globalIndex const > const inputDispDofNumber,
                                   globalIndex const rankOffset,
                                   CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                   arrayView1d< real64 > const inputRhs,
                                   real64 const inputDt,
                                   real64 const (&gravityVector)[3],
                                   string const inputFlowDofKey,
                                   integer const performStressInitialization,
                                   string const fluidModelKey );

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geos::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    static constexpr int numDispDofPerElem =  Base::StackVariables::maxNumRows;

    /// Constructor.
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dLocalResidualMomentum_dTemperature{ {0.0} },
      dLocalResidualMass_dTemperature{ {0.0} },
      localTemperatureDofIndex{}
    {}

    // Storage for helper variables used in the quadrature point kernel

    /// Derivative of body force wrt temperature
    real64 dBodyForce_dTemperature[3]{};
    /// Derivative of mass accumulation wrt temperature
    real64 dFluidMassIncrement_dTemperature{};

    /// Energy accumulation
    real64 energyIncrement{};
    /// Derivative of energy accumulation wrt volumetric strain increment
    real64 dEnergyIncrement_dVolStrainIncrement{};
    /// Derivative of energy accumulation wrt pressure
    real64 dEnergyIncrement_dPressure{};
    /// Derivative of energy accumulation wrt temperature
    real64 dEnergyIncrement_dTemperature{};

    // Storage for mass/energy residual and degrees of freedom

    /// Derivative of linear momentum balance residual wrt temperature
    real64 dLocalResidualMomentum_dTemperature[numDispDofPerElem][1]{};
    /// Derivative of mass balance residual wrt pressure
    real64 dLocalResidualMass_dTemperature[1][1]{};

    /// Energy balance residual
    real64 localResidualEnergy[1]{};
    /// Derivative of energy balance residual wrt displacement
    real64 dLocalResidualEnergy_dDisplacement[1][numDispDofPerElem]{};
    /// Derivative of energy balance residual wrt pressure
    real64 dLocalResidualEnergy_dPressure[1][1]{};
    /// Derivative of energy balance residual wrt temperature
    real64 dLocalResidualEnergy_dTemperature[1][1]{};

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localTemperatureDofIndex{};

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::ImplicitKernelBase::setup
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  /**
   * @brief Helper function to compute 1) the total stress, 2) the body force term, and 3) the fluidMass/EnergyIncrement
   * using quantities returned by the PorousSolid constitutive model.
   * This function also computes the derivatives of these three quantities wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          StackVariables & stack ) const;

  /**
   * @brief Helper function to compute the body force term (\rho g) and its derivatives wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[in] dSolidDensity_dPressure the derivative of solid density wrt pressure
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeBodyForce( localIndex const k,
                         localIndex const q,
                         real64 const & porosity,
                         real64 const & dPorosity_dVolStrain,
                         real64 const & dPorosity_dPressure,
                         real64 const & dPorosity_dTemperature,
                         real64 const & dSolidDensity_dPressure,
                         StackVariables & stack ) const;

  /**
   * @brief Helper function to compute the fluid mass/energy increment and its derivatives wrt primary variables
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] porosity_n the element porosity at the previous converged time step
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computeFluidIncrement( localIndex const k,
                              localIndex const q,
                              real64 const & porosity,
                              real64 const & porosity_n,
                              real64 const & dPorosity_dVolStrain,
                              real64 const & dPorosity_dPressure,
                              real64 const & dPorosity_dTemperature,
                              StackVariables & stack ) const;

  /**
   * @brief Assemble the local linear momentum balance residual and derivatives using total stress and body force terms
   * @param[in] N displacement finite element basis functions
   * @param[in] dNdX basis function derivatives
   * @param[in] detJxW determinant of the Jacobian transformation matrix times the quadrature weight
   * @param[inout] stack the stack variables
   * @detail This function assembles the discretized form of the following equation
   *   divergence( totalStress ) + bodyForce = 0
   * with the following dependencies on the strainIncrement tensor, pressure, and temperature
   *   totalStress = totalStress( strainIncrement, pressure, temperature )
   *   bodyForce   = bodyForce( strainIncrement, pressure, temperature )
   */
  GEOS_HOST_DEVICE
  void assembleMomentumBalanceTerms( real64 const ( &N )[numNodesPerElem],
                                     real64 const ( &dNdX )[numNodesPerElem][3],
                                     real64 const & detJxW,
                                     StackVariables & stack ) const;

  /**
   * @brief Assemble the local mass/energy balance residual and derivatives using fluid mass/energy increment
   * @param[in] dNdX basis function derivatives
   * @param[in] detJxW determinant of the Jacobian transformation matrix times the quadrature weight
   * @param[inout] stack the stack variables
   * @detail This function assembles the discretized form of the temporal derivative in the following equation
   *   dFluidMass_dTime + divergence( fluidMassFlux ) = source  (fluid phase mass balance)
   * with the following dependencies on the strainIncrement tensor, pressure, and temperature
   *   fluidMass = fluidMass( strainIncrement, pressure, temperature )
   *   fluidMassFlux = fluidMassFlux( pressure, temperature )
   */
  GEOS_HOST_DEVICE
  void assembleElementBasedFlowTerms( real64 const ( &dNdX )[numNodesPerElem][3],
                                      real64 const & detJxW,
                                      StackVariables & stack ) const;

  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geos::finiteElement::ImplicitKernelBase::complete
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

  /**
   * @copydoc geos::finiteElement::KernelBase::kernelLaunch
   *
   * ### SinglePhasePoromechancis Description
   * Copy of the KernelBase::kernelLaunch function
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );

protected:

  /// Views on fluid density derivative wrt temperature
  arrayView2d< real64 const > const m_dFluidDensity_dTemperature;

  /// Views on fluid internal energy
  arrayView2d< real64 const > const m_fluidInternalEnergy_n;
  arrayView2d< real64 const > const m_fluidInternalEnergy;
  arrayView2d< real64 const > const m_dFluidInternalEnergy_dPressure;
  arrayView2d< real64 const > const m_dFluidInternalEnergy_dTemperature;

  /// Views on rock internal energy
  arrayView2d< real64 const > const m_rockInternalEnergy_n;
  arrayView2d< real64 const > const m_rockInternalEnergy;
  arrayView2d< real64 const > const m_dRockInternalEnergy_dTemperature;

  /// Views on temperature
  arrayView1d< real64 const > const m_temperature_n;
  arrayView1d< real64 const > const m_temperature;

};

using ThermalSinglePhasePoromechanicsKernelFactory =
  finiteElement::KernelFactory< ThermalSinglePhasePoromechanics,
                                arrayView1d< globalIndex const > const,
                                globalIndex const,
                                CRSMatrixView< real64, globalIndex const > const,
                                arrayView1d< real64 > const,
                                real64 const,
                                real64 const (&)[3],
                                string const,
                                integer const,
                                string const >;

} // namespace thermalPoromechanicsKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICS_HPP_
