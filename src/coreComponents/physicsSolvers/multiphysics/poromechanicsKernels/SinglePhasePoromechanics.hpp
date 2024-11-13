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
 * @file SinglePhasePoromechanics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_HPP_

#include "physicsSolvers/multiphysics/PoromechanicsFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/PoromechanicsBase.hpp"

namespace geos
{

namespace poromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase poromechanics.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 *
 * ### SinglePhasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static single-phase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhasePoromechanics :
  public PoromechanicsBase< SUBREGION_TYPE,
                            CONSTITUTIVE_TYPE,
                            FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = PoromechanicsBase< SUBREGION_TYPE,
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
  using Base::m_solidDensity;
  using Base::m_pressure;
  using Base::m_pressure_n;
  using Base::m_meshData;
  using Base::m_dt;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param gravityVector The gravity vector.
   */
  SinglePhasePoromechanics( NodeManager const & nodeManager,
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

    /// Constructor.
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables()
    {}

    /// Mass accumulation
    real64 fluidMassIncrement{};
    /// Derivative of mass accumulation wrt volumetric strain increment
    real64 dFluidMassIncrement_dVolStrainIncrement{};
    /// Derivative of mass accumulation wrt pressure
    real64 dFluidMassIncrement_dPressure{};

    // Storage for mass residual and degrees of freedom

    /// Mass balance residual
    real64 localResidualMass[1]{};
    /// Derivative of mass balance residual wrt displacement
    real64 dLocalResidualMass_dDisplacement[1][Base::StackVariables::numDispDofPerElem]{};
    /// Derivative of mass balance residual wrt pressure
    real64 dLocalResidualMass_dPressure[1][1]{};

  };
  //*****************************************************************************

  /**
   * @brief Helper function to compute 1) the total stress, 2) the body force term, and 3) the fluidMassIncrement
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
   * @brief Helper function to compute the fluid mass increment and its derivatives wrt primary variables
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
   * with the following dependencies on the strainIncrement tensor and pressure
   *   totalStress = totalStress( strainIncrement, pressure)
   *   bodyForce   = bodyForce( strainIncrement, pressure)
   */
  GEOS_HOST_DEVICE
  void assembleMomentumBalanceTerms( real64 const ( &N )[numNodesPerElem],
                                     real64 const ( &dNdX )[numNodesPerElem][3],
                                     real64 const & detJxW,
                                     StackVariables & stack ) const;

  /**
   * @brief Assemble the local mass balance residual and derivatives using fluid mass/energy increment
   * @param[in] dNdX basis function derivatives
   * @param[in] detJxW determinant of the Jacobian transformation matrix times the quadrature weight
   * @param[inout] stack the stack variables
   * @detail This function assembles the discretized form of the temporal derivative in the following equation
   *   dFluidMass_dTime + divergence( fluidMassFlux ) = source  (fluid phase mass balance)
   * with the following dependencies on the strainIncrement tensor and pressure
   *   fluidMass = fluidMass( strainIncrement, pressure)
   *   fluidMassFlux = fluidMassFlux( pressure)
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

  /// Fluid density
  arrayView2d< real64 const > const m_fluidDensity;
  /// Fluid density at the previous converged time step
  arrayView2d< real64 const > const m_fluidDensity_n;
  /// Derivative of fluid density wrt pressure
  arrayView2d< real64 const > const m_dFluidDensity_dPressure;

  integer const m_performStressInitialization;
};

using SinglePhasePoromechanicsKernelFactory =
  finiteElement::KernelFactory< SinglePhasePoromechanics,
                                arrayView1d< globalIndex const > const,
                                globalIndex const,
                                CRSMatrixView< real64, globalIndex const > const,
                                arrayView1d< real64 > const,
                                real64 const,
                                real64 const (&)[3],
                                string const,
                                integer const,
                                string const >;

/**
 * @class SinglePhaseBulkDensityKernel
 * @brief Kernel to update the bulk density before a mechanics solve in sequential schemes
 */
class SinglePhaseBulkDensityKernel
{
public:

  /**
   * @brief Constructor
   * @param[in] fluid the fluid model
   * @param[in] solid the porous solid model
   * @param[in] subRegion the element subregion
   */
  SinglePhaseBulkDensityKernel( constitutive::SingleFluidBase const & fluid,
                                constitutive::CoupledSolidBase const & solid,
                                ElementSubRegionBase & subRegion )
    : m_bulkDensity( subRegion.getField< fields::poromechanics::bulkDensity >() ),
    m_rockDensity( solid.getDensity() ),
    m_fluidDensity( fluid.density() ),
    m_porosity( solid.getPorosity() )
  {}

  /**
   * @brief Compute the bulk density in an element
   * @param[in] ei the element index
   * @param[in] q the quadrature point index
   */
  GEOS_HOST_DEVICE
  void compute( localIndex const ei,
                localIndex const q ) const
  {
    m_bulkDensity[ei][q] =
      ( 1 - m_porosity[ei][q] ) * m_rockDensity[ei][q] + m_porosity[ei][q] * m_fluidDensity[ei][q];
  }

  /**
   * @brief Performs the kernel launch
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] numElems the number of elements
   * @param[in] numQuad the number of quadrature points
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launch( localIndex const numElems,
          localIndex const numQuad,
          KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex q = 0; q < numQuad; ++q )
      {
        kernelComponent.compute( ei, q );
      }
    } );
  }

protected:

  // the bulk density
  arrayView2d< real64 > const m_bulkDensity;

  // the rock density
  arrayView2d< real64 const > const m_rockDensity;

  // the fluid density
  arrayView2d< real64 const > const m_fluidDensity;

  // the porosity
  arrayView2d< real64 const > const m_porosity;

};

/**
 * @class SinglePhaseBulkDensityKernelFactory
 */
class SinglePhaseBulkDensityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] fluid the fluid model
   * @param[in] solid the porous solid model
   * @param[in] subRegion the element subregion
   */
  template< typename POLICY >
  static void
  createAndLaunch( constitutive::SingleFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   ElementSubRegionBase & subRegion )
  {
    SinglePhaseBulkDensityKernel kernel( fluid, solid, subRegion );
    SinglePhaseBulkDensityKernel::launch< POLICY >( subRegion.size(),
                                                    subRegion.getField< fields::poromechanics::bulkDensity >().size( 1 ),
                                                    kernel );
  }
};


} // namespace poromechanicsKernels

} // namespace geos


#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICS_HPP_
