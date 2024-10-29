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
 * @file MultiphasePoromechanics.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICS_HPP_

#include "codingUtilities/Utilities.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/multiphysics/PoromechanicsFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/PoromechanicsBase.hpp"

namespace geos
{

namespace poromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static multiphase poromechanics.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 *
 * ### MultiphasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static multiphase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class MultiphasePoromechanics :
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
  static constexpr int maxNumComponents = 3;
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
  using Base::m_pressure_n;
  using Base::m_pressure;
  using Base::m_flowDofNumber;
  using Base::m_meshData;
  using Base::m_dt;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param gravityVector The gravity vector.
   */
  MultiphasePoromechanics( NodeManager const & nodeManager,
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
                           localIndex const numComponents,
                           localIndex const numPhases,
                           integer const useSimpleAccumulation,
                           integer const useTotalMassEquation,
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
    real64 compMassIncrement[maxNumComponents]{};
    /// Derivative of mass accumulation wrt volumetric strain increment
    real64 dCompMassIncrement_dVolStrainIncrement[maxNumComponents]{};
    /// Derivative of mass accumulation wrt pressure
    real64 dCompMassIncrement_dPressure[maxNumComponents]{};
    /// Derivative of mass accumulation wrt comp density
    real64 dCompMassIncrement_dComponents[maxNumComponents][maxNumComponents]{};

    /// Pore volume constraint
    real64 poreVolConstraint{};
    /// Derivative of pore volume constraint wrt pressure
    real64 dPoreVolConstraint_dPressure{};
    /// Derivative of pore volume constraint wrt comp density
    real64 dPoreVolConstraint_dComponents[1][maxNumComponents]{};

    /// Derivative of body force wrt comp density
    real64 dBodyForce_dComponents[3][maxNumComponents]{};

    // Storage for mass residual and degrees of freedom

    /// Derivative of momemtum balance residual wrt comp density
    real64 dLocalResidualMomentum_dComponents[Base::StackVariables::numDispDofPerElem][maxNumComponents]{};

    /// Mass balance residual
    real64 localResidualMass[maxNumComponents]{};
    /// Derivative of mass balance residual wrt volumetric strain increment
    real64 dLocalResidualMass_dDisplacement[maxNumComponents][Base::StackVariables::numDispDofPerElem]{};
    /// Derivative of mass balance residual wrt pressure
    real64 dLocalResidualMass_dPressure[maxNumComponents][1]{};
    /// Derivative of mass balance residual wrt components
    real64 dLocalResidualMass_dComponents[maxNumComponents][maxNumComponents]{};

    /// Pore volume constraint residual (sum of saturations = 1)
    real64 localResidualPoreVolConstraint[1]{};
    /// Derivative of pore volume constraint residual wrt pressure
    real64 dLocalResidualPoreVolConstraint_dPressure[1][1]{};
    /// Derivative of pore volume constraint residual wrt components
    real64 dLocalResidualPoreVolConstraint_dComponents[1][maxNumComponents]{};

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localComponentDofIndices[maxNumComponents]{};

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::ImplicitKernelBase::setup
   *
   * For the MultiphasePoromechanics implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;


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
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[in] dSolidDensity_dPressure the derivative of solid density wrt pressure
   * @param[inout] stack the stack variables
   * @param[inout] bodyForceKernelOp the lambda function to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeBodyForce( localIndex const k,
                         localIndex const q,
                         real64 const & porosity,
                         real64 const & dPorosity_dVolStrain,
                         real64 const & dPorosity_dPressure,
                         real64 const & dPorosity_dTemperature,
                         real64 const & dSolidDensity_dPressure,
                         StackVariables & stack,
                         FUNC && bodyForceKernelOp = NoOpFunc{} ) const;

  /**
   * @brief Helper function to compute the fluid mass increment and its derivatives wrt primary variables
   * @tparam FUNC the type of the function that can be used to customize the kernel
   * @param[in] k the element index
   * @param[in] q the quadrature point index
   * @param[in] porosity the element porosity
   * @param[in] porosity_n the element porosity at the previous converged time step
   * @param[in] dPorosity_dVolStrain the derivative of porosity wrt volumetric strain increment
   * @param[in] dPorosity_dPressure the derivative of porosity wrt pressure
   * @param[in] dPorosity_dTemperature the derivative of porosity wrt temperature
   * @param[inout] stack the stack variables
   * @param[inout] fluidIncrementKernelOp the lambda function to customize the kernel
   */
  template< typename FUNC = NoOpFunc >
  GEOS_HOST_DEVICE
  void computeFluidIncrement( localIndex const k,
                              localIndex const q,
                              real64 const & porosity,
                              real64 const & porosity_n,
                              real64 const & dPorosity_dVolStrain,
                              real64 const & dPorosity_dPressure,
                              real64 const & dPorosity_dTemperature,
                              StackVariables & stack,
                              FUNC && fluidIncrementKernelOp = NoOpFunc{} ) const;

  /**
   * @brief Helper function to compute the pore-volume constraint and its derivatives wrt primary variables
   * @param[in] k the element index
   * @param[in] porosity_n the element porosity at the previous converged time step
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  void computePoreVolumeConstraint( localIndex const k,
                                    real64 const & porosity_n,
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
   * ### MultiphasePoromechanics Description
   * Copy of the KernelBase::keranelLaunch function
   * elements.
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );


protected:

  /// Views on phase densities and derivatives
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseDensity;
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseDensity_n;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseDensity;

  /// Views on phase component fractions and derivatives
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > m_fluidPhaseCompFrac;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_COMP > m_fluidPhaseCompFrac_n;
  arrayView5d< real64 const, constitutive::multifluid::USD_PHASE_COMP_DC > m_dFluidPhaseCompFrac;

  /// Views on phase mass densities
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > m_fluidPhaseMassDensity;
  arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > m_dFluidPhaseMassDensity;

  /// Views on phase volume fraction (saturation)
  arrayView2d< real64 const, compflow::USD_PHASE > const m_fluidPhaseVolFrac;
  arrayView2d< real64 const, compflow::USD_PHASE > const m_fluidPhaseVolFrac_n;
  arrayView3d< real64 const, compflow::USD_PHASE_DC > const m_dFluidPhaseVolFrac;

  /// Views on derivatives of global comp fraction wrt global comp density
  arrayView3d< real64 const, compflow::USD_COMP_DC > const m_dGlobalCompFraction_dGlobalCompDensity;

  // Views on component densities
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens;
  arrayView2d< real64 const, compflow::USD_COMP > m_compDens_n;

  /// Number of components
  localIndex const m_numComponents;

  /// Number of phases
  localIndex const m_numPhases;

  /// Use simple accumulation term form
  integer const m_useSimpleAccumulation;

  /// Use total mass equation flag
  integer const m_useTotalMassEquation;

  integer const m_performStressInitialization;
};

using MultiphasePoromechanicsKernelFactory =
  finiteElement::KernelFactory< MultiphasePoromechanics,
                                arrayView1d< globalIndex const > const,
                                globalIndex const,
                                CRSMatrixView< real64, globalIndex const > const,
                                arrayView1d< real64 > const,
                                real64 const,
                                real64 const (&)[3],
                                string const,
                                localIndex const,
                                localIndex const,
                                integer const,
                                integer const,
                                integer const,
                                string const >;

/**
 * @class MultiphaseBulkDensityKernel
 * @brief Kernel to update the bulk density before a mechanics solve in sequential schemes
 */
class MultiphaseBulkDensityKernel
{
public:

  /**
   * @brief Constructor
   * @param[in] numPhases the number of fluid phases
   * @param[in] fluid the fluid model
   * @param[in] solid the porous solid model
   * @param[in] subRegion the element subregion
   */
  MultiphaseBulkDensityKernel( integer const numPhases,
                               constitutive::MultiFluidBase const & fluid,
                               constitutive::CoupledSolidBase const & solid,
                               ElementSubRegionBase & subRegion )
    : m_numPhases( numPhases ),
    m_bulkDensity( subRegion.getField< fields::poromechanics::bulkDensity >() ),
    m_fluidPhaseVolFrac( subRegion.getField< fields::flow::phaseVolumeFraction >() ),
    m_rockDensity( solid.getDensity() ),
    m_fluidPhaseDensity( fluid.phaseMassDensity() ),
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
    m_bulkDensity[ei][q] = 0.0;
    for( integer ip = 0; ip < m_numPhases; ++ip )
    {
      m_bulkDensity[ei][q] += m_fluidPhaseVolFrac[ei][ip] * m_fluidPhaseDensity[ei][q][ip];
    }
    m_bulkDensity[ei][q] *= m_porosity[ei][q];
    m_bulkDensity[ei][q] += ( 1 - m_porosity[ei][q] ) * m_rockDensity[ei][q];
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

  // number of fluid phases
  integer const m_numPhases;

  // the bulk density
  arrayView2d< real64 > const m_bulkDensity;

  // the fluid phase saturation
  arrayView2d< real64 const, compflow::USD_PHASE > const m_fluidPhaseVolFrac;

  // the rock density
  arrayView2d< real64 const > const m_rockDensity;

  // the fluid density
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const m_fluidPhaseDensity;

  // the porosity
  arrayView2d< real64 const > const m_porosity;

};

/**
 * @class MultiphaseBulkDensityKernelFactory
 */
class MultiphaseBulkDensityKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] numPhases number of phases
   * @param[in] fluid the fluid model
   * @param[in] solid the porous solid model
   * @param[in] subRegion the element subregion
   */
  template< typename POLICY >
  static void
  createAndLaunch( integer const numPhases,
                   constitutive::MultiFluidBase const & fluid,
                   constitutive::CoupledSolidBase const & solid,
                   ElementSubRegionBase & subRegion )
  {
    MultiphaseBulkDensityKernel kernel( numPhases, fluid, solid, subRegion );
    MultiphaseBulkDensityKernel::launch< POLICY >( subRegion.size(),
                                                   subRegion.getField< fields::poromechanics::bulkDensity >().size( 1 ),
                                                   kernel );
  }
};


} // namespace poromechanicsKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_MULTIPHASEPOROMECHANICS_HPP_
