/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhasePoromechanicsDamageKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_HPP_

#include "physicsSolvers/multiphysics/PoromechanicsFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"

namespace geos
{

namespace poromechanicsDamageKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase poromechanics with phase-field damage.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhasePoromechanicsDamage :
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
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param gravityVector The gravity vector.
   */
  SinglePhasePoromechanicsDamage( NodeManager const & nodeManager,
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
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
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
      Base::StackVariables()
    {}

    /// Fracture flow term
    real64 fractureFlowTerm[3]{};
    /// Derivative of the fracture flow term wrt pressure
    real64 dFractureFlowTerm_dPressure[3]{};

  };
  //*****************************************************************************

  /**
   * @brief Helper function to compute 1) the total stress, 2) the body force term, and 3) the fluidMassIncrement
   * using quantities returned by the PorousDamageSolid constitutive model.
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

  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
   *
   * ### SinglePhasePoromechancisDamage Description
   * Copy of the KernelBase::kernelLaunch function
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );

protected:

  /// Views on cell-wise fluid pressure gradient
  arrayView2d< real64 const > const m_fluidPressureGradient;

};

using SinglePhasePoromechanicsDamageKernelFactory =
  finiteElement::KernelFactory< SinglePhasePoromechanicsDamage,
                                arrayView1d< globalIndex const > const,
                                globalIndex const,
                                CRSMatrixView< real64, globalIndex const > const,
                                arrayView1d< real64 > const,
                                real64 const,
                                real64 const (&)[3],
                                string const,
                                integer const,
                                string const >;

} // namespace poromechanicsDamageKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSDAMAGE_HPP_
