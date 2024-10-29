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
 * @file ExplicitSmallStrain.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITSMALLSTRAIN_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITSMALLSTRAIN_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"


namespace geos
{

/// Namespace to contain the solid mechanics kernels.
namespace solidMechanicsLagrangianFEMKernels
{

/// If UPDATE_STRESS is undef, uses total displacement and stress is not
/// updated at all.
/// If UPDATE_STRESS 1, uses total displacement to and adds material stress
/// state to integral for nodalforces.
/// If UPDATE_STRESS 2 then velocity*dt is used to update material stress state
#define UPDATE_STRESS 2

/**
 * @brief Implements kernels for solving the equations of motion using the
 *   explicit Newmark method under the small strain assumption.
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### Explicit Small Strain Description
 * Implements the KernelBase interface functions required for explicit time
 * integration of the equations of motion using the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the interface for KernelBase is used, but
 * ExplicitSmallStrain only conforms to the interface set by KernelBase, and
 * does not inherit from KernelBase.
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3`.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitSmallStrain : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              FE_TYPE,
                                                              3,
                                                              3 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          3,
                                          3 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   * @param elementListName The name of the entry that holds the list of
   *   elements to be processed during this kernel launch.
   */
  ExplicitSmallStrain( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       real64 const dt,
                       string const elementListName );

  //*****************************************************************************
  /**
   * @copydoc geos::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitSmallStrain Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      fLocal{ { 0.0} },
      varLocal{ {0.0} },
      xLocal()
    {}

    /// C-array stack storage for the element local force
    real64 fLocal[ numNodesPerElem ][ numDofPerTrialSupportPoint ];

    /// C-array stack storage for element local primary variable values.
    real64 varLocal[ numNodesPerElem ][ numDofPerTestSupportPoint ];

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif
  };
  //***************************************************************************


  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitSmallStrain Description
   * Calculates the shape function derivatives, and the strain tensor. Then
   * calls the constitutive update, and also performs the integration of
   * the stress divergence, rather than using the dedicated component function
   * to allow for some variable reuse.
   */
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc geos::finiteElement::KernelBase::complete
   *
   * ### ExplicitSmallStrain Description
   * Performs the distribution of the nodal force out to the rank local arrays.
   */
  GEOS_HOST_DEVICE
  real64 complete( localIndex const k,
                   StackVariables const & stack ) const;

  /**
   * @copydoc geos::finiteElement::KernelBase::kernelLaunch
   *
   * ### ExplicitSmallStrain Description
   * Copy of the KernelBase::kernelLaunch function without the exclusion of ghost
   * elements.
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );


protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_u;

  /// The array containing the nodal velocity array.
  arrayView2d< real64 const, nodes::VELOCITY_USD > const m_vel;

  /// The array containing the nodal acceleration array, which is used to store
  /// the force.
  arrayView2d< real64, nodes::ACCELERATION_USD > const m_acc;

  /// The time increment for this time integration step.
  real64 const m_dt; ///TODO: Consider moving to finite element kernel base?

  /// The list of elements to process for the kernel launch.
  SortedArrayView< localIndex const > const m_elementList;


};



/// The factory used to construct a ExplicitSmallStrain kernel.
using ExplicitSmallStrainFactory = finiteElement::KernelFactory< ExplicitSmallStrain,
                                                                 real64,
                                                                 string const >;



} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITSMALLSTRAIN_HPP_
