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
 * @file ExplicitFiniteStrain.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_HPP_

#include "ExplicitSmallStrain.hpp"
#include "finiteElement/Kinematics.h"

namespace geos
{

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
 *   explicit Newmark method under the finite strain assumption.
 * @copydoc ExplicitSmallStrain
 *
 * ### Explicit Small Strain Description
 * Finite strain implementation.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitFiniteStrain : public ExplicitSmallStrain< SUBREGION_TYPE,
                                                         CONSTITUTIVE_TYPE,
                                                         FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = ExplicitSmallStrain< SUBREGION_TYPE,
                                    CONSTITUTIVE_TYPE,
                                    FE_TYPE >;

  using Base::numNodesPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

  using Base::m_dt;
  using Base::m_u;
  using Base::m_vel;
  using Base::m_acc;
#if !defined(CALCFEMSHAPE)
  using Base::m_X;
#endif

  /**
   * @copydoc ExplicitSmallStrain
   */
  ExplicitFiniteStrain( NodeManager & nodeManager,
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
   * @copydoc ExplicitSmallStrain::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
    using Base::StackVariables::fLocal;
    using Base::StackVariables::varLocal;
    using Base::StackVariables::xLocal;


    /**
     * @brief constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            uLocal{ {0.0} }
    {}

    /// Local stack storage for nodal displacements.
    real64 uLocal[ numNodesPerElem ][ numDofPerTrialSupportPoint ];
  };
  //*****************************************************************************


  /**
   * @copydoc ExplicitSmallStrain::setup
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  /**
   * @copydoc ExplicitSmallStrain::quadraturePointKernel
   */
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;


  /**
   * @copydoc ExplicitSmallStrain::kernelLaunch
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );



};
/// The factory used to construct a ExplicitFiniteStrain kernel.
using ExplicitFiniteStrainFactory = finiteElement::KernelFactory< ExplicitFiniteStrain,
                                                                  real64,
                                                                  string const >;

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_EXPLICITFINITESTRAIN_HPP_
