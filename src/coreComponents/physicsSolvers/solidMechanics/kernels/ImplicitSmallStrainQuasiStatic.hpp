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
 * @file ImplicitSmallStrainQuasiStatic.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{

namespace solidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### QuasiStatic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static equilibrium equations using one of the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ImplicitSmallStrainQuasiStatic :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
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
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_meshData;
  using Base::m_dt;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  ImplicitSmallStrainQuasiStatic( NodeManager const & nodeManager,
                                  EdgeManager const & edgeManager,
                                  FaceManager const & faceManager,
                                  localIndex const targetRegionIndex,
                                  SUBREGION_TYPE const & elementSubRegion,
                                  FE_TYPE const & finiteElementSpace,
                                  CONSTITUTIVE_TYPE & inputConstitutiveType,
                                  arrayView1d< globalIndex const > const inputDofNumber,
                                  globalIndex const rankOffset,
                                  CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                  arrayView1d< real64 > const inputRhs,
                                  real64 const inputDt,
                                  real64 const (&inputGravityVector)[3] );

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
      Base::StackVariables(),
                                       xLocal(),
                                       u_local(),
                                       uhat_local(),
                                       constitutiveStiffness()
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[ 6 ][ 6 ];
  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::ImplicitKernelBase::setup
   *
   * For the QuasiStatic implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;


  /**
   * @brief Internal struct to provide no-op defaults used in the inclusion
   *   of lambda functions into kernel component functions.
   * @struct NoOpFunctors
   */
  struct NoOpFunctors
  {
    /**
     * @brief operator() no-op used for adding an additional dynamics term
     *   inside the jacobian assembly loop.
     * @param a Node index for the row.
     * @param b Node index for the col.
     */
    GEOS_HOST_DEVICE inline constexpr
    void operator() ( localIndex const a, localIndex const b )
    {
      GEOS_UNUSED_VAR( a );
      GEOS_UNUSED_VAR( b );
    }

    /**
     * @brief operator() no-op used for modifying the stress tensor prior to
     *   integrating the divergence to produce nodal forces.
     * @param stress The stress array.
     */
    GEOS_HOST_DEVICE inline constexpr
    void operator() ( real64 (& stress)[6] )
    {
      GEOS_UNUSED_VAR( stress );
    }
  };

  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   * @tparam STRESS_MODIFIER Type of optional functor to allow for the
   * modification of stress prior to integration.
   * @param stressModifier An optional functor to allow for the modification
   *  of stress prior to integration.
   * For solid mechanics kernels, the strain increment is calculated, and the
   * constitutive update is called. In addition, the constitutive stiffness
   * stack variable is filled by the constitutive model.
   */
  template< typename STRESS_MODIFIER = NoOpFunctors >
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack,
                              STRESS_MODIFIER && stressModifier = NoOpFunctors{} ) const;

  /**
   * @copydoc geos::finiteElement::ImplicitKernelBase::complete
   */
  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const;

  /**
   * @copydoc geos::finiteElement::KernelBase::kernelLaunch
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent );

protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];

  /// The rank global density
  arrayView2d< real64 const > const m_density;

  /**
   * @brief Get a parameter representative of the stiffness, used as physical scaling for the
   * stabilization matrix.
   * @param[in] k Element index.
   * @return A parameter representative of the stiffness matrix dstress/dstrain
   */
  GEOS_HOST_DEVICE
  inline
  real64 computeStabilizationScaling( localIndex const k ) const
  {
    // TODO: generalize this to other constitutive models (currently we assume linear elasticity).
    return 2.0 * m_constitutiveUpdate.getShearModulus( k );
  }
};

/// The factory used to construct a QuasiStatic kernel.
using QuasiStaticFactory = finiteElement::KernelFactory< ImplicitSmallStrainQuasiStatic,
                                                         arrayView1d< globalIndex const > const,
                                                         globalIndex,
                                                         CRSMatrixView< real64, globalIndex const > const,
                                                         arrayView1d< real64 > const,
                                                         real64 const,
                                                         real64 const (&)[3] >;

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINQUASISTATIC_HPP_
