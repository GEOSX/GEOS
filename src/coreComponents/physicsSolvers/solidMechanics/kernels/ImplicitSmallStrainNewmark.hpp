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
 * @file ImplicitSmallStrainNewmark.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINNEWMARK_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINNEWMARK_HPP_

#include "ImplicitSmallStrainQuasiStatic.hpp"


namespace geos
{

namespace solidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving the equations of motion using an
 *   implicit Newmark's method..
 * @copydoc ImplicitSmallStrainQuasiStatic
 *
 * ### Implicit Newmark Description
 * Implements the KernelBase interface functions required for solving the
 * equations of motion using with an Implicit Newmark's Method with one of the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ImplicitSmallStrainNewmark : public ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE,
                                                                          CONSTITUTIVE_TYPE,
                                                                          FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE,
                                               CONSTITUTIVE_TYPE,
                                               FE_TYPE >;

  using Base::numNodesPerElem;
  using Base::maxNumTestSupportPointsPerElem;
  using Base::maxNumTrialSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;

  using Base::m_dofNumber;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_disp;
  using Base::m_uhat;
  using Base::m_density;
  using Base::m_finiteElementSpace;
  using Base::m_dt;

#if !defined(CALCFEMSHAPE)
  using Base::m_X;
#endif
  using Base::m_gravityVector;

  /**
   * @brief Constructor
   * @copydoc QuasiStatic
   * @param inputNewmarkGamma The Gamma parameter of the Newmark method.
   * @param inputNewmarkBeta The Beta parameter for the Newmark method.
   * @param inputMassDamping The mass damping coefficient.
   * @param inputStiffnessDamping The stiffness damping coefficient.
   * @param inputDt The timestep for the physics update.
   */
  ImplicitSmallStrainNewmark( NodeManager const & nodeManager,
                              EdgeManager const & edgeManager,
                              FaceManager const & faceManager,
                              localIndex const targetRegionIndex,
                              SUBREGION_TYPE const & elementSubRegion,
                              FE_TYPE const & finiteElementSpace,
                              CONSTITUTIVE_TYPE & inputConstitutiveType,
                              arrayView1d< globalIndex const > const & inputDofNumber,
                              globalIndex const rankOffset,
                              CRSMatrixView< real64, globalIndex const > const inputMatrix,
                              arrayView1d< real64 > const inputRhs,
                              real64 const inputDt,
                              real64 const (&inputGravityVector)[3],
                              real64 const inputNewmarkGamma,
                              real64 const inputNewmarkBeta,
                              real64 const inputMassDamping,
                              real64 const inputStiffnessDamping );

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc QuasiStatic::StackVariables
   *
   * Adds a stack array for the vtilde, uhattilde, and the
   * Inertial mass damping.
   */
  struct StackVariables : public Base::StackVariables
  {
public:
    using Base::StackVariables::maxNumRows;
    using Base::StackVariables::maxNumCols;

    /// Constructor.
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dRdU_InertiaMassDamping{ {0.0} },
      vtilde_local(),
      uhattilde_local()
    {}

    /// Stack storage for the Inertial damping contributions to the Jacobian
    real64 dRdU_InertiaMassDamping[ maxNumRows ][ maxNumCols ];

    /// Stack storage for the velocity predictor.
    real64 vtilde_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the incremental displacement predictor.
    real64 uhattilde_local[numNodesPerElem][numDofPerTrialSupportPoint];
  };
  //***************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc QuasiStatic::setup
   *
   * For the ImplicitNewmark implementation, global values from the velocity
   * predictor, and the incremental displacement predictor are placed into
   * element local stack storage.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const;

  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * The ImplcitNewmark kernel adds the calculation of the inertia damping,
   * jacobian and residual contributions.
   */
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const;

  /**
   * @copydoc QuasiStatic::complete
   *
   * The ImplicitNewmark implementation adds residual and jacobian
   * contributions from  stiffness based damping.
   */
  GEOS_HOST_DEVICE
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
  /// The rank-global velocity predictor
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_vtilde;

  /// The rank-global incremental displacement predictor
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhattilde;

  /// The Gamma parameter for Newmark's method.
  real64 const m_newmarkGamma;

  /// The Beta parameter for Newmark's method.
  real64 const m_newmarkBeta;

  /// The mass damping coefficient.
  real64 const m_massDamping;

  /// The stiffness damping coefficient.
  real64 const m_stiffnessDamping;


};

/// The factory used to construct a ImplicitNewmark kernel.
using ImplicitNewmarkFactory = finiteElement::KernelFactory< ImplicitSmallStrainNewmark,
                                                             arrayView1d< globalIndex const > const &,
                                                             globalIndex,
                                                             CRSMatrixView< real64, globalIndex const > const,
                                                             arrayView1d< real64 > const,
                                                             real64 const,
                                                             real64 const (&)[3],
                                                             real64,
                                                             real64,
                                                             real64,
                                                             real64 >;

} // namespace solidMechanicsLagrangianFEMKernels

} // namespace geos

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif //GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_IMPLCITSMALLSTRAINNEWMARK_HPP_
