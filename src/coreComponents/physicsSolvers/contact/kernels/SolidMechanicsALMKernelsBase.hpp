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
 * @file SolidMechanicsALMKernelsBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMKERNELSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMKERNELSBASE_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"
#include "SolidMechanicsConformingContactKernelsHelper.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

/**
 * @brief A struct to check for constraint satisfaction
 */
struct ConstraintCheckKernel
{

  /**
   * @brief Launch the kernel function to check the constraint satisfaction
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONTACT_WRAPPER the type of contact wrapper doing the fracture traction updates
   * @param[in] size the size of the subregion
   * @param[in] traction the array containing the current traction
   * @param[in] dispJump the array containing the displacement jump
   * @param[in] deltaDispJump the array containing the delta displacement jump
   * @param[in] normalTractionTolerance Check tolerance (normal traction)
   * @param[in] normalDisplacementTolerance Check tolerance (compenetration)
   * @param[in] slidingTolerance Check tolerance (sliding)
   * @param[in] slidingCheckTolerance Check tolerance (if shear strass exceeds tauLim)
   * @param[in] area interface element area
   * @param[in] fractureState the array containing the fracture state
   * @param[out] condConv the array containing the convergence flag:
   *                      0: Constraint conditions satisfied
   *                      1: Open
   *                      2: Compenetration
   *                      3: Slip exceeds sliding tolerance
   *                      4: Shear stress exceeds tauLim
   */
  template< typename POLICY, typename CONTACT_WRAPPER >
  static void
  launch( localIndex const size,
          CONTACT_WRAPPER const & contactWrapper,
          arrayView1d< integer const > const & ghostRank,
          arrayView2d< real64 > const & traction,
          arrayView2d< real64 const > const & dispJump,
          arrayView2d< real64 const > const & deltaDispJump,
          arrayView1d< real64 const > const & normalTractionTolerance,
          arrayView1d< real64 const > const & normalDisplacementTolerance,
          arrayView1d< real64 const > const & slidingTolerance,
          real64 const slidingCheckTolerance,
          arrayView1d< integer const > const & fractureState,
          arrayView1d< integer > const & condConv )
  {

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      if( ghostRank[k] < 0 )
      {
        contactWrapper.constraintCheck( dispJump[k],
                                        deltaDispJump[k],
                                        traction[k],
                                        fractureState[k],
                                        normalTractionTolerance[k],
                                        normalDisplacementTolerance[k],
                                        slidingTolerance[k],
                                        slidingCheckTolerance,
                                        condConv[k] );
      }

    } );
  }
};

/**
 * @brief A struct to check for constraint satisfaction
 */
struct UpdateStateKernel
{

  /**
   * @brief Launch the kernel function to check the constraint satisfaction
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONTACT_WRAPPER the type of contact wrapper doing the fracture traction updates
   * @param[in] size the size of the subregion
   * @param[in] oldDispJump the array containing the old displacement jump (previous time step)
   * @param[in] dispJump the array containing the displacement jump
   * @param[in] penalty the array containing the penalty coefficients
   * @param[in] symmetric flag to compute symmetric penalty matrix
   * @param[in] normalTractionTolerance Check tolerance (normal traction)
   * @param[in] traction the array containing the current traction
   * @param[in] fractureState the array containing the fracture state
   */
  template< typename POLICY, typename CONTACT_WRAPPER >
  static void
  launch( localIndex const size,
          CONTACT_WRAPPER const & contactWrapper,
          arrayView2d< real64 const > const & oldDispJump,
          arrayView2d< real64 const > const & dispJump,
          arrayView2d< real64 > const & penalty,
          arrayView1d< real64 const > const & faceArea,
          bool const symmetric,
          arrayView1d< real64 const > const & normalTractionTolerance,
          arrayView2d< real64 > const & traction,
          arrayView1d< integer > const & fractureState )

  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      real64 const zero = LvArray::NumericLimits< real64 >::epsilon;

      real64 localPenalty[3][3]{};
      real64 localTractionNew[3]{};
      contactWrapper.updateTraction( oldDispJump[k],
                                     dispJump[k],
                                     penalty[k],
                                     traction[k],
                                     faceArea[k],
                                     symmetric,
                                     false,
                                     normalTractionTolerance[k],
                                     zero,
                                     localPenalty,
                                     localTractionNew,
                                     fractureState[k] );

      if( fractureState[k] == fields::contact::FractureState::Open )
      {

        LvArray::tensorOps::fill< 3 >( localTractionNew, 0.0 );
      }
      else if( LvArray::math::abs( localTractionNew[ 0 ] ) < normalTractionTolerance[k] )
      {
        LvArray::tensorOps::fill< 3 >( localTractionNew, 0.0 );
        fractureState[k] = fields::contact::FractureState::Slip;
      }

      LvArray::tensorOps::copy< 3 >( traction[k], localTractionNew );
      penalty[k][2] = -localPenalty[1][1];
      penalty[k][3] = -localPenalty[2][2];
      penalty[k][4] = -localPenalty[1][2];

    } );
  }

};

} // namespace SolidMechanicsALMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMKERNELSBASE_HPP_ */
