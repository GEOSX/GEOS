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
 * @file SolidMechanicsALMUpdateKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMUPDATEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMUPDATEKERNELS_HPP_

#include "SolidMechanicsALMKernelsBase.hpp"

namespace geos
{

namespace solidMechanicsConformingContactKernels
{

/**
 * @copydoc geos::finiteElement::ImplicitKernelBase
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class DispJumpUpdate :
  public ALMKernelsBase< CONSTITUTIVE_TYPE,
                         FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = ALMKernelsBase< CONSTITUTIVE_TYPE,
                               FE_TYPE >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::m_X;
  using Base::m_finiteElementSpace;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_elemsToFaces;
  using Base::m_faceToNodes;
  using Base::m_rotationMatrix;
  using Base::m_dispJump;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::InterfaceKernelBase::InterfaceKernelBase
   */
  DispJumpUpdate( NodeManager const & nodeManager,
                 EdgeManager const & edgeManager,
                 FaceManager const & faceManager,
                 localIndex const targetRegionIndex,
                 FaceElementSubRegion & elementSubRegion,
                 FE_TYPE const & finiteElementSpace,
                 CONSTITUTIVE_TYPE & inputConstitutiveType,
                 arrayView1d< globalIndex const > const uDofNumber,
                 arrayView1d< globalIndex const > const bDofNumber,
                 globalIndex const rankOffset,
                 CRSMatrixView< real64, globalIndex const > const inputMatrix,
                 arrayView1d< real64 > const inputRhs,
                 real64 const inputDt,
                 arrayView1d< localIndex const > const & faceElementList ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          uDofNumber,
          bDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt,
          faceElementList ),
    m_displacement( nodeManager.getField< fields::solidMechanics::totalDisplacement >()),
    m_bubbleDisp( faceManager.getField< fields::solidMechanics::totalBubbleDisplacement >() ),
    m_incrDisp( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
    m_incrBubbleDisp( faceManager.getField< fields::solidMechanics::incrementalBubbleDisplacement >() ),
    m_deltaDispJump( elementSubRegion.getField< fields::contact::deltaDispJump >().toView() )
  {}

  //***************************************************************************

  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3 * 2;

    /// The number of bubble dofs per element.
    static constexpr int numBdofs = 3 * 2;

    /// The number of lagrange multiplier dofs per element.
    static constexpr int numTdofs = 3;

public:

    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            uLocal{},
            bLocal{},
            duLocal{},
            dbLocal{},
            deltaDispJumpLocal{}
    {}

    /// Stack storage for the element local displacement vector
    real64 uLocal[numUdofs];

    /// Stack storage for the element local bubble displacement vector
    real64 bLocal[numBdofs];

    /// Stack storage for the element local incremental displacement vector
    real64 duLocal[numUdofs];

    /// Stack storage for the element local incremental bubble displacement vector
    real64 dbLocal[numBdofs];

    /// Stack storage for the element local delta displacement jump vector
    real64 deltaDispJumpLocal[numTdofs];

  };

  //***************************************************************************

  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
  }
  //END_kernelLauncher

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::InterfaceKernelBase::setup
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    constexpr int shift = numNodesPerElem * 3;

    int permutation[numNodesPerElem];
    m_finiteElementSpace.getPermutation( permutation );

    localIndex const kf0 = m_elemsToFaces[k][0];
    localIndex const kf1 = m_elemsToFaces[k][1];
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const kn0 = m_faceToNodes( kf0, a );
      localIndex const kn1 = m_faceToNodes( kf1, a );

      for( int i=0; i<3; ++i )
      {
        stack.X[ a ][ i ] = m_X[ m_faceToNodes( kf0, permutation[ a ] ) ][ i ];
        stack.uLocal[a*3+i] = m_displacement[kn0][i];
        stack.uLocal[shift + a*3+i] = m_displacement[kn1][i];
        stack.duLocal[a*3+i] = m_incrDisp[kn0][i];
        stack.duLocal[shift + a*3+i] = m_incrDisp[kn1][i];
      }
    }

    for( int j=0; j<3; ++j )
    {
      for( int i=0; i<3; ++i )
      {
        stack.localRotationMatrix[ i ][ j ] = m_rotationMatrix( k, i, j );
      }
    }

    for( int i=0; i<3; ++i )
    {
      stack.bLocal[ i ] = m_bubbleDisp[ kf0 ][i];
      stack.bLocal[ 3 + i ] = m_bubbleDisp[ kf1 ][i];
      stack.dbLocal[ i ] = m_incrBubbleDisp[ kf0 ][i];
      stack.dbLocal[ 3 + i ] = m_incrBubbleDisp[ kf1 ][i];
    }

  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {

    constexpr int numUdofs = numNodesPerElem * 3 * 2;

    constexpr int numBdofs = 3 * 2;

    real64 matRtAtu[3][numUdofs];
    real64 matRtAtb[3][numBdofs];

    // transp(R) * Atu
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, numUdofs, 3 >( matRtAtu, stack.localRotationMatrix,
                                                         stack.localAtu );
    // transp(R) * Atb
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, numBdofs, 3 >( matRtAtb, stack.localRotationMatrix,
                                                         stack.localAtb );

    // Compute the node contribute of the displacement and delta displacement jump
    LvArray::tensorOps::Ri_eq_AijBj< 3, numUdofs >( stack.dispJumpLocal, matRtAtu, stack.uLocal );
    LvArray::tensorOps::Ri_eq_AijBj< 3, numUdofs >( stack.deltaDispJumpLocal, matRtAtu, stack.duLocal );

    // Compute the bubble contribute of the displacement and delta displacement jump
    LvArray::tensorOps::Ri_add_AijBj< 3, numBdofs >( stack.dispJumpLocal, matRtAtb, stack.bLocal );
    LvArray::tensorOps::Ri_add_AijBj< 3, numBdofs >( stack.deltaDispJumpLocal, matRtAtb, stack.dbLocal );

    // Store the results
    for( int i=0; i<3; ++i )
    {
      m_dispJump[ k ][ i ] = stack.dispJumpLocal[ i ];
      m_deltaDispJump[ k ][ i ] = stack.deltaDispJumpLocal[ i ];
    }

    return 0.0;
  }

protected:

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_displacement;

  /// The rank-global bubble displacement array.
  arrayView2d< real64 const > const m_bubbleDisp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_incrDisp;

  /// The rank-global incremental bubble displacement array.
  arrayView2d< real64 const > const m_incrBubbleDisp;

  /// The rank-global delta displacement jump array.
  arrayView2d< real64 > const m_deltaDispJump;

};

using DispJumpUpdateFactory = finiteElement::InterfaceKernelFactory< DispJumpUpdate,
                                                                    arrayView1d< globalIndex const > const,
                                                                    arrayView1d< globalIndex const > const,
                                                                    globalIndex const,
                                                                    CRSMatrixView< real64, globalIndex const > const,
                                                                    arrayView1d< real64 > const,
                                                                    real64 const,
                                                                    arrayView1d< localIndex const > const >;

} // namespace SolidMechanicsALMKernels

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMUPDATEKERNELS_HPP_ */
