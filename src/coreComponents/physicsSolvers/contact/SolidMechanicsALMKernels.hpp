/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsALMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_

#include "SolidMechanicsALMKernelsBase.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

/**
 * @copydoc geos::finiteElement::ImplicitKernelBase
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ALM :
  public ALMKernelsBase< CONSTITUTIVE_TYPE,
                         FE_TYPE >
{
public:
  /// Alias for the base class.
  using Base = ALMKernelsBase< CONSTITUTIVE_TYPE,
                               FE_TYPE >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::m_elemsToFaces;
  using Base::m_faceToNodes;
  using Base::m_finiteElementSpace;
  using Base::m_dofNumber;
  using Base::m_bDofNumber;
  using Base::m_dofRankOffset;
  using Base::m_X;
  using Base::m_rotationMatrix;
  using Base::m_penalty;
  using Base::m_dispJump;
  using Base::m_oldDispJump;
  using Base::m_matrix;
  using Base::m_rhs;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::InterfaceKernelBase::InterfaceKernelBase
   */
  ALM( NodeManager const & nodeManager,
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
    m_traction( elementSubRegion.getField< fields::contact::traction >().toViewConst())
  {}

  //***************************************************************************

  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3 * 2;

    /// The number of lagrange multiplier dofs per element.
    static constexpr int numTdofs = 3;

    /// The number of bubble dofs per element.
    static constexpr int numBdofs = 3*2;

public:

    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dispEqnRowIndices{},
            dispColIndices{},
            bEqnRowIndices{},
            bColIndices{},
            localRu{},
            localRb{},
            localAutAtu{ {} },
      localAbtAtb{ {} },
      localAbtAtu{ {} },
      localAutAtb{ {} },
      tLocal{}
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex bEqnRowIndices[numBdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex bColIndices[numBdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local Rb residual vector.
    real64 localRb[numBdofs];

    /// C-array storage for the element local AutAtu matrix.
    real64 localAutAtu[numUdofs][numUdofs];

    /// C-array storage for the element local AbtAtb matrix.
    real64 localAbtAtb[numBdofs][numBdofs];

    /// C-array storage for the element local AbtAtu matrix.
    real64 localAbtAtu[numBdofs][numUdofs];

    /// C-array storage for the element local AbtAtu matrix.
    real64 localAutAtb[numUdofs][numBdofs];

    /// Stack storage for the element local lagrange multiplier vector
    real64 tLocal[numTdofs];

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

    constexpr int numTdofs = 3;

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
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[kn0]+i-m_dofRankOffset;
        stack.dispEqnRowIndices[shift + a*3+i] = m_dofNumber[kn1]+i-m_dofRankOffset;
        stack.dispColIndices[a*3+i] = m_dofNumber[kn0]+i;
        stack.dispColIndices[shift + a*3+i] = m_dofNumber[kn1]+i;
        stack.X[ a ][ i ] = m_X[ m_faceToNodes( kf0, permutation[ a ] ) ][ i ];
      }
    }

    for( int j=0; j<3; ++j )
    {
      for( int i=0; i<3; ++i )
      {
        stack.localRotationMatrix[ i ][ j ] = m_rotationMatrix( k, i, j );
      }
    }

    // The minus sign is consistent with the sign of the Jacobian
    stack.localPenalty[0][0] = -m_penalty( k, 0 );
    stack.localPenalty[1][1] = -m_penalty( k, 1 );
    stack.localPenalty[2][2] = -m_penalty( k, 1 );

    for( int i=0; i<numTdofs; ++i )
    {
      stack.tLocal[i] = m_traction( k, i );
      stack.dispJumpLocal[i] = m_dispJump( k, i );
      stack.oldDispJumpLocal[i] = m_oldDispJump( k, i );
    }

    for( int i=0; i<3; ++i )
    {
      // need to grab the index.
      stack.bEqnRowIndices[i]   = m_bDofNumber[kf0] + i - m_dofRankOffset;
      stack.bEqnRowIndices[3+i] = m_bDofNumber[kf1] + i - m_dofRankOffset;
      stack.bColIndices[i]      = m_bDofNumber[kf0] + i;
      stack.bColIndices[3+i]    = m_bDofNumber[kf1] + i;
    }
  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );

    constexpr int numUdofs = numNodesPerElem * 3 * 2;

    constexpr int numBdofs = 3*2;

    real64 matRRtAtu[3][numUdofs], matDRtAtu[3][numUdofs];
    real64 matRRtAtb[3][numBdofs], matDRtAtb[3][numBdofs];

    real64 dispJumpR[numUdofs];
    real64 oldDispJumpR[numUdofs];
    real64 tractionR[numUdofs];
    real64 dispJumpRb[numBdofs];
    real64 oldDispJumpRb[numBdofs];
    real64 tractionRb[numBdofs];

    // transp(R) * Atu
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, numUdofs, 3 >( matRRtAtu, stack.localRotationMatrix,
                                                         stack.localAtu );
    // transp(R) * Atb
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, numBdofs, 3 >( matRRtAtb, stack.localRotationMatrix,
                                                         stack.localAtb );

    // Compute the traction contribute of the local residuals
    LvArray::tensorOps::Ri_eq_AjiBj< numUdofs, 3 >( tractionR, matRRtAtu, stack.tLocal );
    LvArray::tensorOps::Ri_eq_AjiBj< numBdofs, 3 >( tractionRb, matRRtAtb, stack.tLocal );

    // D*RtAtu
    LvArray::tensorOps::Rij_eq_AikBkj< 3, numUdofs, 3 >( matDRtAtu, stack.localPenalty,
                                                         matRRtAtu );
    // D*RtAtb
    LvArray::tensorOps::Rij_eq_AikBkj< 3, numBdofs, 3 >( matDRtAtb, stack.localPenalty,
                                                         matRRtAtb );

    // R*RtAtu
    LvArray::tensorOps::Rij_eq_AikBkj< 3, numUdofs, 3 >( matRRtAtu, stack.localRotationMatrix,
                                                         matDRtAtu );
    // R*RtAtb
    LvArray::tensorOps::Rij_eq_AikBkj< 3, numBdofs, 3 >( matRRtAtb, stack.localRotationMatrix,
                                                         matDRtAtb );

    // transp(Atu)*RRtAtu
    LvArray::tensorOps::Rij_eq_AkiBkj< numUdofs, numUdofs, 3 >( stack.localAutAtu, stack.localAtu,
                                                                matRRtAtu );
    // transp(Atb)*RRtAtb
    LvArray::tensorOps::Rij_eq_AkiBkj< numBdofs, numBdofs, 3 >( stack.localAbtAtb, stack.localAtb,
                                                                matRRtAtb );

    // transp(Atb)*RRtAtu
    LvArray::tensorOps::Rij_eq_AkiBkj< numBdofs, numUdofs, 3 >( stack.localAbtAtu, stack.localAtb,
                                                                matRRtAtu );

    // transp(Atu)*RRtAtb
    LvArray::tensorOps::Rij_eq_AkiBkj< numUdofs, numBdofs, 3 >( stack.localAutAtb, stack.localAtu,
                                                                matRRtAtb );

    // Compute the local residuals
    LvArray::tensorOps::Ri_eq_AjiBj< numUdofs, 3 >( dispJumpR, matDRtAtu, stack.dispJumpLocal );
    LvArray::tensorOps::Ri_eq_AjiBj< numUdofs, 3 >( oldDispJumpR, matDRtAtu, stack.oldDispJumpLocal );
    LvArray::tensorOps::Ri_eq_AjiBj< numBdofs, 3 >( dispJumpRb, matDRtAtb, stack.dispJumpLocal );
    LvArray::tensorOps::Ri_eq_AjiBj< numBdofs, 3 >( oldDispJumpRb, matDRtAtb, stack.oldDispJumpLocal );

    LvArray::tensorOps::scaledAdd< numUdofs >( stack.localRu, tractionR, -1 );
    LvArray::tensorOps::scaledAdd< numUdofs >( stack.localRu, dispJumpR, 1 );
    LvArray::tensorOps::scaledAdd< numUdofs >( stack.localRu, oldDispJumpR, -1 );

    LvArray::tensorOps::scaledAdd< numBdofs >( stack.localRb, tractionRb, -1 );
    LvArray::tensorOps::scaledAdd< numBdofs >( stack.localRb, dispJumpRb, 1 );
    LvArray::tensorOps::scaledAdd< numBdofs >( stack.localRb, oldDispJumpRb, -1 );

    for( localIndex i=0; i < numUdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      // Is it necessary? Each row should be indepenedent
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );

      // Fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              stack.localAutAtu[i],
                                                                              numUdofs );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.bColIndices,
                                                                              stack.localAutAtb[i],
                                                                              numBdofs );
    }

    for( localIndex i=0; i < numBdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.bEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      // Is it necessary? Each row should be indepenedent
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRb[i] );

      // Fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.bColIndices,
                                                                              stack.localAbtAtb[i],
                                                                              numBdofs );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              stack.localAbtAtu[i],
                                                                              numUdofs );
    }

    return 0.0;
  }

protected:

  arrayView2d< real64 const > const m_traction;

};

/// The factory used to construct the kernel.
using ALMFactory = finiteElement::InterfaceKernelFactory< ALM,
                                                          arrayView1d< globalIndex const > const,
                                                          arrayView1d< globalIndex const > const,
                                                          globalIndex const,
                                                          CRSMatrixView< real64, globalIndex const > const,
                                                          arrayView1d< real64 > const,
                                                          real64 const,
                                                          arrayView1d< localIndex const > const >;


/**
 * @brief A struct to compute rotation matrices
 */
struct ComputeRotationMatricesKernel
{

  /**
   * @brief Launch the kernel function to comute rotation matrices
   * @tparam POLICY the type of policy used in the kernel launch
   * @param[in] size the size of the subregion
   * @param[in] faceNormal the array of array containing the face to nodes map
   * @param[in] elemsToFaces the array of array containing the element to faces map
   * @param[out] rotationMatrix the array containing the rotation matrices
   */
  template< typename POLICY >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const > const & faceNormal,
          ArrayOfArraysView< localIndex const > const & elemsToFaces,
          arrayView3d< real64 > const & rotationMatrix )
  {

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      localIndex const & f0 = elemsToFaces[k][0];
      localIndex const & f1 = elemsToFaces[k][1];

      real64 Nbar[3];
      Nbar[0] = faceNormal[f0][0] - faceNormal[f1][0];
      Nbar[1] = faceNormal[f0][1] - faceNormal[f1][1];
      Nbar[2] = faceNormal[f0][2] - faceNormal[f1][2];

      LvArray::tensorOps::normalize< 3 >( Nbar );
      computationalGeometry::RotationMatrix_3D( Nbar, rotationMatrix[k] );

    } );
  }

};

} // namespace SolidMechanicsALMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELS_HPP_ */
