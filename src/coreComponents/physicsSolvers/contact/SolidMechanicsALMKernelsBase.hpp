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

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSBASE_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"
#include "SolidMechanicsALMKernelsHelper.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

/**
 * @brief Implements kernels for ALM.
 * @copydoc geos::finiteElement::InterfaceKernelBase
 *
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ALMKernelsBase :
  public finiteElement::InterfaceKernelBase< CONSTITUTIVE_TYPE,
                                             FE_TYPE,
                                             3, 3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::InterfaceKernelBase< CONSTITUTIVE_TYPE,
                                                   FE_TYPE,
                                                   3, 3 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_finiteElementSpace;
  using Base::m_matrix;
  using Base::m_rhs;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::InterfaceKernelBase::InterfaceKernelBase
   */
  ALMKernelsBase( NodeManager const & nodeManager,
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
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt ),
    m_X( nodeManager.referencePosition()),
    m_faceToNodes( faceManager.nodeList().toViewConst()),
    m_elemsToFaces( elementSubRegion.faceList().toViewConst()),
    m_faceElementList( faceElementList ),
    m_bDofNumber( bDofNumber ),
    m_rotationMatrix( elementSubRegion.getField< fields::contact::rotationMatrix >().toViewConst()),
    m_dispJump( elementSubRegion.getField< fields::contact::dispJump >().toView() ),
    m_oldDispJump( elementSubRegion.getField< fields::contact::oldDispJump >().toViewConst() ),
    m_penalty( elementSubRegion.getField< fields::contact::iterativePenalty >().toViewConst() )
  {}

  //***************************************************************************
  /**
   * @copydoc finiteElement::InterfaceKernelBase::StackVariables
   */
  struct StackVariables
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
      localAtu{ {} },
      localAtb{ {} },
      localRotationMatrix{ {} },
      localPenalty{ {} },
      dispJumpLocal{},
      oldDispJumpLocal{},
      X{ {} }
    {}

    /// C-array storage for the element local Atu matrix.
    real64 localAtu[numTdofs][numUdofs];

    /// C-array storage for the element local Atb matrix.
    real64 localAtb[numTdofs][numBdofs];

    /// C-array storage for rotation matrix
    real64 localRotationMatrix[3][3];

    /// C-array storage for penalty matrix
    real64 localPenalty[3][3];

    /// Stack storage for the element local displacement jump vector
    real64 dispJumpLocal[numTdofs];

    /// Stack storage for the element local old displacement jump vector
    real64 oldDispJumpLocal[numTdofs];

    /// local nodal coordinates
    real64 X[ numNodesPerElem ][ 3 ];

  };

  //***************************************************************************

  /**
   * @copydoc ::geos::finiteElement::InterfaceKernelBase::kernelLaunch
   *
   */
  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    GEOS_UNUSED_VAR( numElems );

    // Define a RAJA reduction variable to get the maximum residual contribution.
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( kernelComponent.m_faceElementList.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const i )
    {

      localIndex k = kernelComponent.m_faceElementList[i];
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      maxResidual.max( kernelComponent.complete( k, stack ) );
    } );

    return maxResidual.get();
  }
  //END_kernelLauncher

  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );

    constexpr int numUdofs = numNodesPerElem * 3 * 2;

    constexpr int numTdofs = 3;

    constexpr int numBdofs = 3*2;

    real64 const detJ = m_finiteElementSpace.transformedQuadratureWeight( q, stack.X );

    real64 N[ numNodesPerElem ];
    m_finiteElementSpace.calcN( q, N );

    real64 BubbleN[1];
    // Next line is needed because I only inserted a placeholder for calcBubbleN in some finite elements
    BubbleN[0]=0.0;  //make 0
    constexpr int bperm[1] = {0};
    m_finiteElementSpace.calcBubbleN( q, BubbleN );

    int permutation[numNodesPerElem];
    m_finiteElementSpace.getPermutation( permutation );

    // TODO: Try using bilinear utilities to perform these two operations
    solidMechanicsALMKernelsHelper::accumulateAtuLocalOperator< numTdofs,
                                                                numUdofs,
                                                                numNodesPerElem >( stack.localAtu,
                                                                                   N,
                                                                                   permutation,
                                                                                   detJ );

    solidMechanicsALMKernelsHelper::accumulateAtuLocalOperator< numTdofs,
                                                                numBdofs,
                                                                1 >( stack.localAtb,
                                                                     BubbleN,
                                                                     bperm,
                                                                     detJ );
  }

protected:

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array of array containing the face to node map.
  ArrayOfArraysView< localIndex const > const m_faceToNodes;

  /// The array of array containing the element to face map.
  ArrayOfArraysView< localIndex const > const m_elemsToFaces;

  /// The array containing the list of face element of the same type.
  arrayView1d< localIndex const > const m_faceElementList;

  /// The global degree of freedom number of bubble.
  arrayView1d< globalIndex const > const m_bDofNumber;

  /// The array containing the rotation matrix for each element.
  arrayView3d< real64 const > const m_rotationMatrix;

  /// The array containing the displacement jump.
  arrayView2d< real64 > const m_dispJump;

  /// The array containing the displacement jump of previus time step.
  arrayView2d< real64 const > const m_oldDispJump;

  /// The array containing the penalty coefficients for each element.
  arrayView2d< real64 const > const m_penalty;

};

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
          arrayView1d< real64 const > const & area,
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
                                        normalDisplacementTolerance[k]*area[k],
                                        slidingTolerance[k]*area[k],
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


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSBASE_HPP_ */
