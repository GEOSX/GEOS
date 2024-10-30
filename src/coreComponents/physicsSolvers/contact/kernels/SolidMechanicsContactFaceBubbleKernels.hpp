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
 * @file SolidMechanicsContactFaceBubbleKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMBUBBLEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMBUBBLEKERNELS_HPP_

#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"
#include "SolidMechanicsALMKernelsHelper.hpp"

// TODO: Use the bilinear form utilities
//#include "finiteElement/BilinearFormUtilities.hpp"

namespace geos
{

namespace solidMechanicsConformingContactKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class FaceBubbleKernels :
  public solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE,
                                                                             CONSTITUTIVE_TYPE,
                                                                             FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE,
                                                                                   CONSTITUTIVE_TYPE,
                                                                                   FE_TYPE >;

  /// Number of nodes per element, which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// Compile time value for the number of faces per element.
  static constexpr int numFacesPerElem = FE_TYPE::numFaces;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::m_X;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;
  using Base::m_constitutiveUpdate;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_disp;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   */
  FaceBubbleKernels( NodeManager const & nodeManager,
                    EdgeManager const & edgeManager,
                    FaceManager const & faceManager,
                    localIndex const targetRegionIndex,
                    SUBREGION_TYPE const & elementSubRegion,
                    FE_TYPE const & finiteElementSpace,
                    CONSTITUTIVE_TYPE & inputConstitutiveType,
                    arrayView1d< globalIndex const > const uDofNumber,
                    arrayView1d< globalIndex const > const bDofNumber,
                    globalIndex const rankOffset,
                    CRSMatrixView< real64, globalIndex const > const inputMatrix,
                    arrayView1d< real64 > const inputRhs,
                    real64 const inputDt,
                    real64 const (&inputGravityVector)[3] ):
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
          inputDt,
          inputGravityVector ),
    m_bubbleDisp( faceManager.getField< fields::solidMechanics::totalBubbleDisplacement >().toViewConst() ),
    m_bDofNumber( bDofNumber ),
    m_bubbleElems( elementSubRegion.bubbleElementsList() ),
    m_elemsToFaces( elementSubRegion.faceElementsList() )
  {}

  //***************************************************************************

  /**
   * @copydoc finiteElement::ImplicitKernelBase::StackVariables
   */
  struct StackVariables
  {
public:
    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;

    /// The number of jump dofs per element.
    static constexpr int numBubbleUdofs = numFacesPerElem * 3;

    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      dispEqnRowIndices{},
      dispColIndices{},
      bEqnRowIndices{},
      bColIndices{},
      localRu{},
      localRb{},
      localAbb{ {} },
      localAbu{ {} },
      localAub{ {} },
      bLocal{},
      uLocal{},
      X{ {} },
      constitutiveStiffness{ {} }
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex bEqnRowIndices[3];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex bColIndices[3];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local Rb residual vector.
    real64 localRb[numBubbleUdofs];

    /// C-array storage for the element local Abb matrix.
    real64 localAbb[numBubbleUdofs][numBubbleUdofs];

    /// C-array storage for the element local Abu matrix.
    real64 localAbu[numBubbleUdofs][numUdofs];

    /// C-array storage for the element local Aub matrix.
    real64 localAub[numUdofs][numBubbleUdofs];

    /// Stack storage for the element local bubble vector
    real64 bLocal[3];

    /// Stack storage for the element displacement vector.
    real64 uLocal[numUdofs];

    /// local nodal coordinates
    real64 X[ numNodesPerElem ][ 3 ];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[ 6 ][ 6 ];

  };

  //***************************************************************************

  /**
   * @copydoc ::geos::finiteElement::KernelBase::kernelLaunch
   *
   * @detail it uses the kernelLaunch interface of KernelBase but it only launches the kernel
   * on the set of elements that have bubble dof within the subregion.
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

    forAll< POLICY >( kernelComponent.m_bubbleElems.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( i, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( i, q, stack );
      }
      maxResidual.max( kernelComponent.complete( i, stack ) );
    } );

    return maxResidual.get();
  }
  //END_kernelLauncher

  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const kk,
              StackVariables & stack ) const
  {

    localIndex k = m_bubbleElems[kk];

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[localNodeIndex]+i-m_dofRankOffset;
        stack.dispColIndices[a*3+i]    = m_dofNumber[localNodeIndex]+i;
        stack.X[ a ][ i ] = m_X[ localNodeIndex ][ i ];
        stack.uLocal[ a*3 + i ] = m_disp[localNodeIndex][i];
      }
    }

    localIndex const localFaceIndex = m_elemsToFaces[kk][0];

    for( int i=0; i<3; ++i )
    {
      // need to grab the index.
      stack.bEqnRowIndices[i] = m_bDofNumber[localFaceIndex] + i - m_dofRankOffset;
      stack.bColIndices[i]    = m_bDofNumber[localFaceIndex] + i;
      stack.bLocal[ i ] = m_bubbleDisp[ localFaceIndex ][i];
    }
  }


  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const kk,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    localIndex k = m_bubbleElems[kk];
    constexpr int nUdof = numNodesPerElem*3;
    constexpr int nBubbleUdof = numFacesPerElem*3;

    real64 dBubbleNdX[ numFacesPerElem ][ 3 ];
    // Next line is needed because I only inserted a placeholder for calcGradFaceBubbleN in some finite elements
    LvArray::tensorOps::fill< numFacesPerElem, 3 >( dBubbleNdX, 0 );  //make 0

    real64 detJ = m_finiteElementSpace.calcGradFaceBubbleN( q, stack.X, dBubbleNdX );

    real64 dNdX[ numNodesPerElem ][ 3 ];
    detJ = m_finiteElementSpace.calcGradN( q, stack.X, dNdX );

    m_constitutiveUpdate.getElasticStiffness( k, q, stack.constitutiveStiffness );

    real64 strainMatrix[6][nUdof];
    solidMechanicsALMKernelsHelper::assembleStrainOperator< 6, nUdof, numNodesPerElem >( strainMatrix, dNdX );

    real64 strainBubbleMatrix[6][nBubbleUdof];
    solidMechanicsALMKernelsHelper::assembleStrainOperator< 6, nBubbleUdof, numFacesPerElem >( strainBubbleMatrix, dBubbleNdX );

    // TODO: It would be nice use BilinearFormUtilities::compute

    real64 matBD[nBubbleUdof][6];
    real64 Abb_gauss[nBubbleUdof][nBubbleUdof], Abu_gauss[nBubbleUdof][nUdof], Aub_gauss[nUdof][nBubbleUdof];

    // transp(Bb)D
    LvArray::tensorOps::Rij_eq_AkiBkj< nBubbleUdof, 6, 6 >( matBD, strainBubbleMatrix, stack.constitutiveStiffness );

    // transp(Bb)DB
    LvArray::tensorOps::Rij_eq_AikBkj< nBubbleUdof, nUdof, 6 >( Abu_gauss, matBD, strainMatrix );

    // transp(Bb)DBb
    LvArray::tensorOps::Rij_eq_AikBkj< nBubbleUdof, nBubbleUdof, 6 >( Abb_gauss, matBD, strainBubbleMatrix );

    // transp(B)DBb
    LvArray::tensorOps::transpose< nUdof, nBubbleUdof >( Aub_gauss, Abu_gauss );

    // multiply by determinant and add to element matrix
    LvArray::tensorOps::scaledAdd< nBubbleUdof, nBubbleUdof >( stack.localAbb, Abb_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< nBubbleUdof, nUdof >( stack.localAbu, Abu_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< nUdof, nBubbleUdof >( stack.localAub, Aub_gauss, -detJ );

    // Compute the initial stress
    // The following block assumes a linear elastic constitutive model
    real64 rb_gauss[nBubbleUdof]{};
    real64 strain[6] = {0};
    LvArray::tensorOps::Ri_eq_AijBj< 6, nUdof >( strain, strainMatrix, stack.uLocal );

    real64 initStressLocal[ 6 ] = {0};
    LvArray::tensorOps::Ri_eq_AijBj< 6, 6 >( initStressLocal, stack.constitutiveStiffness, strain );
    for( localIndex c = 0; c < 6; ++c )
    {
      initStressLocal[ c ] -= m_constitutiveUpdate.m_newStress( k, q, c );
    }

    LvArray::tensorOps::Ri_eq_AjiBj< nBubbleUdof, 6 >( rb_gauss, strainBubbleMatrix, initStressLocal );
    LvArray::tensorOps::scaledAdd< nBubbleUdof >( stack.localRb, rb_gauss, detJ );

  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const kk,
                   StackVariables & stack ) const
  {

    real64 maxForce = 0;
    constexpr int nUdof = numNodesPerElem*3;

    localIndex const parentFaceIndex = m_elemsToFaces[kk][1];

    // Extract only the submatrix and known term corresponding to the index of the local face
    // on which the bubble function was applied.
    real64 localAub[nUdof][3];
    real64 localAbu[3][nUdof];
    real64 localAbb[3][3];
    real64 localRb[3];
    for( localIndex i = 0; i < 3; ++i )
    {
      for( localIndex j = 0; j < nUdof; ++j )
      {
        localAub[j][i] = stack.localAub[j][parentFaceIndex*3+i];
      }
      for( localIndex j = 0; j < 3; ++j )
      {
        localAbb[j][i] = stack.localAbb[parentFaceIndex*3+j][parentFaceIndex*3+i];
      }
      localRb[i] = stack.localRb[parentFaceIndex*3+i];
    }

    for( localIndex j = 0; j < nUdof; ++j )
    {
      for( localIndex i = 0; i < 3; ++i )
      {
        localAbu[i][j] = stack.localAbu[parentFaceIndex*3+i][j];
      }
    }

    // Compute the local residuals
    LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( localRb, localAbb, stack.bLocal );
    LvArray::tensorOps::Ri_add_AijBj< 3, nUdof >( localRb, localAbu, stack.uLocal );
    LvArray::tensorOps::Ri_add_AijBj< nUdof, 3 >( stack.localRu, localAub, stack.bLocal );

    for( localIndex i = 0; i < nUdof; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );

      // fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.bColIndices,
                                                                              localAub[i],
                                                                              3 );

    }

    for( localIndex i=0; i < 3; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.bEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], localRb[i] );

      // fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.bColIndices,
                                                                              localAbb[i],
                                                                              3 );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              localAbu[i],
                                                                              numNodesPerElem*3 );
    }

    return maxForce;
  }

protected:

  /// The array containing the bubble displacement.
  arrayView2d< real64 const > const m_bubbleDisp;

  /// The global degree of freedom number of bubble
  arrayView1d< globalIndex const > const m_bDofNumber;

  /// The array containing the list of bubble elements.
  arrayView1d< localIndex const > const m_bubbleElems;

  /// The array containing the element to bubble face map.
  /// The bubble face is the face of the element to which the bubble is applied.
  /// Both the local and global face index are stored in this array.
  arrayView2d< localIndex const > const m_elemsToFaces;

};

/// The factory used to construct a QuasiStatic kernel.
using FaceBubbleFactory = finiteElement::KernelFactory< FaceBubbleKernels,
                                                       arrayView1d< globalIndex const > const,
                                                       arrayView1d< globalIndex const > const,
                                                       globalIndex const,
                                                       CRSMatrixView< real64, globalIndex const > const,
                                                       arrayView1d< real64 > const,
                                                       real64 const,
                                                       real64 const (&) [3] >;

} // namespace SolidMechanicsContactFaceBubbleKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMBUBBLEKERNELS_HPP_ */
