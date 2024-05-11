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
 * @file SolidMechanicsALMBubbleKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMBUBBLEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMBUBBLEKERNELS_HPP_

#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"
#include "SolidMechanicsALMKernelsHelper.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ALMBubbleKernels :
  public solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE,
                                                                             CONSTITUTIVE_TYPE,
                                                                             FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE,
                                                                                   CONSTITUTIVE_TYPE,
                                                                                   FE_TYPE >;

  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  static constexpr int numFacesPerElem = FE_TYPE::numFaces;

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

  ALMBubbleKernels( NodeManager const & nodeManager,
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
    m_bubbleDisp( faceManager.getField< fields::solidMechanics::totalBubbleDisplacement >() ),
    m_bDofNumber( bDofNumber ),
    m_bubbleElems( elementSubRegion.bubbleElementsList() ),
    m_elemsToFaces(elementSubRegion.faceElementsList() )
  {}

  struct StackVariables  // it's better not to inherit all the stack variable. There is a lot of unused ones.
  {
public:

    static constexpr int numUdofs = numNodesPerElem * 3;

    static constexpr int numBubbleUdofs = numFacesPerElem * 3;

    GEOS_HOST_DEVICE
    StackVariables():
      dispEqnRowIndices{},
      dispColIndices{},
      bEqnRowIndices{},
      bColIndices{},
      localRu{},
      localRb{},
      localAbb{{}},
      localAbu{{}},
      localAub{{}},
      bLocal{},
      uLocal{},
      X{{}},
      constitutiveStiffness{{}}
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

  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const kk,
              StackVariables & stack ) const
  {
    GEOS_MARK_FUNCTION;
    //std::cout << "Start Setup" << std::endl;

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
    //std::cout << "localFaceIndex" << localFaceIndex << std::endl;

    for( int i=0; i<3; ++i )
    {
      // need to grab the index.
      stack.bEqnRowIndices[i] = m_bDofNumber[localFaceIndex] + i - m_dofRankOffset;
      stack.bColIndices[i]    = m_bDofNumber[localFaceIndex] + i;
      stack.bLocal[ i ] = m_bubbleDisp[ localFaceIndex ][i];
    }
    //std::cout << "End Setup" << std::endl;
  }


  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const kk,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    GEOS_MARK_FUNCTION;

    localIndex k = m_bubbleElems[kk];
    constexpr int nUdof = numNodesPerElem*3;
    constexpr int nBubbleUdof = numFacesPerElem*3;

    real64 dBubbleNdX[ numFacesPerElem ][ 3 ];
    // It is needed only because I inserted a placeholder for calcGradFaceBubbleN in some finite elements
    LvArray::tensorOps::fill< numFacesPerElem, 3 >( dBubbleNdX, 0 );  //make 0
    //

    real64 detJ = m_finiteElementSpace.template calcGradFaceBubbleN( q, stack.X, dBubbleNdX );

    real64 dNdX[ numNodesPerElem ][ 3 ];
    detJ = m_finiteElementSpace.template calcGradN( q, stack.X, dNdX );

    m_constitutiveUpdate.getElasticStiffness( k, q, stack.constitutiveStiffness );

    real64 strainMatrix[6][nUdof];
    solidMechanicsALMKernelsHelper::assembleStrainOperator< 6, nUdof, numNodesPerElem >( strainMatrix, dNdX );

    real64 strainBubbleMatrix[6][nBubbleUdof];
    solidMechanicsALMKernelsHelper::assembleStrainOperator< 6, nBubbleUdof, numFacesPerElem >( strainBubbleMatrix, dBubbleNdX );

    // TODO: Use the following functions
    //BilinearFormUtilities::compute< displacementTestSpace,
    //                                displacementTrialSpace,
    //                                DifferentialOperator::SymmetricGradient,
    //                                DifferentialOperator::SymmetricGradient >
    //(
    //  stack.dLocalResidualMomentum_dDisplacement,
    //  dNdX,
    //  stack.stiffness, // fourth-order tensor handled via DiscretizationOps
    //  dNdX,
    //  -detJxW );

   
    //LinearFormUtilities::compute< displacementTestSpace,
    //                            DifferentialOperator::Identity >
    //(
    //stack.localResidualMomentum,
    //N,
    //stack.bodyForce,
    //detJxW );


    //for (int i=0; i<6; ++i)
    //{
    //  for (int j=0; j<nBubbleUdof; ++j)
    //  {
    //    std::cout << strainBubbleMatrix[i][j] << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //abort();

    real64 matBD[nBubbleUdof][6];
    real64 Abb_gauss[nBubbleUdof][nBubbleUdof], Abu_gauss[nBubbleUdof][nUdof], Aub_gauss[nUdof][nBubbleUdof];

    // transp(Bb)D
    LvArray::tensorOps::Rij_eq_AkiBkj< nBubbleUdof, 6, 6 >( matBD, strainBubbleMatrix, stack.constitutiveStiffness );

    // transp(Bb)DB
    LvArray::tensorOps::Rij_eq_AikBkj< nBubbleUdof, nUdof, 6 >( Abu_gauss, matBD, strainMatrix );

    // transp(Bb)DBb
    LvArray::tensorOps::Rij_eq_AikBkj< nBubbleUdof, nBubbleUdof, 6 >( Abb_gauss, matBD, strainBubbleMatrix );

    // transp(B)DBb
    tensorOps::transpose< nUdof, nBubbleUdof >( Aub_gauss, Abu_gauss );

    //GEOS_UNUSED_VAR(Aub_gauss, Abu_gauss, Abb_gauss, matBD, detJ);

    //real64 strainInc[6] = {0};
    //FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainInc );

    // multiply by determinant and add to element matrix
    LvArray::tensorOps::scaledAdd< nBubbleUdof, nBubbleUdof >( stack.localAbb, Abb_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< nBubbleUdof, nUdof >( stack.localAbu, Abu_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< nUdof, nBubbleUdof >( stack.localAub, Aub_gauss, -detJ );
  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const kk,
                   StackVariables & stack ) const
  {
    GEOS_MARK_FUNCTION;

    real64 maxForce = 0;
    constexpr int nUdof = numNodesPerElem*3;

    localIndex const parentFaceIndex = m_elemsToFaces[kk][1];


    real64 localAub[nUdof][3];
    real64 localAbu[3][nUdof];
    real64 localAbb[3][3];
    real64 localRb[3];
    real64 localRu[nUdof];
    for( localIndex i = 0; i < 3; ++i )
    {
      for( localIndex j = 0; j < nUdof; ++j )
      {
        localAub[j][i] = stack.localAub[parentFaceIndex*3+j][i];
      }
      for( localIndex j = 0; j < 3; ++j )
      {
        localAbb[j][i] = stack.localAbu[parentFaceIndex*3+j][parentFaceIndex*3+i];
      }
      localRb[i] = stack.localRb[i];
    }

    for( localIndex j = 0; j < nUdof; ++j )
    {
      for( localIndex i = 0; i < 3; ++i )
      {
        localAbu[i][j] = stack.localAbu[i][parentFaceIndex*3+j];

      }
      localRu[j] = stack.localRu[j];
    }


    // Compute the local residuals
    LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( localRb, localAbb, stack.bLocal );
    LvArray::tensorOps::Ri_add_AijBj< 3, nUdof >( localRb, localAbu, stack.uLocal );
    LvArray::tensorOps::Ri_add_AijBj< nUdof, 3 >( localRu, localAub, stack.bLocal );

    //real64 localAub[nUdof][3];
    //for( localIndex i = 0; i < nUdof; ++i )
    //{
    //  for( localIndex j = 0; j < 3; ++j )
    //  {
    //    localAub[i][j] = stack.localAub[i][parentFaceIndex*3+j];
    //  }
    //}

    for( localIndex i = 0; i < nUdof; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], localRu[i] );

      //real64 rowAub[3];
      ///for( localIndex j = 0; j < 3; ++j )
      //{
      //  rowAub[j] = stack.localAub[i][parentFaceIndex*3+j];
      //}

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.bColIndices,
                                                                              localAub[i],
                                                                              3 );
      //GEOS_UNUSED_VAR( rowAub, dof );

    }

    for( localIndex i=0; i < 3; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.bEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], localRb[i] );

      //real64 rowAbb[3];
      //for( localIndex j = 0; j < 3; ++j )
      //{
      //  rowAbb[j] = stack.localAbb[parentFaceIndex*3+i][parentFaceIndex*3+j];
      //}

      //GEOS_UNUSED_VAR( rowAbb, dof );

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

  arrayView2d< real64 const > const m_bubbleDisp;

  /// The global degree of freedom number of bubble
  arrayView1d< globalIndex const > const m_bDofNumber;

  arrayView1d< localIndex const > const m_bubbleElems;

  arrayView2d< localIndex const > const  m_elemsToFaces;

};

/// The factory used to construct a QuasiStatic kernel.
using ALMBubbleFactory = finiteElement::KernelFactory< ALMBubbleKernels,
                                                       arrayView1d< globalIndex const > const,
                                                       arrayView1d< globalIndex const > const,
                                                       globalIndex const,
                                                       CRSMatrixView< real64, globalIndex const > const,
                                                       arrayView1d< real64 > const,
                                                       real64 const,
                                                       real64 const (&) [3] >;

} // namespace SolidMechanicsALMBubbleKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMBUBBLEKERNELS_HPP_ */
