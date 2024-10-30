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

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSLAGRANGECONTACTKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSLAGRANGECONTACTKERNELS_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"
#include "SolidMechanicsConformingContactKernelsBase.hpp"

namespace geos
{

namespace solidMechanicsLagrangeContactKernels
{

/**
 * @copydoc geos::finiteElement::ImplicitKernelBase
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class LagrangeContact :
  public solidMechanicsConformingContactKernels::ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                                                               FE_TYPE >
{
public:
  /// Alias for the base class.
  using Base = solidMechanicsConformingContactKernels::ConformingContactKernelsBase< CONSTITUTIVE_TYPE,
                                                                                     FE_TYPE >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  using Base::numUdofs;
  using Base::numTdofs;
  using Base::numBdofs;

  using Base::m_elemsToFaces;
  using Base::m_faceToNodes;
  using Base::m_finiteElementSpace;
  using Base::m_constitutiveUpdate;
  using Base::m_dofNumber;
  using Base::m_bDofNumber;
  using Base::m_dofRankOffset;
  using Base::m_X;
  using Base::m_rotationMatrix;
  using Base::m_dispJump;
  using Base::m_oldDispJump;
  using Base::m_matrix;
  using Base::m_rhs;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::InterfaceKernelBase::InterfaceKernelBase
   */
  LagrangeContact( NodeManager const & nodeManager,
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
                   arrayView1d< localIndex const > const & faceElementList,
                   string const tractionDofKey ):
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
    m_traction( elementSubRegion.getField< fields::contact::traction >().toViewConst() ),
    m_tDofNumber( elementSubRegion.getReference< globalIndex_array >( tractionDofKey ) )
  {}

  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
                                       dispEqnRowIndices{},
                                       dispColIndices{},
                                       bEqnRowIndices{},
                                       bColIndices{},
                                       tColIndices{},
                                       localRu{},
                                       localRb{},
                                       localRt{},
                                       localAtt{ {} },
      localAut{ {} },
      localAbt{ {} }
    {}

    /// C-array storage for the element local row degrees of freedom.
    localIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local row degrees of freedom.
    localIndex bEqnRowIndices[numBdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex bColIndices[numBdofs];

    /// C-array storage for the traction local row degrees of freedom.
    localIndex tEqnRowIndices[numTdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex tColIndices[numTdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local Rb residual vector.
    real64 localRb[numBdofs];

    /// C-array storage for the element local Rt residual vector.
    real64 localRt[numTdofs];

    /// C-array storage for the element local Att matrix.
    real64 localAtt[numTdofs][numTdofs];

    /// C-array storage for the element local Aut matrix.
    real64 localAut[numUdofs][numTdofs];

    /// C-array storage for the element local Abt matrix.
    real64 localAbt[numBdofs][numTdofs];
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

    for( int i=0; i<numTdofs; ++i )
    {
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

    for( int i=0; i<3; ++i )
    {
      stack.tEqnRowIndices[i]   = m_tDofNumber[k] + i - m_dofRankOffset;
      stack.tColIndices[i]      = m_tDofNumber[k] + i;
    }
  }

  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    Base::quadraturePointKernel( k, q, stack, [ =, &stack ] GEOS_HOST_DEVICE ( real64 const detJ )
    {
      // This will vary depending on the state.
      stack.localRt[0] += detJ * m_dispJump[k][0];
      stack.localRt[1] += detJ * ( m_dispJump[k][1] - m_oldDispJump[k][1] );
      stack.localRt[2] += detJ * ( m_dispJump[k][2] - m_oldDispJump[k][2] );
    } );
  }


  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );
    constexpr real64 zero = LvArray::NumericLimits< real64 >::epsilon;

    real64 tractionR[numUdofs];
    real64 tractionRb[numBdofs];

    // transp(R) * Atu
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, numUdofs, 3 >( stack.localAtu, stack.localRotationMatrix, stack.localAtu );
    // transp(R) * Atb
    LvArray::tensorOps::Rij_eq_AkiBkj< 3, numBdofs, 3 >( stack.localAtb, stack.localRotationMatrix, stack.localAtb );

    LvArray::tensorOps::transpose< numUdofs, numTdofs >( stack.localAut, stack.localAtu );
    LvArray::tensorOps::transpose< numBdofs, numTdofs >( stack.localAbt, stack.localAtb );

    // Compute the traction contribute of the local residuals
    LvArray::tensorOps::Ri_eq_AijBj< numUdofs, 3 >( tractionR, stack.localAut, m_traction[k] );
    LvArray::tensorOps::Ri_eq_AijBj< numBdofs, 3 >( tractionRb, stack.localAbt, m_traction[k] );

    // Compute the local residuals
    // Force Balance for nodal displacement dofs
    LvArray::tensorOps::scaledAdd< numUdofs >( stack.localRu, tractionR, 1.0 );
    // Force Balance for the bubble dofs
    LvArray::tensorOps::scaledAdd< numBdofs >( stack.localRb, tractionRb, 1.0 );

    fillGlobalMatrix( k, stack );

    return 0.0;
  }

protected:

  arrayView2d< real64 const > const m_traction;

  arrayView1d< localIndex const > const m_tDofNumber;

  /**
   * @brief Create the list of finite elements of the same type
   *   for each FaceElementSubRegion (Triangle or Quadrilateral)
   *   and of the same fracture state (Stick or Slip).
   * @param domain The physical domain object
   */
  void updateStickSlipList( DomainPartition const & domain );

  /**
   * @brief Create the list of finite elements of the same type
   *   for each FaceElementSubRegion (Triangle or Quadrilateral).
   * @param domain The physical domain object
   */
  void createFaceTypeList( DomainPartition const & domain );

  /**
   * @brief Create the list of elements belonging to CellElementSubRegion
   *  that are enriched with the bubble basis functions
   * @param domain The physical domain object
   */
  void createBubbleCellList( DomainPartition & domain ) const;

private:

  void fillGlobalMatrix( localIndex const k,
                         StackVariables & stack ) const
  {

    for( localIndex i=0; i < numTdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.tEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      // TODO: May not need to be an atomic operation
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRt[i] );

      // Fill in matrix block Att
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.tColIndices,
                                                                              stack.localAtt[i],
                                                                              numTdofs );

      // Fill in matrix block Atu
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              stack.localAtu[i],
                                                                              numUdofs );

      // Fill in matrix block Atb
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.bColIndices,
                                                                              stack.localAtb[i],
                                                                              numBdofs );
    }

    for( localIndex i=0; i < numUdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      // Is it necessary? Each row should be indepenedent
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );

      // Fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.tColIndices,
                                                                              stack.localAut[i],
                                                                              numTdofs );

    }

    for( localIndex i=0; i < numBdofs; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.bEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      // Is it necessary? Each row should be indepenedent
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRb[i] );

      // Fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.tColIndices,
                                                                              stack.localAbt[i],
                                                                              numTdofs );
    }
  }

};

/// The factory used to construct the kernel.
using LagrangeContactFactory = finiteElement::InterfaceKernelFactory< LagrangeContact,
                                                                      arrayView1d< globalIndex const > const,
                                                                      arrayView1d< globalIndex const > const,
                                                                      globalIndex const,
                                                                      CRSMatrixView< real64, globalIndex const > const,
                                                                      arrayView1d< real64 > const,
                                                                      real64 const,
                                                                      arrayView1d< localIndex const > const,
                                                                      string const >;

} // namespace solidMechanicsLagrangeContactKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSLAGRANGECONTACTKERNELS_HPP_ */
