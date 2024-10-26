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
 * @file SolidMechanicsEFEMKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEFEMKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEFEMKERNELS_HPP_

#include "SolidMechanicsEFEMKernelsBase.hpp"

namespace geos
{

namespace solidMechanicsEFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class EFEM :
  public EFEMKernelsBase< SUBREGION_TYPE,
                          CONSTITUTIVE_TYPE,
                          FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = EFEMKernelsBase< SUBREGION_TYPE,
                                CONSTITUTIVE_TYPE,
                                FE_TYPE >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_X;
  using Base::m_disp;
  using Base::m_uhat;
  using Base:: m_w;
  using Base::m_tractionVec;
  using Base::m_dTraction_dJump;
  using Base::m_nVec;
  using Base::m_tVec1;
  using Base::m_tVec2;
  using Base::m_surfaceCenter;
  using Base::m_surfaceArea;
  using Base::m_elementVolume;
  using Base::m_fracturedElems;
  using Base::m_cellsToEmbeddedSurfaces;
  using Base::m_dt;


  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  EFEM( NodeManager const & nodeManager,
        EdgeManager const & edgeManager,
        FaceManager const & faceManager,
        localIndex const targetRegionIndex,
        SUBREGION_TYPE const & elementSubRegion,
        FE_TYPE const & finiteElementSpace,
        CONSTITUTIVE_TYPE & inputConstitutiveType,
        EmbeddedSurfaceSubRegion & embeddedSurfSubRegion,
        arrayView1d< globalIndex const > const uDofNumber,
        arrayView1d< globalIndex const > const wDofNumber,
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
          embeddedSurfSubRegion,
          uDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt,
          inputGravityVector ),
    m_wDofNumber( wDofNumber )
  {}

  //***************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    /// The number of jump dofs per element.
    static constexpr int numWdofs = 3;

    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            jumpEqnRowIndices{ 0 },
      jumpColIndices{ 0 }
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex jumpEqnRowIndices[numWdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex jumpColIndices[numWdofs];
  };
  //***************************************************************************


  /**
   * @copydoc ::geos::finiteElement::KernelBase::kernelLaunch
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
    return Base::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernelComponent );
  }

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::ImplicitKernelBase::setup
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {

    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

    stack.hInv = m_surfaceArea[embSurfIndex] / m_elementVolume[k];
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

    for( int i=0; i<3; ++i )
    {
      // need to grab the index.
      stack.jumpEqnRowIndices[i] = m_wDofNumber[embSurfIndex] + i - m_dofRankOffset;
      stack.jumpColIndices[i]    = m_wDofNumber[embSurfIndex] + i;
      stack.wLocal[ i ] = m_w[ embSurfIndex ][i];
      stack.tractionVec[ i ] = m_tractionVec[ embSurfIndex ][i] * m_surfaceArea[embSurfIndex];
      for( int ii=0; ii < 3; ++ii )
      {
        stack.dTractiondw[ i ][ ii ] = m_dTraction_dJump[embSurfIndex][i][ii] * m_surfaceArea[embSurfIndex];
      }
    }
  }

  /**
   * @copydoc geos::finiteElement::ImplicitKernelBase::complete
   */
  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );
    real64 maxForce = 0;
    constexpr int nUdof = numNodesPerElem*3;

    // Compute the local residuals
    LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( stack.localRw, stack.localKww, stack.wLocal );
    LvArray::tensorOps::Ri_add_AijBj< nUdof, 3 >( stack.localRu, stack.localKuw, stack.wLocal );
    // add EqM * effStress into the residual of enrichment nodes
    LvArray::tensorOps::add< 3 >( stack.localRw, stack.localEqMStress);

    // Add traction contribution
    LvArray::tensorOps::scaledAdd< 3 >( stack.localRw, stack.tractionVec, -1 );
    LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, stack.dTractiondw, -1 );

    for( localIndex i = 0; i < nUdof; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRu[i] );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.jumpColIndices,
                                                                              stack.localKuw[i],
                                                                              3 );

    }

    for( localIndex i=0; i < 3; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.jumpEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localRw[i] );

      // fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.jumpColIndices,
                                                                              stack.localKww[i],
                                                                              3 );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              stack.localKwu[i],
                                                                              numNodesPerElem*3 );
    }


    return maxForce;

  }

protected:

  arrayView1d< globalIndex const > const m_wDofNumber;

};

/// The factory used to construct a QuasiStatic kernel.
using EFEMFactory = finiteElement::KernelFactory< EFEM,
                                                  EmbeddedSurfaceSubRegion &,
                                                  arrayView1d< globalIndex const > const,
                                                  arrayView1d< globalIndex const > const,
                                                  globalIndex const,
                                                  CRSMatrixView< real64, globalIndex const > const,
                                                  arrayView1d< real64 > const,
                                                  real64 const,
                                                  real64 const (&) [3] >;
/**
 * @brief A struct to update fracture traction
 */
struct StateUpdateKernel
{

  /**
   * @brief Launch the kernel function doing fracture traction updates
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam FRICTION_WRAPPER the type of contact wrapper doing the fracture traction updates
   * @param[in] size the size of the subregion
   * @param[in] frictionWrapper the wrapper implementing the contact relationship
   * @param[in] jump the displacement jump
   * @param[out] fractureTraction the fracture traction
   * @param[out] dFractureTraction_dJump the derivative of the fracture traction wrt displacement jump
   */
  template< typename POLICY, typename FRICTION_WRAPPER >
  static void
  launch( localIndex const size,
          FRICTION_WRAPPER const & frictionWrapper,
          real64 const contactPenaltyStiffness,
          arrayView2d< real64 const > const & oldJump,
          arrayView2d< real64 const > const & jump,
          arrayView2d< real64 > const & fractureTraction,
          arrayView3d< real64 > const & dFractureTraction_dJump,
          arrayView1d< integer const > const & fractureState,
          arrayView1d< real64 > const & slip )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // Initialize traction and derivatives to 0
      LvArray::forValuesInSlice( fractureTraction[k], []( real64 & val ){ val = 0.0; } );
      LvArray::forValuesInSlice( dFractureTraction_dJump[k], []( real64 & val ){ val = 0.0; } );

      // If the fracture is open the traction is 0 and so are its derivatives so there is nothing to do
      bool const isOpen = fractureState[k] == fields::contact::FractureState::Open;

      if( !isOpen )
      {
        // normal component of the traction
        fractureTraction[k][0] = contactPenaltyStiffness * jump[k][0];

        // derivative of the normal component w.r.t. to the normal dispJump
        dFractureTraction_dJump[k][0][0] = contactPenaltyStiffness;

        frictionWrapper.computeShearTraction( k, oldJump[k], jump[k],
                                              fractureState[k],
                                              fractureTraction[k],
                                              dFractureTraction_dJump[k] );
      }

      slip[ k ] = LvArray::math::sqrt( LvArray::math::square( jump( k, 1 ) ) +
                                       LvArray::math::square( jump( k, 2 ) ) );
    } );
  }

};


} // namespace solidMechanicsEFEMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEFEMKERNELS_HPP_ */
