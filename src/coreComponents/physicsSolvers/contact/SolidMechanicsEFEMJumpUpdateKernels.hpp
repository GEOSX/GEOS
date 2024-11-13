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

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEFEMJUMPUPDATEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEFEMJUMPUPDATEKERNELS_HPP_

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
class EFEMJumpUpdate :
  public EFEMKernelsBase< SUBREGION_TYPE,
                          CONSTITUTIVE_TYPE,
                          FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = EFEMKernelsBase< SUBREGION_TYPE,
                                CONSTITUTIVE_TYPE,
                                FE_TYPE >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
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
  using Base::m_w;
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
  EFEMJumpUpdate( NodeManager const & nodeManager,
                  EdgeManager const & edgeManager,
                  FaceManager const & faceManager,
                  localIndex const targetRegionIndex,
                  SUBREGION_TYPE const & elementSubRegion,
                  FE_TYPE const & finiteElementSpace,
                  CONSTITUTIVE_TYPE & inputConstitutiveType,
                  EmbeddedSurfaceSubRegion & embeddedSurfSubRegion,
                  arrayView1d< globalIndex const > const uDofNumber,
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
          inputGravityVector )
  {}

  //***************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;

    /// Stack storage for the element displacement vector.
    real64 dUlocal[numUdofs];
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
        stack.dUlocal[ a*3 + i ] = m_uhat[localNodeIndex][i];
      }
    }

    for( int i=0; i<3; ++i )
    {
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
    // LvArray::tensorOps::Ri_add_AijBj< 3, nUdof >( stack.localRw, stack.localKwu, stack.uLocal );

    // Add traction contribution tranction
    LvArray::tensorOps::scaledAdd< 3 >( stack.localRw, stack.tractionVec, -1 );
    LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, stack.dTractiondw, -1 );

    // Compute dw given du
    real64 dWlocal[3], InvKww[3][3], InvKwwKwu[3][nUdof];

    // compute inv(Kww)
    LvArray::tensorOps::invert< 3 >( InvKww, stack.localKww );
    // compute inv(kww)Kwu
    LvArray::tensorOps::Rij_eq_AikBkj< 3, nUdof, 3 >( InvKwwKwu, InvKww, stack.localKwu );
    // add invKww Rw to dw
    LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( dWlocal, InvKww, stack.localRw );
    // add inv(kww)Kwu * du to dw
    LvArray::tensorOps::Ri_add_AijBj< 3, nUdof >( dWlocal, InvKwwKwu, stack.dUlocal );
    LvArray::tensorOps::scale< 3 >( dWlocal, -1 );

    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

    for( int i=0; i<3; ++i )
    {
      m_w[ embSurfIndex ][i] = stack.wLocal[ i ] + dWlocal[i];
    }

    return maxForce;

  }

};

/// The factory used to construct a QuasiStatic kernel.
using EFEMJumpUpdateFactory = finiteElement::KernelFactory< EFEMJumpUpdate,
                                                            EmbeddedSurfaceSubRegion &,
                                                            arrayView1d< globalIndex const > const,
                                                            globalIndex const,
                                                            CRSMatrixView< real64, globalIndex const > const,
                                                            arrayView1d< real64 > const,
                                                            real64 const,
                                                            real64 const (&) [3] >;

} // namespace SolidMechanicsEFEMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSEFEMJUMPUPDATEKERNELS_HPP_ */
