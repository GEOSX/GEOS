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

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSEFEMKERNELSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSEFEMKERNELSBASE_HPP_

#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"
#include "SolidMechanicsEFEMKernelsHelper.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"

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
class EFEMKernelsBase :
  public solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE,
                                                                             CONSTITUTIVE_TYPE,
                                                                             FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = solidMechanicsLagrangianFEMKernels::ImplicitSmallStrainQuasiStatic< SUBREGION_TYPE,
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
  using Base::m_dt;


  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  EFEMKernelsBase( NodeManager const & nodeManager,
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
          uDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt,
          inputGravityVector ),
    m_w( embeddedSurfSubRegion.getField< fields::contact::dispJump >().toView() ),
    m_tractionVec( embeddedSurfSubRegion.getField< fields::contact::traction >().toViewConst() ),
    m_dTraction_dJump( embeddedSurfSubRegion.getField< fields::contact::dTraction_dJump >().toViewConst() ),
    m_nVec( embeddedSurfSubRegion.getNormalVector().toViewConst() ),
    m_tVec1( embeddedSurfSubRegion.getTangentVector1().toViewConst() ),
    m_tVec2( embeddedSurfSubRegion.getTangentVector2().toViewConst() ),
    m_surfaceCenter( embeddedSurfSubRegion.getElementCenter().toViewConst() ),
    m_surfaceArea( embeddedSurfSubRegion.getElementArea().toViewConst() ),
    m_elementVolume( elementSubRegion.getElementVolume().toViewConst() ),
    m_fracturedElems( elementSubRegion.fracturedElementsList().toViewConst()),
    m_cellsToEmbeddedSurfaces( elementSubRegion.embeddedSurfacesList().toViewConst() )
  {}

  //***************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables  // it's better not to inherit all the stack variable. There is a lot of unused ones.
  {
public:
    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;

    /// The number of jump dofs per element.
    static constexpr int numWdofs = 3;

    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      dispEqnRowIndices{ 0 },
      dispColIndices{ 0 },
      localRu{ 0.0 },
      localRw{ 0.0 },
      localKww{ { 0.0 } },
      localKwu{ { 0.0 } },
      localKuw{ { 0.0 } },
      wLocal(),
      uLocal(),
      hInv(),
      X()
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localRu[numUdofs];

    /// C-array storage for the element local Rw residual vector.
    real64 localRw[numWdofs];

    /// C-array storage for the element local Kww matrix.
    real64 localKww[numWdofs][numWdofs];

    /// C-array storage for the element local Kwu matrix.
    real64 localKwu[numWdofs][numUdofs];

    /// C-array storage for the element local Kuw matrix.
    real64 localKuw[numUdofs][numWdofs];

    /// Stack storage for the element local jump vector
    real64 wLocal[3];

    /// Stack storage for the element displacement vector.
    real64 uLocal[numUdofs];

    /// Stack storage for Area/Volume
    real64 hInv;

    /// local nodal coordinates
    real64 X[ numNodesPerElem ][ 3 ];

    /// Stack storage for the traction
    real64 tractionVec[3];

    /// Stack storage for the derivative of the traction
    real64 dTractiondw[3][3];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[ 6 ][ 6 ];
  };
  //***************************************************************************

  /**
   * @copydoc ::geos::finiteElement::KernelBase::kernelLaunch
   *
   * @detail it uses the kernelLaunch interface of KernelBase but it only launches the kernel
   * on the set of fractured elements within the subregion.
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

    forAll< POLICY >( kernelComponent.m_fracturedElems.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex k = kernelComponent.m_fracturedElems[i];
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
    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.X, dNdX );

    constexpr int nUdof = numNodesPerElem*3;

    // Gauss contribution to Kww, Kwu and Kuw blocks
    real64 Kww_gauss[3][3], Kwu_gauss[3][nUdof], Kuw_gauss[nUdof][3];

    //  Compatibility, equilibrium and strain operators. The compatibility operator is constructed as
    //  a 3 x 6 because it is more convenient for construction purposes (reduces number of local var).
    real64 compMatrix[3][6], strainMatrix[6][nUdof], eqMatrix[3][6];
    real64 matBD[nUdof][6], matED[3][6];

    int Heaviside[ numNodesPerElem ];

    // TODO: asking for the stiffness here will only work for elastic models.  most other models
    //       need to know the strain increment to compute the current stiffness value.

    m_constitutiveUpdate.getElasticStiffness( k, q, stack.constitutiveStiffness );

    solidMechanicsEFEMKernelsHelper::computeHeavisideFunction< numNodesPerElem >( Heaviside,
                                                                                  stack.X,
                                                                                  m_nVec[embSurfIndex],
                                                                                  m_surfaceCenter[embSurfIndex] );


    solidMechanicsEFEMKernelsHelper::assembleEquilibriumOperator( eqMatrix,
                                                                  m_nVec[embSurfIndex],
                                                                  m_tVec1[embSurfIndex],
                                                                  m_tVec2[embSurfIndex],
                                                                  stack.hInv );

    solidMechanicsEFEMKernelsHelper::assembleCompatibilityOperator< numNodesPerElem >( compMatrix,
                                                                                       m_nVec[embSurfIndex],
                                                                                       m_tVec1[embSurfIndex],
                                                                                       m_tVec2[embSurfIndex],
                                                                                       Heaviside,
                                                                                       dNdX );

    solidMechanicsEFEMKernelsHelper::assembleStrainOperator< 6, nUdof, numNodesPerElem >( strainMatrix, dNdX );

    // transp(B)D
    LvArray::tensorOps::Rij_eq_AkiBkj< nUdof, 6, 6 >( matBD, strainMatrix, stack.constitutiveStiffness );
    // ED
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 6, 6 >( matED, eqMatrix, stack.constitutiveStiffness );
    // EDC
    LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 6 >( Kww_gauss, matED, compMatrix );
    // EDB
    LvArray::tensorOps::Rij_eq_AikBkj< 3, nUdof, 6 >( Kwu_gauss, matED, strainMatrix );
    // transp(B)DB
    LvArray::tensorOps::Rij_eq_AikBjk< nUdof, 3, 6 >( Kuw_gauss, matBD, compMatrix );

    // multiply by determinant and add to element matrix
    LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, Kww_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< 3, nUdof >( stack.localKwu, Kwu_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< nUdof, 3 >( stack.localKuw, Kuw_gauss, -detJ );
  }

protected:

  arrayView2d< real64 > const m_w;

  arrayView2d< real64 const > const m_tractionVec;

  arrayView3d< real64 const > const m_dTraction_dJump;

  arrayView2d< real64 const > const m_nVec;

  arrayView2d< real64 const > const m_tVec1;

  arrayView2d< real64 const > const m_tVec2;

  arrayView2d< real64 const > const m_surfaceCenter;

  arrayView1d< real64 const > const m_surfaceArea;

  arrayView1d< real64 const > const m_elementVolume;

  SortedArrayView< localIndex const > const m_fracturedElems;

  ArrayOfArraysView< localIndex const > const m_cellsToEmbeddedSurfaces;
};

} // namespace SolidMechanicsEFEMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSEFEMKERNELSBASE_HPP_ */
