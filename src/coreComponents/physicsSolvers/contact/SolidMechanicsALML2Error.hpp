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
 * @file SolidMechanicsALML2Error.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALML2ERROR_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALML2ERROR_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

#define TEST_CASE 1

/**
 * @brief l2-error.
 * @copydoc geos::finiteElement::InterfaceKernelBase
 *
 */
template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ALML2Error :
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

  using Base::m_finiteElementSpace;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::InterfaceKernelBase::InterfaceKernelBase
   */
  ALML2Error( NodeManager const & nodeManager,
              EdgeManager const & edgeManager,
              FaceManager const & faceManager,
              localIndex const targetRegionIndex,
              FaceElementSubRegion & elementSubRegion,
              FE_TYPE const & finiteElementSpace,
              CONSTITUTIVE_TYPE & inputConstitutiveType,
              arrayView1d< globalIndex const > const uDofNumber,
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
    m_traction( elementSubRegion.getField< fields::contact::traction >().toViewConst()),
    m_dispJump( elementSubRegion.getField< fields::contact::dispJump >().toViewConst()),
    m_oldDispJump( elementSubRegion.getField< fields::contact::oldDispJump >().toViewConst()),
    m_faceNormal( faceManager.faceNormal() )
  {}

  //***************************************************************************
  /**
   * @copydoc finiteElement::InterfaceKernelBase::StackVariables
   */
  struct StackVariables
  {

    /// The number of lagrange multiplier dofs per element.
    static constexpr int numTdofs = 3;

  public:

    //GEOS_HOST_DEVICE
    StackVariables():
      tLocal{},
      slipLocal{},
      X{ {} }
      //l2ErrLocal
    {}

    /// Stack storage for the element local lagrange multiplier vector
    real64 tLocal[numTdofs];

    /// Stack storage for the element local slip vector
    real64 slipLocal[numTdofs];

    /// local nodal coordinates
    real64 X[ numNodesPerElem ][ 3 ];

    real64 l2ErrLocal;

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
    RAJA::ReduceSum< ReducePolicy< POLICY >, real64 > l2Error( 0 );

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
      l2Error += kernelComponent.complete( k, stack );
    } );

    return static_cast<real64>(l2Error.get());
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
    constexpr int numTdofs = 3;

    int permutation[numNodesPerElem];
    m_finiteElementSpace.getPermutation( permutation );

    localIndex const kf0 = m_elemsToFaces[k][0];
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      for( int i=0; i<3; ++i )
      {
        stack.X[ a ][ i ] = m_X[ m_faceToNodes( kf0, permutation[ a ] ) ][ i ];
      }
    }

    for( int i=0; i<numTdofs; ++i )
    {
      stack.tLocal[i] = m_traction( k, i );
    }

    stack.slipLocal[1] = m_dispJump( k, 1 ) - m_oldDispJump( k, 1 );
    stack.slipLocal[2] = m_dispJump( k, 2 ) - m_oldDispJump( k, 2 );

    stack.l2ErrLocal = 0.0;

  }

  GEOS_HOST_DEVICE
  inline
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );

    real64 const detJ = m_finiteElementSpace.transformedQuadratureWeight( q, stack.X );


    // Test case: verticalFaultOffset                                           
#if TEST_CASE == 0
    real64 N[ numNodesPerElem ];
    m_finiteElementSpace.calcN( q, N );

    real64 y = 0.0;
    for (int a=0; a < numNodesPerElem; ++a )
    {
      y += N[ a ] * stack.X[ a ][ 1 ];
    }

    real64 const slip[2] = { stack.slipLocal[1],
                             stack.slipLocal[2] };
    real64 const slipNorm = LvArray::tensorOps::l2Norm< 2 >( slip );

    real64 const C = ((1 - 2 * m_nu) / (2 * m_pi * (1 - m_nu) ) ) * m_alpha * m_p; 


    //real64 tt = 0.5 * C * std::log((pow((y-m_a),2)*pow((y+m_a),2))/(pow((y-m_b),2)*pow((y+m_b),2)));
    //stack.l2ErrLocal += pow( ( LvArray::math::abs(tt) - LvArray::math::abs( stack.tLocal[2] ) ), 2 ) * detJ;

    real64 const A = m_G / ( 2.0 * m_pi *( 1.0 - m_nu));

    real64 slipA=0.0;
    if (y <= -m_b)
    {
      slipA = 0.0;
    }
    else if ( (y > -m_b) && (y <= -m_a) )
    {
      slipA = -(y + m_b);
    }
    else if ( (y > -m_a) && (y <= m_a) )
    {
      slipA = m_a - m_b;
    }
    else if  ( (y > m_a) && (y <= m_b) )
    {
      slipA = y - m_b;
    }
    else if ( (y >= m_b) )
    {
      slipA = 0.0;
    }
    slipA *= C/A;

    stack.l2ErrLocal += pow( slipNorm - LvArray::math::abs(slipA), 2 ) * detJ;
#endif

    // Test case: singleCrackCompression
#if TEST_CASE == 1
    real64 ym = 0.0;
    for (int a=0; a < numNodesPerElem; ++a )
    {
      ym +=  stack.X[ a ][ 1 ];
    }
    ym /= numNodesPerElem;

    ////real64 tn = m_p * (pow(std::sin(m_alpha*(m_pi/180.0)), 2) + 
    ////                   (m_nu/(1.0 - m_nu)) * pow(std::cos(m_alpha*(m_pi / 180.0 )), 2));

    //real64 tn = m_p * (pow(std::sin(m_alpha*(m_pi/180.0)), 2));
    
    ////if ( ( ym <= 0.60 ) && ( ym >= -0.60 ) )
    //if ( ( ym <= 0.50*1.0 ) && ( ym >= -0.50*1.0 ) )
    //{
    //  stack.l2ErrLocal += pow( tn - stack.tLocal[0], 2 ) * detJ;
    //}

    real64 N[ numNodesPerElem ];
    m_finiteElementSpace.calcN( q, N );

    real64 y = 0.0;
    real64 x = 0.0;
    for (int a=0; a < numNodesPerElem; ++a )
    {
      x += N[ a ] * stack.X[ a ][ 0 ];
      y += N[ a ] * stack.X[ a ][ 1 ];
    }

    real64 const xsiVec[2] = { x, y };
    real64 xsi = LvArray::tensorOps::l2Norm< 2 >( xsiVec );
    xsi = ( y < 0 ) ? -xsi : xsi;
    xsi += m_b;

    real64 const slip[2] = { stack.slipLocal[1],
                             stack.slipLocal[2] };
    real64 const slipNorm = LvArray::tensorOps::l2Norm< 2 >( slip );

    real64 const slipA = 4.0 * ( ( 1.0 - pow(m_nu,2) ) / m_E ) * 
          ( -m_p * std::sin(m_alpha*(m_pi/180.0)) * 
          ( std::cos(m_alpha*(m_pi/180.0)) - std::sin(m_alpha*(m_pi/180.0) ) * std::tan(m_phi*(m_pi/180.0)) ) )   
          * sqrt( pow(m_b,2) - pow((m_b-xsi),2)); 

    stack.l2ErrLocal += pow( slipNorm - slipA, 2 ) * detJ;

#endif

    //stack.l2ErrLocal += pow( ( LvArray::math::abs(tt_m) - LvArray::math::abs( stack.tLocal[2] ) ), 2 ) * detJ;
    //stack.l2ErrLocal += pow( ( tt - tt_m ), 2 ) * detJ;

  }

  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );

    //real64 ym = 0.0;
    //for (int a=0; a < numNodesPerElem; ++a )
    //{
    //  ym +=  stack.X[ a ][ 1 ];
    //}
    //ym /= numNodesPerElem;
    //
    //real64 diam = 0.0;
    //for (int a=0; a < numNodesPerElem; ++a )
    //{
    //  diam +=  LvArray::math::abs(ym - stack.X[ a ][ 1 ] );
    //}
    //diam = 2.0 * diam / numNodesPerElem;

    //return diam*stack.l2ErrLocal;
    //return stack.l2ErrLocal/diam;

    return stack.l2ErrLocal;
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

  arrayView2d< real64 const > const m_traction;

  arrayView2d< real64 const > const m_dispJump;

  arrayView2d< real64 const > const m_oldDispJump;

  arrayView2d< real64 const > const m_faceNormal;

  real64 const m_pi = 3.14159265358979;

#if TEST_CASE == 0
  real64 const m_G = 6500000000;

  real64 const m_nu = 0.15;

  real64 const m_alpha = 0.9;
  
  real64 const m_a = 75;

  real64 const m_b = 150;

  real64 const m_p = -25000000;
#endif

#if TEST_CASE == 1
  real64 const m_E = 1.666666666666666e4;

  real64 const m_nu = 0.25;

  real64 const m_b = 1.0;

  real64 const m_p = -100;

  real64 const m_alpha = 20;

  real64 const m_phi = 30;
#endif


};

/// The factory used to construct the kernel.
using ALML2ErrorFactory = finiteElement::InterfaceKernelFactory< ALML2Error,
                                                          arrayView1d< globalIndex const > const,
                                                          globalIndex const,
                                                          CRSMatrixView< real64, globalIndex const > const,
                                                          arrayView1d< real64 > const,
                                                          real64 const,
                                                          arrayView1d< localIndex const > const >;

#undef TEST_CASE

} // namespace SolidMechanicsALMKernels

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALML2ERROR_HPP_ */
