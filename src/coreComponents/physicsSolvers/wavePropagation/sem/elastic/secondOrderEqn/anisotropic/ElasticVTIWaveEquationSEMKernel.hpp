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
 * @file ElasticVTIWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ElasticVTIWaveEquationSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ElasticVTIWaveEquationSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#include "physicsSolvers/wavePropagation/sem/elastic/shared/ElasticFields.hpp"
#include "ElasticVTIFields.hpp"

namespace geos
{

/// Namespace to contain the elastic wave kernels.
namespace elasticVTIWaveEquationSEMKernels
{

template< typename FE_TYPE >
struct ComputeTimeStep
{

  ComputeTimeStep( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Compute timestep using power iteration method
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  real64
  launch( localIndex const sizeElem,
          localIndex const sizeNode,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView1d< real32 const > const density,
          arrayView1d< real32 const > const Vp,
          arrayView1d< real32 const > const Vs,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView1d< real32 > const mass )
  {

    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

    real64 const epsilon = 0.00001;
    localIndex const nIterMax = 10000;
    localIndex numberIter = 0;
    localIndex counter = 0;
    real64 lambdaNew = 0.0;

    array1d< real32 > const ux( sizeNode );
    array1d< real32 > const uy( sizeNode );
    array1d< real32 > const uz( sizeNode );
    array1d< real32 > const uxAux( sizeNode );
    array1d< real32 > const uyAux( sizeNode );
    array1d< real32 > const uzAux( sizeNode );



    arrayView1d< real32 > const uxView = ux;
    arrayView1d< real32 > const uyView = uy;
    arrayView1d< real32 > const uzView = uz;
    arrayView1d< real32 > const uxAuxView = uxAux;
    arrayView1d< real32 > const uyAuxView = uyAux;
    arrayView1d< real32 > const uzAuxView = uzAux;

    //Randomize p values
    srand( time( NULL ));
    for( localIndex a = 0; a < sizeNode; ++a )
    {
      uxView[a] = (real64)rand()/(real64) RAND_MAX;
      uyView[a] = (real64)rand()/(real64) RAND_MAX;
      uzView[a] = (real64)rand()/(real64) RAND_MAX;
    }

    //Step 1: Normalize randomized pressure
    real64 normUx= 0.0;
    real64 normUy= 0.0;
    real64 normUz= 0.0;
    WaveSolverUtils::dotProduct( sizeNode, uxView, uxView, normUx );
    WaveSolverUtils::dotProduct( sizeNode, uyView, uyView, normUy );
    WaveSolverUtils::dotProduct( sizeNode, uzView, uzView, normUz );
    real64 normUtot = normUx+normUy+normUz;


    forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      uxView[a]/= sqrt( normUtot );
      uyView[a]/= sqrt( normUtot );
      uzView[a]/= sqrt( normUtot );
    } );

    //Step 2: Initial iteration of (M^{-1}K)p
    uxAuxView.zero();
    uyAuxView.zero();
    uzAuxView.zero();
    forAll< EXEC_POLICY >( sizeElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      real64 mu = density[k]*Vs[k]*Vs[k];
      real64 lambda = density[k]*Vp[k]*Vp[k]-2*mu;

      real64 xLocal[8][3];
      for( localIndex a=0; a< 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i=0; i<3; ++i )
        {
          xLocal[a][i] = nodeCoords[ nodeIndex ][ i ];
        }
      }
      for( localIndex q = 0; q < numNodesPerElem; ++q )
      {
        m_finiteElement.computeFirstOrderStiffnessTerm( q, xLocal, [&] ( int i, int j, real64 val, real64 J[3][3], int p, int r )
        {
          real32 const Rxx_ij = val*((lambda+2.0*mu)*J[p][0]*J[r][0]+mu*(J[p][1]*J[r][1]+J[p][2]*J[r][2]));
          real32 const Ryy_ij = val*((lambda+2.0*mu)*J[p][1]*J[r][1]+mu*(J[p][0]*J[r][0]+J[p][2]*J[r][2]));
          real32 const Rzz_ij = val*((lambda+2.0*mu)*J[p][2]*J[r][2]+mu*(J[p][0]*J[r][0]+J[p][1]*J[r][1]));
          real32 const Rxy_ij = val*(lambda*J[p][0]*J[r][1]+mu*J[p][1]*J[r][0]);
          real32 const Ryx_ij = val*(mu*J[p][0]*J[r][1]+lambda*J[p][1]*J[r][0]);
          real32 const Rxz_ij = val*(lambda*J[p][0]*J[r][2]+mu*J[p][2]*J[r][0]);
          real32 const Rzx_ij = val*(mu*J[p][0]*J[r][2]+lambda*J[p][2]*J[r][0]);
          real32 const Ryz_ij = val*(lambda*J[p][1]*J[r][2]+mu*J[p][2]*J[r][1]);
          real32 const Rzy_ij = val*(mu*J[p][1]*J[r][2]+lambda*J[p][2]*J[r][1]);

          real32 const localIncrementx = (Rxx_ij * uxView[elemsToNodes[k][j]] + Rxy_ij*uyView[elemsToNodes[k][j]] + Rxz_ij*uzView[elemsToNodes[k][j]]);
          real32 const localIncrementy = (Ryx_ij * uxView[elemsToNodes[k][j]] + Ryy_ij*uyView[elemsToNodes[k][j]] + Ryz_ij*uzView[elemsToNodes[k][j]]);
          real32 const localIncrementz = (Rzx_ij * uxView[elemsToNodes[k][j]] + Rzy_ij*uyView[elemsToNodes[k][j]] + Rzz_ij*uzView[elemsToNodes[k][j]]);

          RAJA::atomicAdd< parallelDeviceAtomic >( &uxAuxView[elemsToNodes[k][i]], localIncrementx );
          RAJA::atomicAdd< parallelDeviceAtomic >( &uyAuxView[elemsToNodes[k][i]], localIncrementy );
          RAJA::atomicAdd< parallelDeviceAtomic >( &uzAuxView[elemsToNodes[k][i]], localIncrementz );
        } );

      }

    } );

    forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      uxAuxView[a]/= mass[a];
      uyAuxView[a]/= mass[a];
      uzAuxView[a]/= mass[a];
    } );
    real64 lambdaOld = lambdaNew;

    //Compute lambdaNew using two dotProducts
    real64 dotProductUxUxaux = 0.0;
    real64 dotProductUyUyaux = 0.0;
    real64 dotProductUzUzaux = 0.0;

    WaveSolverUtils::dotProduct( sizeNode, uxView, uxAuxView, dotProductUxUxaux );
    WaveSolverUtils::dotProduct( sizeNode, uyView, uyAuxView, dotProductUyUyaux );
    WaveSolverUtils::dotProduct( sizeNode, uzView, uzAuxView, dotProductUzUzaux );
    real64 dotProductUtotUtotAux = dotProductUxUxaux+dotProductUyUyaux+dotProductUzUzaux;

    normUx = 0.0;
    normUy = 0.0;
    normUz = 0.0;

    WaveSolverUtils::dotProduct( sizeNode, uxView, uxView, normUx );
    WaveSolverUtils::dotProduct( sizeNode, uyView, uyView, normUy );
    WaveSolverUtils::dotProduct( sizeNode, uzView, uzView, normUz );
    normUtot = normUx+normUy+normUz;


    lambdaNew = dotProductUtotUtotAux/normUtot;

    real64 normUxaux = 0.0;
    real64 normUyaux = 0.0;
    real64 normUzaux = 0.0;
    WaveSolverUtils::dotProduct( sizeNode, uxAuxView, uxAuxView, normUxaux );
    WaveSolverUtils::dotProduct( sizeNode, uyAuxView, uyAuxView, normUyaux );
    WaveSolverUtils::dotProduct( sizeNode, uzAuxView, uzAuxView, normUzaux );

    real64 normUtotAux = normUxaux+normUyaux+normUzaux;


    forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      uxView[a]= uxAuxView[a]/sqrt( normUtotAux );
      uyView[a]= uyAuxView[a]/sqrt( normUtotAux );
      uzView[a]= uzAuxView[a]/sqrt( normUtotAux );
    } );

    //Step 3: Do previous algorithm until we found the max eigenvalues
    do
    {
      uxAuxView.zero();
      uyAuxView.zero();
      uzAuxView.zero();
      forAll< EXEC_POLICY >( sizeElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
      {

        real64 mu = density[k]*Vs[k]*Vs[k];
        real64 lambda = density[k]*Vp[k]*Vp[k]-2*mu;

        real64 xLocal[8][3];
        for( localIndex a=0; a< 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i=0; i<3; ++i )
          {
            xLocal[a][i] = nodeCoords[ nodeIndex ][ i ];
          }
        }
        for( localIndex q = 0; q < numNodesPerElem; ++q )
        {
          m_finiteElement.computeFirstOrderStiffnessTerm( q, xLocal, [&] ( int i, int j, real64 val, real64 J[3][3], int p, int r )
          {
            real32 const Rxx_ij = val*((lambda+2.0*mu)*J[p][0]*J[r][0]+mu*(J[p][1]*J[r][1]+J[p][2]*J[r][2]));
            real32 const Ryy_ij = val*((lambda+2.0*mu)*J[p][1]*J[r][1]+mu*(J[p][0]*J[r][0]+J[p][2]*J[r][2]));
            real32 const Rzz_ij = val*((lambda+2.0*mu)*J[p][2]*J[r][2]+mu*(J[p][0]*J[r][0]+J[p][1]*J[r][1]));
            real32 const Rxy_ij = val*(lambda*J[p][0]*J[r][1]+mu*J[p][1]*J[r][0]);
            real32 const Ryx_ij = val*(mu*J[p][0]*J[r][1]+lambda*J[p][1]*J[r][0]);
            real32 const Rxz_ij = val*(lambda*J[p][0]*J[r][2]+mu*J[p][2]*J[r][0]);
            real32 const Rzx_ij = val*(mu*J[p][0]*J[r][2]+lambda*J[p][2]*J[r][0]);
            real32 const Ryz_ij = val*(lambda*J[p][1]*J[r][2]+mu*J[p][2]*J[r][1]);
            real32 const Rzy_ij = val*(mu*J[p][1]*J[r][2]+lambda*J[p][2]*J[r][1]);

            real32 const localIncrementx = (Rxx_ij * uxView[elemsToNodes[k][j]] + Rxy_ij*uyView[elemsToNodes[k][j]] + Rxz_ij*uzView[elemsToNodes[k][j]]);
            real32 const localIncrementy = (Ryx_ij * uxView[elemsToNodes[k][j]] + Ryy_ij*uyView[elemsToNodes[k][j]] + Ryz_ij*uzView[elemsToNodes[k][j]]);
            real32 const localIncrementz = (Rzx_ij * uxView[elemsToNodes[k][j]] + Rzy_ij*uyView[elemsToNodes[k][j]] + Rzz_ij*uzView[elemsToNodes[k][j]]);

            RAJA::atomicAdd< parallelDeviceAtomic >( &uxAuxView[elemsToNodes[k][i]], localIncrementx );
            RAJA::atomicAdd< parallelDeviceAtomic >( &uyAuxView[elemsToNodes[k][i]], localIncrementy );
            RAJA::atomicAdd< parallelDeviceAtomic >( &uzAuxView[elemsToNodes[k][i]], localIncrementz );
          } );

        }

      } );

      forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        uxAuxView[a]/= mass[a];
        uyAuxView[a]/= mass[a];
        uzAuxView[a]/= mass[a];
      } );
      lambdaOld = lambdaNew;

      //Compute lambdaNew using two dotProducts
      dotProductUxUxaux = 0.0;
      dotProductUyUyaux = 0.0;
      dotProductUzUzaux = 0.0;

      WaveSolverUtils::dotProduct( sizeNode, uxView, uxAuxView, dotProductUxUxaux );
      WaveSolverUtils::dotProduct( sizeNode, uyView, uyAuxView, dotProductUyUyaux );
      WaveSolverUtils::dotProduct( sizeNode, uzView, uzAuxView, dotProductUzUzaux );
      dotProductUtotUtotAux = dotProductUxUxaux+dotProductUyUyaux+dotProductUzUzaux;

      normUx = 0.0;
      normUy = 0.0;
      normUz = 0.0;

      WaveSolverUtils::dotProduct( sizeNode, uxView, uxView, normUx );
      WaveSolverUtils::dotProduct( sizeNode, uyView, uyView, normUy );
      WaveSolverUtils::dotProduct( sizeNode, uzView, uzView, normUz );
      normUtot = normUx+normUy+normUz;


      lambdaNew = dotProductUtotUtotAux/normUtot;

      normUxaux = 0.0;
      normUyaux = 0.0;
      normUzaux = 0.0;
      WaveSolverUtils::dotProduct( sizeNode, uxAuxView, uxAuxView, normUxaux );
      WaveSolverUtils::dotProduct( sizeNode, uyAuxView, uyAuxView, normUyaux );
      WaveSolverUtils::dotProduct( sizeNode, uzAuxView, uzAuxView, normUzaux );

      normUtotAux = normUxaux+normUyaux+normUzaux;


      forAll< EXEC_POLICY >( sizeNode, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        uxView[a]= uxAuxView[a]/sqrt( normUtotAux );
        uyView[a]= uyAuxView[a]/sqrt( normUtotAux );
        uzView[a]= uzAuxView[a]/sqrt( normUtotAux );
      } );

      if( LvArray::math::abs( lambdaNew-lambdaOld )/LvArray::math::abs( lambdaNew )<= epsilon )
      {
        counter++;
      }
      else
      {
        counter=0;
      }

      numberIter++;


    }
    while (counter < 10 && numberIter < nIterMax);

    GEOS_THROW_IF( numberIter> nIterMax, "Power Iteration algorithm does not converge", std::runtime_error );

    real64 dt = 1.99/sqrt( LvArray::math::abs( lambdaNew ));

    return dt;

  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;
};

/**
 * @brief Implements kernels for solving the elastic wave equations
 *   explicit central FD method and SEM in the Vertical Transverse Isotropic (VTI) case
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### ElasticVTIWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the acoustic wave equations using the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitElasticVTISEM : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                                CONSTITUTIVE_TYPE,
                                                                FE_TYPE,
                                                                1,
                                                                1 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

//*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param dt The time interval for the step.
   */
  ExplicitElasticVTISEM( NodeManager & nodeManager,
                         EdgeManager const & edgeManager,
                         FaceManager const & faceManager,
                         localIndex const targetRegionIndex,
                         SUBREGION_TYPE const & elementSubRegion,
                         FE_TYPE const & finiteElementSpace,
                         CONSTITUTIVE_TYPE & inputConstitutiveType,
                         real64 const dt ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
    m_nodeCoords( nodeManager.getField< fields::referencePosition32 >() ),
    m_ux_n( nodeManager.getField< fields::elasticfields::Displacementx_n >() ),
    m_uy_n( nodeManager.getField< fields::elasticfields::Displacementy_n >() ),
    m_uz_n( nodeManager.getField< fields::elasticfields::Displacementz_n >() ),
    m_stiffnessVectorx( nodeManager.getField< fields::elasticfields::StiffnessVectorx >() ),
    m_stiffnessVectory( nodeManager.getField< fields::elasticfields::StiffnessVectory >() ),
    m_stiffnessVectorz( nodeManager.getField< fields::elasticfields::StiffnessVectorz >() ),
    m_density( elementSubRegion.template getField< fields::elasticfields::ElasticDensity >() ),
    m_velocityVp( elementSubRegion.template getField< fields::elasticfields::ElasticVelocityVp >() ),
    m_velocityVs( elementSubRegion.template getField< fields::elasticfields::ElasticVelocityVs >() ),
    m_gamma( elementSubRegion.template getField< fields::elasticvtifields::Gamma >()),
    m_epsilon( elementSubRegion.template getField< fields::elasticvtifields::Epsilon >()),
    m_delta( elementSubRegion.template getField< fields::elasticvtifields::Delta >()),
    m_dt( dt )
  {
    GEOS_UNUSED_VAR( edgeManager );
    GEOS_UNUSED_VAR( faceManager );
    GEOS_UNUSED_VAR( targetRegionIndex );
  }



  //*****************************************************************************
  /**
   * @copydoc geos::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitElasticVTISEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal(),
      stiffnessVectorxLocal(),
      stiffnessVectoryLocal(),
      stiffnessVectorzLocal()
    {}
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ 8 ][ 3 ]{};
    real32 stiffnessVectorxLocal[ numNodesPerElem ];
    real32 stiffnessVectoryLocal[ numNodesPerElem ];
    real32 stiffnessVectorzLocal[ numNodesPerElem ];
    real32 Cvti[6];

  };
  //***************************************************************************


  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a< 8; a++ )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_nodeCoords[ nodeIndex ][ i ];
      }
    }

    stack.Cvti[0] = m_density[k] * pow( m_velocityVp[k], 2 ) * (1.0 + 2.0*m_epsilon[k]);
    stack.Cvti[1] = m_density[k] * pow( m_velocityVp[k], 2 );
    stack.Cvti[2] = m_density[k] *
                    sqrt( pow((pow( m_velocityVp[k],
                                    2 ) - pow( m_velocityVs[k], 2 )),
                              2 ) + 2.0 * pow( m_velocityVp[k], 2 ) * m_delta[k] * (pow( m_velocityVp[k], 2 ) - pow( m_velocityVs[k], 2 )) ) - m_density[k] * pow(
      m_velocityVs[k], 2 );
    stack.Cvti[3] = m_density[k] * pow( m_velocityVs[k], 2 );
    stack.Cvti[4] = m_density[k] * pow( m_velocityVs[k], 2 )*(1.0 + 2.0 * m_gamma[k]);
    stack.Cvti[5] = stack.Cvti[0] - 2.0 * stack.Cvti[4];

  }

  /**
   * @copydoc geos::finiteElement::KernelBase::complete
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    for( int i=0; i<numNodesPerElem; i++ )
    {
      const localIndex nodeIndex = m_elemsToNodes( k, i );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectorx[ nodeIndex ], stack.stiffnessVectorxLocal[ i ] );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectory[ nodeIndex ], stack.stiffnessVectoryLocal[ i ] );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVectorz[ nodeIndex ], stack.stiffnessVectorzLocal[ i ] );
    }
    return 0;
  }

  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitElasticVTISEM Description
   * Calculates stiffness vector
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    m_finiteElementSpace.template computeFirstOrderStiffnessTerm( q, stack.xLocal, [&] ( int i, int j, real64 val, real64 J[3][3], int p, int r )
    {
      real32 const Rxx_ij = val*(stack.Cvti[0]*J[p][0]*J[r][0]+stack.Cvti[4]*(J[p][1]*J[r][1])+stack.Cvti[3]*(J[p][2]*J[r][2]));
      real32 const Ryy_ij = val*(stack.Cvti[0]*J[p][1]*J[r][1]+stack.Cvti[4]*(J[p][0]*J[r][0])+stack.Cvti[3]*(J[p][2]*J[r][2]));
      real32 const Rzz_ij = val*(stack.Cvti[1]*J[p][2]*J[r][2]+stack.Cvti[3]*(J[p][0]*J[r][0]+J[p][1]*J[r][1]));
      real32 const Rxy_ij = val*(stack.Cvti[5]*J[p][0]*J[r][1]+stack.Cvti[4]*J[p][1]*J[r][0]);
      real32 const Ryx_ij = val*(stack.Cvti[4]*J[p][0]*J[r][1]+stack.Cvti[5]*J[p][1]*J[r][0]);
      real32 const Rxz_ij = val*(stack.Cvti[2]*J[p][0]*J[r][2]+stack.Cvti[3]*J[p][2]*J[r][0]);
      real32 const Rzx_ij = val*(stack.Cvti[3]*J[p][0]*J[r][2]+stack.Cvti[2]*J[p][2]*J[r][0]);
      real32 const Ryz_ij = val*(stack.Cvti[2]*J[p][1]*J[r][2]+stack.Cvti[3]*J[p][2]*J[r][1]);
      real32 const Rzy_ij = val*(stack.Cvti[3]*J[p][1]*J[r][2]+stack.Cvti[2]*J[p][2]*J[r][1]);

      real32 const localIncrementx = (Rxx_ij * m_ux_n[m_elemsToNodes( k, j )] + Rxy_ij*m_uy_n[m_elemsToNodes( k, j )] + Rxz_ij*m_uz_n[m_elemsToNodes( k, j )]);
      real32 const localIncrementy = (Ryx_ij * m_ux_n[m_elemsToNodes( k, j )] + Ryy_ij*m_uy_n[m_elemsToNodes( k, j )] + Ryz_ij*m_uz_n[m_elemsToNodes( k, j )]);
      real32 const localIncrementz = (Rzx_ij * m_ux_n[m_elemsToNodes( k, j )] + Rzy_ij*m_uy_n[m_elemsToNodes( k, j )] + Rzz_ij*m_uz_n[m_elemsToNodes( k, j )]);

      stack.stiffnessVectorxLocal[ i ] += localIncrementx;
      stack.stiffnessVectoryLocal[ i ] += localIncrementy;
      stack.stiffnessVectorzLocal[ i ] += localIncrementz;
    } );
  }


protected:
  /// The array containing the nodal position array.
  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const m_nodeCoords;

  /// The array containing the nodal displacement array in x direction.
  arrayView1d< real32 > const m_ux_n;

  /// The array containing the nodal displacement array in y direction.
  arrayView1d< real32 > const m_uy_n;

  /// The array containing the nodal displacement array in z direction.
  arrayView1d< real32 > const m_uz_n;

  /// The array containing the product of the stiffness matrix and the nodal displacement.
  arrayView1d< real32 > const m_stiffnessVectorx;

  /// The array containing the product of the stiffness matrix and the nodal displacement.
  arrayView1d< real32 > const m_stiffnessVectory;

  /// The array containing the product of the stiffness matrix and the nodal displacement.
  arrayView1d< real32 > const m_stiffnessVectorz;

  /// The array containing the density of the medium
  arrayView1d< real32 const > const m_density;

  /// The array containing the P-wavespeed
  arrayView1d< real32 const > const m_velocityVp;

  /// The array containing the S-wavespeed
  arrayView1d< real32 const > const m_velocityVs;

  ///The array containing the Thomsen constant gamma
  arrayView1d< real32 const > const m_gamma;

  ///The array containing the Thomsen constant epsilon
  arrayView1d< real32 const > const m_epsilon;

  ///The array containing the Thomsen constant delta
  arrayView1d< real32 const > const m_delta;

  /// The time increment for this time integration step.
  real64 const m_dt;


};


/// The factory used to construct a ExplicitAcousticWaveEquation kernel.
using ExplicitElasticVTISEMFactory = finiteElement::KernelFactory< ExplicitElasticVTISEM,
                                                                   real64 >;

} // namespace elasticVTIWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ElasticVTIWaveEquationSEMKERNEL_HPP_
