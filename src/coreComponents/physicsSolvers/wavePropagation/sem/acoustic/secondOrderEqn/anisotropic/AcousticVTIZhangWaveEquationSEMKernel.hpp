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
 * @file AcousticVTIZhangWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIZHANGWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIZHANGWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticFields.hpp"
#include "AcousticVTIFields.hpp"

namespace geos
{

/// Namespace to contain the acoustic wave kernels.
namespace acousticVTIZhangWaveEquationSEMKernels
{


/**
 * @brief Implements kernels for solving the pseudo-acoustic VTI wave equations
 *   explicit central FD method and SEM
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticVTIZhangWaveEquationSEMKernel Description
 * Implements the KernelBase interface functions required for solving
 * the VTI pseudo-acoustic wave Zhang's set of equations using the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1`.
 */


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitAcousticVTIZhangSEM : public finiteElement::KernelBase< SUBREGION_TYPE,
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

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
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
   *   elements to be processed during this kernel launch.
   */
  ExplicitAcousticVTIZhangSEM( NodeManager & nodeManager,
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
    m_p_n( nodeManager.getField< geos::fields::acousticvtifields::Pressure_p_n >() ),
    m_q_n( nodeManager.getField< geos::fields::acousticvtifields::Pressure_q_n >() ),
    m_stiffnessVector_p( nodeManager.getField< geos::fields::acousticvtifields::StiffnessVector_p >() ),
    m_stiffnessVector_q( nodeManager.getField< geos::fields::acousticvtifields::StiffnessVector_q >() ),
    m_density( elementSubRegion.template getField< geos::fields::acousticfields::AcousticDensity >() ),
    m_vti_epsilon( elementSubRegion.template getField< geos::fields::acousticvtifields::AcousticEpsilon >() ),
    m_vti_delta( elementSubRegion.template getField< geos::fields::acousticvtifields::AcousticDelta >() ),
    m_vti_DofEpsilon( nodeManager.getField< geos::fields::acousticvtifields::AcousticDofEpsilon >() ),
    m_vti_DofDelta( nodeManager.getField< geos::fields::acousticvtifields::AcousticDofDelta >() ),
    m_vti_GradzDelta( elementSubRegion.template getField< geos::fields::acousticvtifields::AcousticGradzDelta >() ),
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
   * ### ExplicitAcousticSEM Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
public:
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal(),
      stiffnessVectorLocal_p(),
      stiffnessVectorLocal_q()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ 8 ][ 3 ];
    real32 stiffnessVectorLocal_p[ numNodesPerElem ]{};
    real32 stiffnessVectorLocal_q[ numNodesPerElem ]{};
    real32 invDensity;
    //Debug
#if 1
    real32 vti_epsi;     // (1 + 2*epsilon)
    real32 vti_sqrtDelta;     // sqrt(1 + 2*delta)
    real32 vti_GradzDelta;
#else
#endif
    //End Debug
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
#if 1
    real32 epsi = std::fabs( m_vti_epsilon[k] );
    real32 delt = std::fabs( m_vti_delta[k] );
    if( std::fabs( epsi ) < 1e-5 )
      epsi = 0;
    if( std::fabs( delt ) < 1e-5 )
      delt = 0;
    if( delt > epsi )
      delt = epsi;
    stack.vti_epsi = (1 + 2 * epsi);
    stack.vti_sqrtDelta = sqrt( 1 + 2 * delt );
    stack.vti_GradzDelta = m_vti_GradzDelta[k];
#endif

    stack.invDensity = 1./m_density[k];
    for( localIndex a=0; a< 8; a++ )
    {
      localIndex const nodeIndex =  m_elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_nodeCoords[ nodeIndex ][ i ];
      }
    }
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
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_p[m_elemsToNodes( k, i )], stack.stiffnessVectorLocal_p[i] );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_q[m_elemsToNodes( k, i )], stack.stiffnessVectorLocal_q[i] );
    }
    return 0;
  }


  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticSEM Description
   * Calculates stiffness vector
   *
   */

#if 0
  // Original version with delta and epsilon constant per element.
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Pseudo Stiffness xy
    m_finiteElementSpace.template computeStiffnessxyTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {

      real32 const localIncrement_p = -val * stack.invDensity * stack.vti_epsi * m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[i] += localIncrement_p;
      real32 const localIncrement_q = -val * stack.invDensity * stack.vti_sqrtDelta * m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_q[i] += localIncrement_q;
    } );

    // Pseudo-Stiffness z
    m_finiteElementSpace.template computeStiffnesszTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 const localIncrement_p = -val * stack.invDensity * stack.vti_sqrtDelta* m_q_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[i] += localIncrement_p;

      real32 const localIncrement_q = -val * stack.invDensity * m_q_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_q[i] += localIncrement_q;
    } );
  }

#else
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Pseudo Stiffness xy
    m_finiteElementSpace.template computeStiffnessxyTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 epsi = std::fabs( m_vti_DofEpsilon[m_elemsToNodes( k, q )] ); // value on control point
      real32 delt = std::fabs( m_vti_DofDelta[m_elemsToNodes( k, q )] ); // value on control point

      if( std::fabs( epsi ) < 1e-5 )
        epsi = 0;
      if( std::fabs( delt ) < 1e-5 )
        delt = 0;
      if( delt > epsi )
        delt = epsi;
      real32 vti_sqrtDelta = std::sqrt( 1 + 2 *delt );

      real32 const localIncrement_p = -val * stack.invDensity * (1 + 2 * epsi) * m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[i] += localIncrement_p;
      real32 const localIncrement_q = -val * stack.invDensity * vti_sqrtDelta * m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_q[i] += localIncrement_q;
    } );

    // Pseudo-Stiffness z
    m_finiteElementSpace.template computeStiffnesszTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 epsi = std::fabs( m_vti_DofEpsilon[m_elemsToNodes( k, q )] ); // value on control point
      real32 delt = std::fabs( m_vti_DofDelta[m_elemsToNodes( k, q )] ); // value on control point

      if( std::fabs( epsi ) < 1e-5 )
        epsi = 0;
      if( std::fabs( delt ) < 1e-5 )
        delt = 0;
      if( delt > epsi )
        delt = epsi;
      real32 vti_sqrtDelta = sqrt( 1 + 2 *delt );

      real32 const localIncrement_p = -val * stack.invDensity * vti_sqrtDelta* m_q_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[i] += localIncrement_p;

      real32 const localIncrement_q = -val * stack.invDensity * m_q_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_q[i] += localIncrement_q;
    } );

#if 0
    // Missing dz term with precomputed dz(f)
    m_finiteElementSpace.template computeMissingzVolumeTerm_precompDzF( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 GradzsqrtDelta = stack.vti_GradzDelta / stack.vti_sqrtDelta;

      real32 const localIncrement_p = -val * stack.invDensity * GradzsqrtDelta* m_q_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[i] += localIncrement_p;
    } );

#else
    // Missing dz term
    m_finiteElementSpace.template computeMissingzVolumeTerm( q, stack.xLocal, [&] ( int iVertice, int j, real64 val )
    {
      //iVertice is the "Qr" index of the vertice in the element (so 0 < iVertice < (r+1)^3)
      real32 epsi = std::fabs( m_vti_DofEpsilon[m_elemsToNodes( k, iVertice )] ); // value on k
      real32 delt = std::fabs( m_vti_DofDelta[m_elemsToNodes( k, iVertice )] ); // value on
      if( std::fabs( epsi ) < 1e-5 )
        epsi = 0;
      if( std::fabs( delt ) < 1e-5 )
        delt = 0;
      if( delt > epsi )
        delt = epsi;

      // Two options for defining "f"
      real32 vti_sqrtDelta = delt / sqrt( 1 + 2 *delt );
      //real32 vti_sqrtDelta = sqrt( 1 + 2 *delt );

      real32 const localIncrement_p = -val * stack.invDensity * vti_sqrtDelta * m_q_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[q] += localIncrement_p;
    } );
#endif


#if 1
    // Missing dxy term
    m_finiteElementSpace.template computeMissingxyVolumeTerm( q, stack.xLocal, [&] ( int iVertice, int j, real64 val )
    {
      //iVertice is the "Qr" index of the vertice in the element (so 0 < iVertice < (r+1)^3)
      real32 epsi = std::fabs( m_vti_DofEpsilon[m_elemsToNodes( k, iVertice )] ); // value on k
      real32 delt = std::fabs( m_vti_DofDelta[m_elemsToNodes( k, iVertice )] ); // value on
      if( std::fabs( epsi ) < 1e-5 )
        epsi = 0;
      if( std::fabs( delt ) < 1e-5 )
        delt = 0;
      if( delt > epsi )
        delt = epsi;



      real32 const localIncrement_p = -val * stack.invDensity * (1 + 2 * epsi) * m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_p[q] += localIncrement_p;



      // Two options for defining "f"
      real32 vti_sqrtDelta = delt / sqrt( 1 + 2 *delt );
      //real32 vti_sqrtDelta = sqrt( 1 + 2 *delt );


      real32 const localIncrement_q = -val * stack.invDensity * vti_sqrtDelta * m_p_n[m_elemsToNodes( k, j )];
      stack.stiffnessVectorLocal_q[q] += localIncrement_q;
    } );



#endif

  }
#endif

protected:
  /// The array containing the nodal position array.
  arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const m_nodeCoords;

  /// The array containing the nodal pressure array.
  arrayView1d< real32 const > const m_p_n;

  /// The array containing the nodal auxiliary variable array.
  arrayView1d< real32 const > const m_q_n;

  /// The array containing the product of the stiffness matrix and the nodal pressure for the equation in p.
  arrayView1d< real32 > const m_stiffnessVector_p;

  /// The array containing the product of the stiffness matrix and the nodal pressure for the equation in q.
  arrayView1d< real32 > const m_stiffnessVector_q;

  /// The array containing the medium density.
  arrayView1d< real32 const > const m_density;


  /// The array containing the epsilon Thomsen parameter.
  arrayView1d< real32 const > const m_vti_epsilon;
  /// The array containing the delta Thomsen parameter.
  arrayView1d< real32 const > const m_vti_delta;

  /// The array containing the epsilon Thomsen parameter.
  arrayView1d< real32 const > const m_vti_DofEpsilon;
  /// The array containing the delta Thomsen parameter.
  arrayView1d< real32 const > const m_vti_DofDelta;



  /// dz delta.
  arrayView1d< real32 const > const m_vti_GradzDelta;

  /// The time increment for this time integration step.
  real64 const m_dt;
};



/// The factory used to construct a ExplicitAcousticVTIZhang kernel.
using ExplicitAcousticVTIZhangSEMFactory = finiteElement::KernelFactory< ExplicitAcousticVTIZhangSEM,
                                                                         real64 >;


} // namespace acousticVTIZhangWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIZHANGWAVEEQUATIONSEMKERNEL_HPP_
