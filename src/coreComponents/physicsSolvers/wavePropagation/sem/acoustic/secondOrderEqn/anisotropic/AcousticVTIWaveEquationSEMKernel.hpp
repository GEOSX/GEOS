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
 * @file AcousticVTIWaveEquationSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIWAVEEQUATIONSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIWAVEEQUATIONSEMKERNEL_HPP_

#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#include "AcousticVTIFields.hpp"

namespace geos
{

/// Namespace to contain the acoustic wave kernels.
namespace acousticVTIWaveEquationSEMKernels
{

struct PrecomputeSourceAndReceiverKernel
{
  /**
   * @brief Launches the precomputation of the source and receiver terms
   * @tparam EXEC_POLICY execution policy
   * @tparam FE_TYPE finite element type
   * @param[in] size the number of cells in the subRegion
   * @param[in] baseFacesToNodes face to node map
   * @param[in] baseNodeCoords coordinates of the nodes
   * @param[in] baseNodeLocalToGlobal local to global index map for nodes
   * @param[in] elementLocalToGlobal local to global index map for elements
   * @param[in] baseNodesToElements node to element map for the base mesh
   * @param[in] baseElemsToNodes element to node map for the base mesh
   * @param[in] elemGhostRank rank of the ghost element
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] elemsToFaces map from element to faces
   * @param[in] elemCenter coordinates of the element centers
   * @param[in] sourceCoordinates coordinates of the source terms
   * @param[out] sourceIsAccessible flag indicating whether the source is accessible or not
   * @param[out] sourceNodeIds indices of the nodes of the element where the source is located
   * @param[out] sourceConstants constant part of the source terms
   * @param[in] receiverCoordinates coordinates of the receiver terms
   * @param[out] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[out] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[out] receiverConstants constant part of the receiver term
   * @param[out] sourceValue value of the temporal source (eg. Ricker)
   * @param[in] dt time-step
   * @param[in] timeSourceFrequency the central frequency of the source
   * @param[in] timeSourceDelay the time delay of the source
   * @param[in] rickerOrder order of the Ricker wavelet
   */
  template< typename EXEC_POLICY, typename FE_TYPE >
  static void
  launch( localIndex const size,
          ArrayOfArraysView< localIndex const > const baseFacesToNodes,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const baseNodeCoords,
          arrayView1d< globalIndex const > const baseNodeLocalToGlobal,
          arrayView1d< globalIndex const > const elementLocalToGlobal,
          ArrayOfArraysView< localIndex const > const baseNodesToElements,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & baseElemsToNodes,
          arrayView1d< integer const > const elemGhostRank,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes,
          arrayView2d< localIndex const > const elemsToFaces,
          arrayView2d< real64 const > const & elemCenter,
          arrayView2d< real64 const > const sourceCoordinates,
          arrayView1d< localIndex > const sourceIsAccessible,
          arrayView2d< localIndex > const sourceNodeIds,
          arrayView2d< real64 > const sourceConstants,
          arrayView2d< real64 const > const receiverCoordinates,
          arrayView1d< localIndex > const receiverIsLocal,
          arrayView2d< localIndex > const receiverNodeIds,
          arrayView2d< real64 > const receiverConstants,
          arrayView2d< real32 > const sourceValue,
          real64 const dt,
          real32 const timeSourceFrequency,
          real32 const timeSourceDelay,
          localIndex const rickerOrder )
  {
    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      real64 const center[3] = { elemCenter[k][0],
                                 elemCenter[k][1],
                                 elemCenter[k][2] };

      // Step 1: locate the sources, and precompute the source term

      /// loop over all the source that haven't been found yet
      for( localIndex isrc = 0; isrc < sourceCoordinates.size( 0 ); ++isrc )
      {
        if( sourceIsAccessible[isrc] == 0 )
        {
          real64 const coords[3] = { sourceCoordinates[isrc][0],
                                     sourceCoordinates[isrc][1],
                                     sourceCoordinates[isrc][2] };

          bool const sourceFound =
            computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                        baseNodeCoords,
                                                                        elemsToFaces,
                                                                        baseFacesToNodes,
                                                                        baseNodesToElements,
                                                                        baseNodeLocalToGlobal,
                                                                        elementLocalToGlobal,
                                                                        center,
                                                                        coords );
          if( sourceFound )
          {
            real64 coordsOnRefElem[3]{};


            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              baseElemsToNodes[k],
                                                                              baseNodeCoords,
                                                                              coordsOnRefElem );

            sourceIsAccessible[isrc] = 1;
            real64 Ntest[numNodesPerElem];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              sourceNodeIds[isrc][a] = elemsToNodes[k][a];
              sourceConstants[isrc][a] = Ntest[a];
            }

            for( localIndex cycle = 0; cycle < sourceValue.size( 0 ); ++cycle )
            {
              sourceValue[cycle][isrc] = WaveSolverUtils::evaluateRicker( cycle * dt, timeSourceFrequency, timeSourceDelay, rickerOrder );
            }
          }
        }
      } // end loop over all sources


      // Step 2: locate the receivers, and precompute the receiver term

      /// loop over all the receivers that haven't been found yet
      for( localIndex ircv = 0; ircv < receiverCoordinates.size( 0 ); ++ircv )
      {
        if( receiverIsLocal[ircv] == 0 )
        {
          real64 const coords[3] = { receiverCoordinates[ircv][0],
                                     receiverCoordinates[ircv][1],
                                     receiverCoordinates[ircv][2] };

          real64 coordsOnRefElem[3]{};

          bool const receiverFound =
            computationalGeometry::isPointInsideConvexPolyhedronRobust( k,
                                                                        baseNodeCoords,
                                                                        elemsToFaces,
                                                                        baseFacesToNodes,
                                                                        baseNodesToElements,
                                                                        baseNodeLocalToGlobal,
                                                                        elementLocalToGlobal,
                                                                        center,
                                                                        coords );
          if( receiverFound && elemGhostRank[k] < 0 )
          {
            WaveSolverUtils::computeCoordinatesOnReferenceElement< FE_TYPE >( coords,
                                                                              baseElemsToNodes[k],
                                                                              baseNodeCoords,
                                                                              coordsOnRefElem );

            receiverIsLocal[ircv] = 1;

            real64 Ntest[numNodesPerElem];
            FE_TYPE::calcN( coordsOnRefElem, Ntest );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              receiverNodeIds[ircv][a] = elemsToNodes[k][a];
              receiverConstants[ircv][a] = Ntest[a];
            }
          }
        }
      } // end loop over receivers

    } );

  }
};

template< typename FE_TYPE >
struct MassMatrixKernel
{

  MassMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the mass matrices
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToNodes map from element to nodes
   * @param[in] velocity cell-wise velocity
   * @param[out] mass diagonal of the mass matrix
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
          arrayView1d< real32 const > const velocity,
          arrayView1d< real32 > const mass )

  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
    {


      real32 const invC2 = 1.0 / pow( velocity[e], 2 );
      // only the eight corners of the mesh cell are needed to compute the Jacobian
      real64 xLocal[ 8 ][ 3 ];
      for( localIndex a = 0; a < 8; ++a )
      {
        localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
        for( localIndex i = 0; i < 3; ++i )
        {
          xLocal[a][i] = nodeCoords( nodeIndex, i );
        }
      }

      constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
      for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
      {
        real32 const localIncrement = invC2 * m_finiteElement.computeMassTerm( q, xLocal );
        RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes( e, q )], localIncrement );
      }
    } ); // end loop over element
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};

template< typename FE_TYPE >
struct DampingMatrixKernel
{

  DampingMatrixKernel( FE_TYPE const & finiteElement )
    : m_finiteElement( finiteElement )
  {}

  /**
   * @brief Launches the precomputation of the damping matrices
   * @tparam EXEC_POLICY the execution policy
   * @tparam ATOMIC_POLICY the atomic policy
   * @param[in] size the number of cells in the subRegion
   * @param[in] nodeCoords coordinates of the nodes
   * @param[in] elemsToFaces map from faces to elements
   * @param[in] facesToNodes map from face to nodes
   * @param[in] facesDomainBoundaryIndicator flag equal to 1 if the face is on the boundary, and to 0 otherwise
   * @param[in] freeSurfaceFaceIndicator flag equal to 1 if the face is on the free surface, and to 0 otherwise
   * @param[in] lateralSurfaceFaceIndicator flag equal to 1 if the face is on a lateral surface, and to 0 otherwise
   * @param[in] bottomSurfaceFaceIndicator flag equal to 1 if the face is on the bottom surface, and to 0 otherwise
   * @param[in] velocity cell-wise velocity
   * @param[in] epsilon cell-wise Thomsen epsilon parameter
   * @param[in] delta cell-wise Thomsen delta parameter
   * @param[in] vti_f cell-wise f parameter in Fletcher's equation
   * @param[out] damping_p diagonal of the damping matrix for quantities in p in p-equation
   * @param[out] damping_q diagonal of the damping matrix for quantities in q in q-equation
   * @param[out] damping_pq diagonal of the damping matrix for quantities in q in p-equation
   * @param[out] damping_qp diagonal of the damping matrix for quantities in p in q-equation
   */
  template< typename EXEC_POLICY, typename ATOMIC_POLICY >
  void
  launch( localIndex const size,
          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const facesToNodes,
          arrayView1d< integer const > const facesDomainBoundaryIndicator,
          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
          arrayView1d< localIndex const > const lateralSurfaceFaceIndicator,
          arrayView1d< localIndex const > const bottomSurfaceFaceIndicator,
          arrayView1d< real32 const > const velocity,
          arrayView1d< real32 const > const epsilon,
          arrayView1d< real32 const > const delta,
          arrayView1d< real32 const > const vti_f,
          arrayView1d< real32 > const damping_p,
          arrayView1d< real32 > const damping_q,
          arrayView1d< real32 > const damping_pq,
          arrayView1d< real32 > const damping_qp )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
    {
      for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
      {
        localIndex const f = elemsToFaces( e, i );
        // face on the domain boundary and not on free surface
        if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
        {
          // ABC coefficients
          real32 alpha = 1.0 / velocity[e];
          // VTI coefficients
          real32 vti_p_xy= 0, vti_p_z = 0, vti_pq_z = 0;
          real32 vti_q_xy= 0, vti_q_z = 0, vti_qp_xy= 0;
          if( lateralSurfaceFaceIndicator[f] == 1 )
          {
            // ABC coefficients updated to fit horizontal velocity
            alpha /= sqrt( 1+2*epsilon[e] );

            vti_p_xy  = (1+2*epsilon[e]);
            vti_q_xy  = -(vti_f[e]-1);
            vti_qp_xy = (vti_f[e]+2*delta[e]);
          }
          if( bottomSurfaceFaceIndicator[f] == 1 )
          {
            // ABC coefficients updated to fit horizontal velocity
            alpha /= sqrt( 1+2*delta[e] );
            vti_p_z  = -(vti_f[e]-1);
            vti_pq_z = vti_f[e];
            vti_q_z  = 1;
          }
          // only the four corners of the mesh face are needed to compute the Jacobian
          real64 xLocal[ 4 ][ 3 ];
          for( localIndex a = 0; a < 4; ++a )
          {
            localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
            for( localIndex d = 0; d < 3; ++d )
            {
              xLocal[a][d] = nodeCoords( nodeIndex, d );
            }
          }

          constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
          for( localIndex q = 0; q < numNodesPerFace; ++q )
          {
            real32 const aux = m_finiteElement.computeDampingTerm( q, xLocal );
            real32 const localIncrement_p = alpha*(vti_p_xy + vti_p_z) * aux;
            RAJA::atomicAdd< ATOMIC_POLICY >( &damping_p[facesToNodes( f, q )], localIncrement_p );

            real32 const localIncrement_pq = alpha*vti_pq_z * aux;
            RAJA::atomicAdd< ATOMIC_POLICY >( &damping_pq[facesToNodes( f, q )], localIncrement_pq );

            real32 const localIncrement_q = alpha*(vti_q_xy + vti_q_z) * aux;
            RAJA::atomicAdd< ATOMIC_POLICY >( &damping_q[facesToNodes( f, q )], localIncrement_q );

            real32 const localIncrement_qp = alpha*vti_qp_xy * aux;
            RAJA::atomicAdd< ATOMIC_POLICY >( &damping_qp[facesToNodes( f, q )], localIncrement_qp );
          }
        }
      }
    } );
  }

  /// The finite element space/discretization object for the element type in the subRegion
  FE_TYPE const & m_finiteElement;

};



/**
 * @brief Implements kernels for solving the acoustic wave equations
 *   explicit central FD method and SEM
 * @copydoc geos::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 *
 * ### AcousticVTIWaveEquationSEMKernel Description
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
          typename FE_TYPE,
          typename Sp = fields::acousticvtifields::StiffnessVector_p, typename Sq = fields::acousticvtifields::StiffnessVector_q >
class ExplicitAcousticVTISEMBase : public finiteElement::KernelBase< SUBREGION_TYPE,
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
  ExplicitAcousticVTISEMBase( NodeManager & nodeManager,
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
    m_p_n( nodeManager.getField< fields::acousticvtifields::Pressure_p_n >() ),
    m_q_n( nodeManager.getField< fields::acousticvtifields::Pressure_q_n >() ),
    m_stiffnessVector_p( nodeManager.getField< Sp >() ),
    m_stiffnessVector_q( nodeManager.getField< Sq >() ),
    m_epsilon( elementSubRegion.template getField< fields::acousticvtifields::Epsilon >() ),
    m_delta( elementSubRegion.template getField< fields::acousticvtifields::Delta >() ),
    m_vti_f( elementSubRegion.template getField< fields::acousticvtifields::F >() ),
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
   * ### ExplicitAcousticVTISEMBase Description
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
    /// only the eight corners of the mesh cell are needed to compute the Jacobian
    real64 xLocal[ 8 ][ 3 ];
    real32 factor;
    /// local (to this element) stiffness vectors
    real32 stiffnessVectorLocal_p[ numNodesPerElem ]{};
    real32 stiffnessVectorLocal_q[ numNodesPerElem ]{};
  };
  //***************************************************************************


  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    /// numDofPerTrialSupportPoint = 1
    for( localIndex a=0; a< 8; a++ )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, FE_TYPE::meshIndexToLinearIndex3D( a ) );
      for( int i=0; i< 3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_nodeCoords[ nodeIndex ][ i ];
      }
    }
    stack.factor = 1.0;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    for( int i=0; i<numNodesPerElem; i++ )
    {
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_p[m_elemsToNodes[k][i]], stack.stiffnessVectorLocal_p[i] );
      RAJA::atomicAdd< parallelDeviceAtomic >( &m_stiffnessVector_q[m_elemsToNodes[k][i]], stack.stiffnessVectorLocal_q[i] );
    }
    return 0;
  }
  /**
   * @copydoc geos::finiteElement::KernelBase::quadraturePointKernel
   *
   * ### ExplicitAcousticVTISEMBase Description
   * Calculates stiffness vector
   *
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Pseudo Stiffness xy
    m_finiteElementSpace.template computeStiffnessxyTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 const localIncrement_p = val*(-1-2*m_epsilon[k])*m_p_n[m_elemsToNodes[k][j]];
      stack.stiffnessVectorLocal_p[ i ] += localIncrement_p * stack.factor;

      real32 const localIncrement_q = val*((-2*m_delta[k]-m_vti_f[k])*m_p_n[m_elemsToNodes[k][j]] +(m_vti_f[k]-1)*m_q_n[m_elemsToNodes[k][j]]);
      stack.stiffnessVectorLocal_q[ i ] += localIncrement_q * stack.factor;
    } );

    // Pseudo-Stiffness z

    m_finiteElementSpace.template computeStiffnesszTerm( q, stack.xLocal, [&] ( int i, int j, real64 val )
    {
      real32 const localIncrement_p = val*((m_vti_f[k]-1)*m_p_n[m_elemsToNodes[k][j]] - m_vti_f[k]*m_q_n[m_elemsToNodes[k][j]]);
      stack.stiffnessVectorLocal_p[ i ] += localIncrement_p * stack.factor;

      real32 const localIncrement_q = -val*m_q_n[m_elemsToNodes[k][j]];
      stack.stiffnessVectorLocal_q[ i ] += localIncrement_q * stack.factor;
    } );
  }

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

  /// The array containing the epsilon Thomsen parameter.
  arrayView1d< real32 const > const m_epsilon;

  /// The array containing the delta Thomsen parameter.
  arrayView1d< real32 const > const m_delta;

  /// The array containing the f parameter.
  arrayView1d< real32 const > const m_vti_f;

  /// The time increment for this time integration step.
  real64 const m_dt;


};


/// Specialization for standard iso elastic kernel
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
using ExplicitAcousticVTISEM = ExplicitAcousticVTISEMBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >;
using ExplicitAcousticVTISEMFactory = finiteElement::KernelFactory< ExplicitAcousticVTISEM,
                                                                    real64 >;
/// Specialization for attenuation kernel
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ExplicitAcousticVTIAttenuativeSEM : public ExplicitAcousticVTISEMBase< SUBREGION_TYPE,
                                                                             CONSTITUTIVE_TYPE,
                                                                             FE_TYPE,
                                                                             fields::acousticvtifields::StiffnessVectorA_p,
                                                                             fields::acousticvtifields::StiffnessVectorA_q >
{
public:

  /// Alias for the base class;
  using Base = ExplicitAcousticVTISEMBase< SUBREGION_TYPE,
                                           CONSTITUTIVE_TYPE,
                                           FE_TYPE,
                                           fields::acousticvtifields::StiffnessVectorA_p,
                                           fields::acousticvtifields::StiffnessVectorA_q >;

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
  ExplicitAcousticVTIAttenuativeSEM( NodeManager & nodeManager,
                                     EdgeManager const & edgeManager,
                                     FaceManager const & faceManager,
                                     localIndex const targetRegionIndex,
                                     SUBREGION_TYPE const & elementSubRegion,
                                     FE_TYPE const & finiteElementSpace,
                                     CONSTITUTIVE_TYPE & inputConstitutiveType,
                                     real64 const dt ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          dt ),
    m_qualityFactor( elementSubRegion.template getField< fields::acousticfields::AcousticQualityFactor >() )
  {}

  /**
   * @copydoc geos::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              typename Base::StackVariables & stack ) const
  {
    Base::setup( k, stack );
    stack.factor = stack.factor / m_qualityFactor[ k ];
  }

protected:

  /// The array containing the acoustic attenuation quality factor
  arrayView1d< real32 const > const m_qualityFactor;

};

using ExplicitAcousticVTIAttenuativeSEMFactory = finiteElement::KernelFactory< ExplicitAcousticVTIAttenuativeSEM,
                                                                               real64 >;

} // namespace acousticVTIWaveEquationSEMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICVTIWAVEEQUATIONSEMKERNEL_HPP_
