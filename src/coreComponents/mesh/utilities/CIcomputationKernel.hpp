/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file CIcomputationKernel.hpp
 */
#ifndef GEOS_MESH_UTILITIES_CICOMPUTATIONKERNEL_HPP_
#define GEOS_MESH_UTILITIES_CICOMPUTATIONKERNEL_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/EmbeddedSurfaceSubRegion.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/ElementType.hpp"


namespace geos
{

/**
 * @brief Kernel to compute EDFM connectivity index
 *
 * @tparam FE_TYPE finite element space (element type) on which CI are computed
 */
template< typename FE_TYPE >
class CIcomputationKernel
{
public:
  /**
   * @brief Construct a new CIcomputationKernel object
   *
   * @param finiteElementSpace the finite element space
   * @param nodeManager the nodeManager
   * @param elementSubRegion the element subRegion
   * @param embeddedSurfSubRegion the embeddedSurfaceSubRegion
   */
  CIcomputationKernel( FE_TYPE const & finiteElementSpace,
                       NodeManager const & nodeManager,
                       CellElementSubRegion const & elementSubRegion,
                       EmbeddedSurfaceSubRegion & embeddedSurfSubRegion ):
    m_finiteElementSpace( finiteElementSpace ),
    m_elementType( elementSubRegion.getElementType() ),
    m_X( nodeManager.referencePosition() ),
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_fracturedElems( elementSubRegion.fracturedElementsList().toViewConst()),
    m_cellsToEmbeddedSurfaces( elementSubRegion.embeddedSurfacesList().toViewConst() ),
    m_normalVector( embeddedSurfSubRegion.getNormalVector().toViewConst()),
    m_fracCenter( embeddedSurfSubRegion.getElementCenter().toViewConst()),
    m_fractureSurfaceArea( embeddedSurfSubRegion.getElementArea().toViewConst() ),
    m_connectivityIndex( embeddedSurfSubRegion.getConnectivityIndex() )
  {}

  /// number of nodes per element
  static constexpr int numNodesPerElem = FE_TYPE::maxSupportPoints;

  /// number of sampling points per element
  static constexpr int numSamplingPoints = FE_TYPE::numSamplingPoints;

  /**
   * @brief stack variables
   *
   */
  struct StackVariables
  {
public:

    /// Constructor.
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal(),
      samplingPointCoord()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ numNodesPerElem ][ 3 ];

    /// C-array stack storage for sampling points coordinates.
    real64 samplingPointCoord[3];
  };

  /**
   * @brief launch of CI calculation
   *
   * @tparam POLICY the exectution policy
   * @tparam KERNEL_TYPE the type of kernel
   * @param kernelComponent the kernel object
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  void launchCIComputationKernel( KERNEL_TYPE & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( kernelComponent.m_fracturedElems.size(),
                      [=] GEOS_HOST_DEVICE ( localIndex const i )
    {

      localIndex k = kernelComponent.m_fracturedElems[i];

      StackVariables stack;
      kernelComponent.setup( k, stack );

      real64 averageDistance = 0.0;
      for( integer np=0; np<numSamplingPoints; ++np )
      {
        kernelComponent.samplingPointCoord( np, stack );
        averageDistance += kernelComponent.computeDistance( k, stack.samplingPointCoord );
      }
      averageDistance /= numSamplingPoints;
      kernelComponent.setConnectivityIndex( k, averageDistance );
    } );
  }

/**
 * @brief set up the kernel object by copying global values in the stack.
 *
 * @param k embedded surface index.
 * @param stack stack variables
 */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
      }
    }
  }

  /**
   * @brief computes the distance between a point inside the element and the cut surface \p k.
   *
   * @param k the index of the embedded surface
   * @param point the coordinates of a point inside the cell element
   * @return the distance
   */
  GEOS_HOST_DEVICE
  real64 computeDistance( localIndex const k,
                          real64 const (&point)[3] ) const
  {
    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];
    real64 pointToFracCenter[3];
    LvArray::tensorOps::copy< 3 >( pointToFracCenter, point );
    LvArray::tensorOps::subtract< 3 >( pointToFracCenter, m_fracCenter[embSurfIndex] );
    return LvArray::math::abs( LvArray::tensorOps::AiBi< 3 >( pointToFracCenter, m_normalVector[embSurfIndex] ));
  }

  /**
   * @brief computes coordinates of the sampling point in the physical space.
   *
   * @param np sampling point linear index
   * @param stack stack variables
   */
  GEOS_HOST_DEVICE
  void samplingPointCoord( integer const np,
                           StackVariables & stack ) const
  {
    // Get sampling point coord in parent space.
    real64 parentSamplingPointCoord[3] = {0.0, 0.0, 0.0};
    FE_TYPE::getSamplingPointCoordInParentSpace( np, parentSamplingPointCoord );

    // Compute shape function values at sampling point
    real64 N[numNodesPerElem];
    FE_TYPE::calcN( parentSamplingPointCoord, N );

    LvArray::tensorOps::fill< 3 >( stack.samplingPointCoord, 0.0 );

    // Compute sampling point coord in the physical space
    for( localIndex a=0; a<numNodesPerElem; a++ )
    {
      for( int i =0; i < 3; i++ )
      {
        stack.samplingPointCoord[i] += stack.xLocal[a][i] * N[a];
      }
    }
  }

  /**
   * @brief Set the Connectivity Index
   *
   * @param k embedded surface element index
   * @param averageDistance the average distance
   */
  GEOS_HOST_DEVICE
  void setConnectivityIndex( localIndex const k,
                             real64 const averageDistance ) const
  {
    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];
    m_connectivityIndex[embSurfIndex] = 2. * m_fractureSurfaceArea[embSurfIndex] / averageDistance;
  }

private:

  /// the finite element space
  FE_TYPE const & m_finiteElementSpace;

  /// the element type
  ElementType const m_elementType;

  /// the reference position of the nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The element to nodes map.
  traits::ViewTypeConst< typename CellElementSubRegion::NodeMapType::base_type > const m_elemsToNodes;

  /// set of fractured cell elements
  SortedArrayView< localIndex const > const m_fracturedElems;

  /// cell to embedded surfaces map
  ArrayOfArraysView< localIndex const > const m_cellsToEmbeddedSurfaces;

  /// normal vector of the embedded surface
  arrayView2d< real64 const > const m_normalVector;

  /// center of the embedded surface
  arrayView2d< real64 const > const m_fracCenter;

  /// area of the embedded surface
  arrayView1d< real64 const > const m_fractureSurfaceArea;

  /// the connectivity index
  arrayView1d< real64 > const m_connectivityIndex;
};

}

#endif /* GEOS_MESH_UTILITIES_CICOMPUTATIONKERNEL_HPP_ */
