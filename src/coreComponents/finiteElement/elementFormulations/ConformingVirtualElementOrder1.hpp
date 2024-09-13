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
 * @file ConformingVirtualElementOrder1.hpp
 */

#ifndef GEOS_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
#define GEOS_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_

#include "FiniteElementBase.hpp"
#include "codingUtilities/traits.hpp"

namespace geos
{
namespace finiteElement
{
/**
 * This class contains the kernel accessible functions specific to the H1-conforming nodal Virtual
 * Element Method of order 1, with a 1-point Gaussian quadrature rule.
 *
 * @tparam MAXCELLNODES The maximum number of nodes per cell that the class expects (used to
 * pre-allocate stack arrays).
 * @tparam MAXFACENODES The maximum number of nodes per face that the class expects.
 */
template< localIndex MAXCELLNODES, localIndex MAXFACENODES >
class ConformingVirtualElementOrder1 final : public FiniteElementBase
{
public:
  /// Type of MeshData::nodesCoords.
  using InputNodeCoords = arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD >;
  /// Type of MeshData::cellToNodeMap.
  template< typename SUBREGION_TYPE >
  using InputCellToNodeMap = traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType >;
  /// Type of MeshData::cellToFaceMap.
  template< typename SUBREGION_TYPE >
  using InputCellToFaceMap = traits::ViewTypeConst< typename SUBREGION_TYPE::FaceMapType >;
  /// Type of MeshData::faceToNodeMap.
  using InputFaceToNodeMap = ArrayOfArraysView< localIndex const >;
  /// Type of MeshData::faceToEdgeMap.
  using InputFaceToEdgeMap = ArrayOfArraysView< localIndex const >;
  /// Type of MeshData::edgeToNodeMap.
  using InputEdgeToNodeMap = arrayView2d< localIndex const >;

  /// The maximum number of support points per element.
  static constexpr localIndex maxSupportPoints = MAXCELLNODES;
  /// Static property kept for consistency with other finite element classes.
  static constexpr localIndex numNodes = MAXCELLNODES;
  /// The number of quadrature points per element.
  static constexpr localIndex numQuadraturePoints = 1;

  ConformingVirtualElementOrder1() = default;

  virtual ~ConformingVirtualElementOrder1() = default;

  /**
   * @struct StackVariables
   * @brief Kernel variables allocated on the stack.
   *
   * It holds the computed projections of basis functions and basis function derivatives and the
   * stabilization matrix. Arrays are pre-allocated using @p MAXCELLNODES.
   * @sa setupStack.
   */
  struct StackVariables : public FiniteElementBase::StackVariables
  {
    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables()
    {}

    /// The number of support points.
    localIndex numSupportPoints;
    /// The quadrature weight.
    real64 quadratureWeight;
    /// Array holding the integral mean of basis functions in the first @ref numSupportPoints
    /// positions.
    real64 basisFunctionsIntegralMean[MAXCELLNODES];
    /// The stabilization matrix. Valid values will be in the upper left @ref numSupportPoints x
    /// @ref numSupportPoints block.
    real64 stabilizationMatrix[MAXCELLNODES][MAXCELLNODES];
    /// Array holding the integral mean of derivatives of basis functions in the first @ref
    /// numSupportPoints position of the first dimension.
    real64 basisDerivativesIntegralMean[MAXCELLNODES][3];
  };

  /**
   * @struct MeshData
   * @brief Variables used to call the @ref setupStack method.
   * @tparam SUBREGION_TYPE The type of mesh sub-region.
   */
  template< typename SUBREGION_TYPE >
  struct MeshData : public FiniteElementBase::MeshData< SUBREGION_TYPE >
  {
    /**
     * Constructor
     */
    MeshData()
    {}

    /// View to the array containing nodes coordinates.
    InputNodeCoords nodesCoords;
    /// View to the cell-to-node map in the sub-region.
    InputCellToNodeMap< SUBREGION_TYPE > cellToNodeMap;
    /// View to the cell-to-face map in the sub-region.
    InputCellToFaceMap< SUBREGION_TYPE > cellToFaceMap;
    /// View to the face-to-node map in the sub-region.
    InputFaceToNodeMap faceToNodeMap;
    /// View to the face-to-edge map in the sub-region.
    InputFaceToEdgeMap faceToEdgeMap;
    /// View to the edge-to-node map in the sub-region.
    InputEdgeToNodeMap edgeToNodeMap;
    /// View to the array of face centers.
    arrayView2d< real64 const > faceCenters;
    /// View to the array of face normals.
    arrayView2d< real64 const > faceNormals;
    /// View to the array of face areas.
    arrayView1d< real64 const > faceAreas;
    /// View to the array of cell centers.
    arrayView2d< real64 const > cellCenters;
    /// View to the array of cell volumes.
    arrayView1d< real64 const > cellVolumes;
  };

  GEOS_HOST_DEVICE
  inline
  localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  GEOS_HOST_DEVICE
  inline
  virtual localIndex getMaxSupportPoints() const override
  {
    return maxSupportPoints;
  }

  /**
   * @brief Get the number of support points.
   * @param stack Object that holds stack variables.
   * @return The number of support points.
   */
  GEOS_HOST_DEVICE
  inline
  static localIndex getNumSupportPoints( StackVariables const & stack )
  {
    return stack.numSupportPoints;
  }

  /**
   * @brief Calculate the shape functions projected derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param gradN Array to contain the shape function projected derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  inline
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[maxSupportPoints][3],
                           StackVariables const & stack,
                           real64 ( & gradN )[maxSupportPoints][3] )
  {
    for( localIndex i = 0; i < stack.numSupportPoints; ++i )
    {
      gradN[i][0] = stack.basisDerivativesIntegralMean[i][0];
      gradN[i][1] = stack.basisDerivativesIntegralMean[i][1];
      gradN[i][2] = stack.basisDerivativesIntegralMean[i][2];
    }
    return transformedQuadratureWeight( q, X, stack );
  }

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOS_HOST_DEVICE
  inline
  static void calcN( localIndex const q,
                     StackVariables const & stack,
                     real64 ( & N )[maxSupportPoints] )
  {
    GEOS_UNUSED_VAR( q );
    for( localIndex i = 0; i < stack.numSupportPoints; ++i )
    {
      N[i] = stack.basisFunctionsIntegralMean[i];
    }
  }

  /**
   * @brief Adds a grad-grad stabilization to @p matrix.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @tparam MAXSUPPORTPOINTS Maximum number of support points allowed for this element.
   * @tparam UPPER If true only the upper triangular part of @p matrix is modified.
   * @param stack Stack variables as filled by @ref setupStack.
   * @param matrix The matrix that needs to be stabilized.
   * @param scaleFactor Scaling of the stabilization matrix.
   */
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT, localIndex MAXSUPPORTPOINTS, bool UPPER >
  GEOS_HOST_DEVICE
  inline
  static void addGradGradStabilization( StackVariables const & stack,
                                        real64 ( & matrix )
                                        [MAXSUPPORTPOINTS * NUMDOFSPERTRIALSUPPORTPOINT]
                                        [MAXSUPPORTPOINTS * NUMDOFSPERTRIALSUPPORTPOINT],
                                        real64 const & scaleFactor )
  {
    for( localIndex i = 0; i < stack.numSupportPoints; ++i )
    {
      localIndex startCol = (UPPER) ? i : 0;
      for( localIndex j = startCol; j < stack.numSupportPoints; ++j )
      {
        real64 const contribution = scaleFactor * stack.stabilizationMatrix[i][j];
        for( localIndex d = 0; d < NUMDOFSPERTRIALSUPPORTPOINT; ++d )
        {
          matrix[i*NUMDOFSPERTRIALSUPPORTPOINT + d][j*NUMDOFSPERTRIALSUPPORTPOINT + d] += contribution;
        }
      }
    }
  }

  /**
   * @brief Adds a grad-grad stabilization evaluated at @p dofs to @p targetVector.
   * @details This method is intended to be used with @p targetVector being the residual and @p dofs
   * being the degrees of freedom of the previous solution.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @param stack Stack variables as filled by @ref setupStack.
   * @param dofs The degrees of freedom of the function where the stabilization operator has to be
   * evaluated.
   * @param targetVector The input vector to which values have to be added, seen in chunks of length
   * @p NUMDOFSPERTRIALSUPPORTPOINT.
   * @param scaleFactor Scaling of the stabilization matrix.
   */
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT, localIndex MAXSUPPORTPOINTS >
  GEOS_HOST_DEVICE
  inline
  static void addEvaluatedGradGradStabilization( StackVariables const & stack,
                                                 real64 const ( &dofs )[MAXSUPPORTPOINTS][NUMDOFSPERTRIALSUPPORTPOINT],
                                                 real64 ( & targetVector )[MAXSUPPORTPOINTS][NUMDOFSPERTRIALSUPPORTPOINT],
                                                 real64 const scaleFactor )
  {
    for( localIndex d = 0; d < NUMDOFSPERTRIALSUPPORTPOINT; ++d )
    {
      for( localIndex i = 0; i < stack.numSupportPoints; ++i )
      {
        for( localIndex j = 0; j < stack.numSupportPoints; ++j )
        {
          targetVector[i][d] += scaleFactor * stack.stabilizationMatrix[i][j] * dofs[j][d];
        }
      }
    }
  }

  /**
   * @brief Calculate the integration weights for a quadrature point.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param stack Stack variables as filled by @ref setupStack.
   * @return The product of the quadrature rule weight and the determinate of
   *   the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  inline
  static real64 transformedQuadratureWeight( localIndex const q,
                                             real64 const ( &X )[maxSupportPoints][3],
                                             StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( q, X );
    return stack.quadratureWeight;
  }

  /**
   * @brief Method to fill a MeshData object.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param cellSubRegion The cell sub-region for which the element has to be initialized.
   * @param meshData MeshData struct to be filled.
   */
  template< typename SUBREGION_TYPE >
  inline
  static void fillMeshData( NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            SUBREGION_TYPE const & cellSubRegion,
                            MeshData< SUBREGION_TYPE > & meshData
                            )
  {
    meshData.nodesCoords = nodeManager.referencePosition();
    meshData.cellToNodeMap = cellSubRegion.nodeList().toViewConst();
    meshData.cellToFaceMap = cellSubRegion.faceList().toViewConst();
    meshData.faceToNodeMap = faceManager.nodeList().toViewConst();
    meshData.faceToEdgeMap = faceManager.edgeList().toViewConst();
    meshData.edgeToNodeMap = edgeManager.nodeList().toViewConst();
    meshData.faceCenters = faceManager.faceCenter();
    meshData.faceNormals = faceManager.faceNormal();
    meshData.faceAreas = faceManager.faceArea();
    meshData.cellCenters = cellSubRegion.getElementCenter();
    meshData.cellVolumes = cellSubRegion.getElementVolume();
  }

  /**
   * @brief Setup method.
   * @param cellIndex The index of the cell with respect to the cell sub region to which the element
   * has been initialized previously (see @ref fillMeshData).
   * @param meshData Object previously initialized by @ref fillMeshData.
   * @param stack Object that holds stack variables.
   */
  template< typename SUBREGION_TYPE >
  GEOS_HOST_DEVICE
  inline
  static void setupStack( localIndex const & cellIndex,
                          MeshData< SUBREGION_TYPE > const & meshData,
                          StackVariables & stack )
  {
    real64 const cellCenter[3] { meshData.cellCenters( cellIndex, 0 ),
                                 meshData.cellCenters( cellIndex, 1 ),
                                 meshData.cellCenters( cellIndex, 2 ) };
    real64 const cellVolume = meshData.cellVolumes( cellIndex );

    computeProjectors< SUBREGION_TYPE >( cellIndex,
                                         meshData.nodesCoords,
                                         meshData.cellToNodeMap,
                                         meshData.cellToFaceMap,
                                         meshData.faceToNodeMap,
                                         meshData.faceToEdgeMap,
                                         meshData.edgeToNodeMap,
                                         meshData.faceCenters,
                                         meshData.faceNormals,
                                         meshData.faceAreas,
                                         cellCenter,
                                         cellVolume,
                                         stack.numSupportPoints,
                                         stack.quadratureWeight,
                                         stack.basisFunctionsIntegralMean,
                                         stack.stabilizationMatrix,
                                         stack.basisDerivativesIntegralMean );
  }

  /**
   * @defgroup DeprecatedSyntax VEM functions with deprecated syntax.
   *
   * Functions that are implemented for consistency with other FEM classes but will issue an error
   * if called.
   *
   * @{
   */

  /**
   * @brief This function returns an error, since to get the number of support points with VEM you
   * have to use the StackVariables version of this function.
   * @return Zero.
   */
  GEOS_HOST_DEVICE
  inline
  localIndex getNumSupportPoints() const override
  {
    GEOS_ERROR( "VEM functions have to be called with the StackVariables syntax" );
    return 0;
  }


  /**
   * @brief Get the Sampling Point Coord In the Parent Space
   *
   * @param linearIndex linear index of the sampling point
   * @param samplingPointCoord coordinates of the sampling point
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void getSamplingPointCoordInParentSpace( int const & linearIndex,
                                                  real64 (& samplingPointCoord)[3] )
  {
    GEOS_UNUSED_VAR( linearIndex, samplingPointCoord );
    GEOS_ERROR( "Element type not supported." );
  }

  /**
   * @brief This function returns an error, since to get projection of basis functions with VEM you
   * have to use the StackVariables version of this function.
   * @param q The quadrature point index in 3d space.
   * @param N Array to store the values of shape functions, that is actually set to zero.
   */
  GEOS_HOST_DEVICE
  inline
  static void calcN( localIndex const q,
                     real64 ( & N )[maxSupportPoints] )
  {
    GEOS_ERROR( "VEM functions have to be called with the StackVariables syntax" );
    GEOS_UNUSED_VAR( q );
    for( localIndex i = 0; i < maxSupportPoints; ++i )
    {
      N[i] = 0.0;
    }
  }

  /**
   * @brief This function returns an error, since to get projection of gradients with VEM you have
   * to use the StackVariables version of this function.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value
   * @param[out] N The shape function values. It is set to zero.
   */
  GEOS_HOST_DEVICE
  inline
  static void calcN( real64 const (&coords)[3],
                     real64 (& N)[maxSupportPoints] )
  {
    GEOS_ERROR( "VEM functions have to be called with the StackVariables syntax" );
    GEOS_UNUSED_VAR( coords );
    for( localIndex i = 0; i < maxSupportPoints; ++i )
    {
      N[i] = 0.0;
    }
  }

  /**
   * @brief This function returns an error, since there is no reference element defined in the VEM
   * context. It is kept for consistency with other finite element classes.
   * @param q The quadrature point index in 3d space.
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   * @return A zero matrix.
   */
  GEOS_HOST_DEVICE
  inline
  static real64 invJacobianTransformation( int const q,
                                           real64 const ( &X )[numNodes][3],
                                           real64 ( & J )[3][3] )
  {
    GEOS_ERROR( "No reference element map is defined for VEM classes" );
    GEOS_UNUSED_VAR( q, X );
    for( localIndex i = 0; i < 3; ++i )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        J[i][j] = 0.0;
      }
    }
    return 0.0;
  }

  /**
   * @brief This function returns an error, since to get projection of gradients with VEM you have
   * to use the StackVariables version of this function.
   * @param q The quadrature point index in 3d space.
   * @param X Array containing the coordinates of the support points.
   * @param gradN Array to store the gradients of shape functions. It is set to zero.
   * @return Zero.
   */
  GEOS_HOST_DEVICE
  inline
  static real64 calcGradN( localIndex const q,
                           real64 const ( &X )[maxSupportPoints][3],
                           real64 ( & gradN )[maxSupportPoints][3] )
  {
    GEOS_ERROR( "VEM functions have to be called with the StackVariables syntax" );
    GEOS_UNUSED_VAR( q, X );
    for( localIndex i = 0; i < maxSupportPoints; ++i )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        gradN[i][j] = 0.0;
      }
    }
    return 0.0;
  }

  /**
   * @brief This function returns an error, since to get the quadrature weight with VEM you have to
   * use the StackVariables version of this function.
   * @param q The quadrature point index in 3d space.
   * @param X Array containing the coordinates of the support points.
   * @return Zero.
   */
  GEOS_HOST_DEVICE
  real64 transformedQuadratureWeight( localIndex const q,
                                      real64 const ( &X )[maxSupportPoints][3] ) const
  {
    GEOS_ERROR( "VEM functions have to be called with the StackVariables syntax" );
    GEOS_UNUSED_VAR( q, X );
    return 0.0;
  }

  /** @} */

private:

  GEOS_HOST_DEVICE
  static void
    computeFaceIntegrals( InputNodeCoords const & nodesCoords,
                          localIndex const (&faceToNodes)[MAXFACENODES],
                          localIndex const (&faceToEdges)[MAXFACENODES],
                          localIndex const & numFaceVertices,
                          real64 const & faceArea,
                          real64 const (&faceCenter)[3],
                          real64 const (&faceNormal)[3],
                          InputEdgeToNodeMap const & edgeToNodes,
                          real64 const & invCellDiameter,
                          real64 const (&cellCenter)[3],
                          real64 ( &basisIntegrals )[MAXFACENODES],
                          real64 ( &threeDMonomialIntegrals )[3] );

  /**
   * @brief This function contains the kernel computations required to compute projections of
   * basis functions and basis function derivatives.
   *
   * It is called by @ref setupStack and outputs are written on the suitable fields of a
   * StackVariables struct.
   * @tparam SUBREGION_TYPE The type of subregion to be processed.
   * @param cellIndex The index of the cell to be processed
   */
  template< typename SUBREGION_TYPE >
  GEOS_HOST_DEVICE
  inline
  static void
    computeProjectors( localIndex const & cellIndex,
                       InputNodeCoords const & nodesCoords,
                       InputCellToNodeMap< SUBREGION_TYPE > const & cellToNodeMap,
                       InputCellToFaceMap< SUBREGION_TYPE > const & elementToFaceMap,
                       InputFaceToNodeMap const & faceToNodeMap,
                       InputFaceToEdgeMap const & faceToEdgeMap,
                       InputEdgeToNodeMap const & edgeToNodeMap,
                       arrayView2d< real64 const > const faceCenters,
                       arrayView2d< real64 const > const faceNormals,
                       arrayView1d< real64 const > const faceAreas,
                       real64 const (&cellCenter)[3],
                       real64 const & cellVolume,
                       localIndex & numSupportPoints,
                       real64 & quadratureWeight,
                       real64 ( &basisFunctionsIntegralMean )[MAXCELLNODES],
                       real64 ( &stabilizationMatrix )[MAXCELLNODES][MAXCELLNODES],
                       real64 ( &basisDerivativesIntegralMean )[MAXCELLNODES][3]
                       );

  template< localIndex DIMENSION, typename POINT_COORDS_TYPE >
  GEOS_HOST_DEVICE
  inline
  static real64 computeDiameter( POINT_COORDS_TYPE points,
                                 localIndex const & numPoints )
  {
    real64 diameter = 0;
    for( localIndex numPoint = 0; numPoint < numPoints; ++numPoint )
    {
      for( localIndex numOthPoint = 0; numOthPoint < numPoint; ++numOthPoint )
      {
        real64 candidateDiameter = 0.0;
        for( localIndex i = 0; i < DIMENSION; ++i )
        {
          real64 coordDiff = points[numPoint][i] - points[numOthPoint][i];
          candidateDiameter += coordDiff * coordDiff;
        }
        if( diameter < candidateDiameter )
        {
          diameter = candidateDiameter;
        }
      }
    }
    return LvArray::math::sqrt< real64 >( diameter );
  }

  template< localIndex DIMENSION, typename POINT_COORDS_TYPE, typename POINT_SELECTION_TYPE >
  GEOS_HOST_DEVICE
  inline
  static real64 computeDiameter( POINT_COORDS_TYPE points,
                                 POINT_SELECTION_TYPE selectedPoints,
                                 localIndex const & numSelectedPoints )
  {
    real64 diameter = 0;
    for( localIndex numPoint = 0; numPoint < numSelectedPoints; ++numPoint )
    {
      for( localIndex numOthPoint = 0; numOthPoint < numPoint; ++numOthPoint )
      {
        real64 candidateDiameter = 0.0;
        for( localIndex i = 0; i < DIMENSION; ++i )
        {
          real64 coordDiff = points[selectedPoints[numPoint]][i] -
                             points[selectedPoints[numOthPoint]][i];
          candidateDiameter += coordDiff * coordDiff;
        }
        if( diameter < candidateDiameter )
        {
          diameter = candidateDiameter;
        }
      }
    }
    return LvArray::math::sqrt< real64 >( diameter );
  }
};

/// Convenience typedef for VEM on tetrahedra.
using H1_Tetrahedron_VEM_Gauss1 = ConformingVirtualElementOrder1< 4, 3 >;
#if !defined( GEOS_USE_HIP )
/// Convenience typedef for VEM on hexahedra.
using H1_Hexahedron_VEM_Gauss1 = ConformingVirtualElementOrder1< 8, 4 >;
#endif
/// Convenience typedef for VEM on pyramids.
using H1_Pyramid_VEM_Gauss1 = ConformingVirtualElementOrder1< 5, 4 >;
#if !defined( GEOS_USE_HIP )
/// Convenience typedef for VEM on wedges.
using H1_Wedge_VEM_Gauss1 = ConformingVirtualElementOrder1< 6, 4 >;
#endif
/// Convenience typedef for VEM on prism5.
using H1_Prism5_VEM_Gauss1 = ConformingVirtualElementOrder1< 10, 5 >;
/// Convenience typedef for VEM on prism6.
using H1_Prism6_VEM_Gauss1 = ConformingVirtualElementOrder1< 12, 6 >;
/// Convenience typedef for VEM on prism7.
using H1_Prism7_VEM_Gauss1 = ConformingVirtualElementOrder1< 14, 7 >;
/// Convenience typedef for VEM on prism8.
using H1_Prism8_VEM_Gauss1 = ConformingVirtualElementOrder1< 16, 8 >;
/// Convenience typedef for VEM on prism9.
using H1_Prism9_VEM_Gauss1 = ConformingVirtualElementOrder1< 18, 9 >;
/// Convenience typedef for VEM on prism10.
using H1_Prism10_VEM_Gauss1 = ConformingVirtualElementOrder1< 20, 10 >;
/// Convenience typedef for VEM on prism11.
#if !defined( GEOS_USE_HIP )
using H1_Prism11_VEM_Gauss1 = ConformingVirtualElementOrder1< 22, 11 >;
#endif
}
}

#include "ConformingVirtualElementOrder1_impl.hpp"

#endif // GEOS_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
