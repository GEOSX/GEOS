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
 * @file ComputationalGeometry.hpp
 */

#ifndef GEOS_MESH_UTILITIES_COMPUTATIONALGEOMETRY_HPP_
#define GEOS_MESH_UTILITIES_COMPUTATIONALGEOMETRY_HPP_

#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "finiteElement/elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp"
#include "finiteElement/elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp"
#include "finiteElement/elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp"
#include "finiteElement/elementFormulations/H1_Wedge_Lagrange1_Gauss6.hpp"
#include "LvArray/src/output.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{
namespace computationalGeometry
{

/// Machine epsilon for double-precision calculations
constexpr real64 machinePrecision = LvArray::NumericLimits< real64 >::epsilon;

/**
 * @brief Calculate the intersection between a line and a plane.
 * @tparam LINEDIR_TYPE the type of @p lineDir
 * @tparam POINT_TYPE the type of @p linePoint
 * @tparam NORMAL_TYPE the type of @p planeNormal
 * @tparam ORIGIN_TYPE the type of @p planeOrigin
 * @tparam INTPOINT_TYPE the type of @p instersectionPoint
 * @param[in] lineDir vector defining direction of the line
 * @param[in] linePoint one point of the line
 * @param[in] planeNormal normal to plane
 * @param[in] planeOrigin plane origin
 * @param[out] intersectionPoint the intersection point
 */
template< typename LINEDIR_TYPE,
          typename POINT_TYPE,
          typename NORMAL_TYPE,
          typename ORIGIN_TYPE,
          typename INTPOINT_TYPE >
void LinePlaneIntersection( LINEDIR_TYPE const & lineDir,
                            POINT_TYPE const & linePoint,
                            NORMAL_TYPE const & planeNormal,
                            ORIGIN_TYPE const & planeOrigin,
                            INTPOINT_TYPE & intersectionPoint )
{
  /* Find intersection line plane
   * line equation: p - (d*lineDir + linePoing) = 0;
   * plane equation: ( p - planeOrigin) * planeNormal = 0;
   * d = (planeOrigin - linePoint) * planeNormal / (lineDir * planeNormal )
   * pInt = d*lineDir+linePoint;
   */
  real64 dummy[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( planeOrigin );
  LvArray::tensorOps::subtract< 3 >( dummy, linePoint );
  real64 const d = LvArray::tensorOps::AiBi< 3 >( dummy, planeNormal ) /
                   LvArray::tensorOps::AiBi< 3 >( lineDir, planeNormal );

  LvArray::tensorOps::copy< 3 >( intersectionPoint, linePoint );
  LvArray::tensorOps::scaledAdd< 3 >( intersectionPoint, lineDir, d );
}

/**
 * @brief Reorder a set of points counter-clockwise.
 * @tparam NORMAL_TYPE the type of @p normal
 * @param[in] points coordinates of the points
 * @param[in] normal vector normal to the plane
 * @return an std::vector containing the original indices of the reordered points.
 */
template< typename NORMAL_TYPE >
array1d< int >  orderPointsCCW( arrayView2d< real64 > const & points,
                                NORMAL_TYPE const & normal )
{
  localIndex const numPoints = points.size( 0 );

  array2d< real64 > orderedPoints( numPoints, 3 );

  array1d< int > indices( numPoints );
  array1d< real64 > angle( numPoints );

  // compute centroid of the set of points
  real64 centroid[3];
  LvArray::tensorOps::fill< 3 >( centroid, 0 );
  for( localIndex a = 0; a < numPoints; ++a )
  {
    LvArray::tensorOps::add< 3 >( centroid, points[ a ] );
    indices[ a ] = a;
  }

  LvArray::tensorOps::scale< 3 >( centroid, 1.0 / numPoints );

  real64 v0[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( centroid );
  LvArray::tensorOps::subtract< 3 >( v0, points[ 0 ] );
  LvArray::tensorOps::normalize< 3 >( v0 );

  // compute angles
  angle[ 0 ] = 0;
  for( localIndex a = 1; a < numPoints; ++a )
  {
    real64 v[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( centroid );
    LvArray::tensorOps::subtract< 3 >( v, points[ a ] );
    real64 const dot = LvArray::tensorOps::AiBi< 3 >( v, v0 );

    real64 crossProduct[ 3 ];
    LvArray::tensorOps::crossProduct( crossProduct, v, v0 );
    real64 const det = LvArray::tensorOps::AiBi< 3 >( normal, crossProduct );

    angle[ a ] = std::atan2( det, dot );
  }

  // sort the indices
  std::sort( indices.begin(), indices.end(), [&]( int i, int j ) { return angle[ i ] < angle[ j ]; } );

  // copy the points in the reorderedPoints array.
  for( localIndex a=0; a < numPoints; a++ )
  {
    // fill in with ordered
    LvArray::tensorOps::copy< 3 >( orderedPoints[ a ], points[ indices[ a ] ] );
  }

  for( localIndex a = 0; a < numPoints; a++ )
  {
    LvArray::tensorOps::copy< 3 >( points[a], orderedPoints[a] );
  }

  return indices;
}

/**
 * @brief Calculate the area of a polygon given the set of points in ccw order defining it.
 * @tparam NORMAL_TYPE the type of @p normal
 * @param[in] points coordinates of the points
 * @param[in] normal vector normal to the plane
 * @return the area of the polygon
 */
template< typename NORMAL_TYPE >
real64 ComputeSurfaceArea( arrayView2d< real64 const > const & points,
                           NORMAL_TYPE const && normal )
{
  real64 surfaceArea = 0.0;

  array2d< real64 > orderedPoints( points.size( 0 ), 3 );

  for( localIndex a = 0; a < points.size( 0 ); a++ )
  {
    LvArray::tensorOps::copy< 3 >( orderedPoints[a], points[a] );
  }

  orderPointsCCW( orderedPoints, normal );

  for( localIndex a = 0; a < points.size( 0 ) - 2; ++a )
  {
    real64 v1[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( orderedPoints[ a + 1 ] );
    real64 v2[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( orderedPoints[ a + 2 ] );

    LvArray::tensorOps::subtract< 3 >( v1, orderedPoints[ 0 ] );
    LvArray::tensorOps::subtract< 3 >( v2, orderedPoints[ 0 ] );

    real64 triangleNormal[ 3 ];
    LvArray::tensorOps::crossProduct( triangleNormal, v1, v2 );
    surfaceArea += LvArray::tensorOps::l2Norm< 3 >( triangleNormal );
  }

  return surfaceArea * 0.5;
}

/**
 * @brief Calculate the diameter of a set of points in a given dimension.
 * @tparam DIMENSION The dimensionality of the points.
 * @tparam POINT_COORDS_TYPE The type of the container holding the point coordinates.
 * @param[in] points The container holding the coordinates of the points.
 * @param[in] numPoints The number of points in the container.
 * @return The diameter of the set of points.
 */
template< localIndex DIMENSION, typename POINT_COORDS_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 computeDiameter( POINT_COORDS_TYPE points,
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

/**
 * @brief Calculate the centroid of a convex 3D polygon as well as the normal and the rotation matrix.
 * @tparam CENTER_TYPE The type of @p center.
 * @tparam NORMAL_TYPE The type of @p normal.
 * @param[in] pointsIndices list of index references for the points array in
 *   order (CW or CCW) about the polygon loop
 * @param[in] points 3D point list
 * @param[out] center 3D center of the given ordered polygon point list
 * @param[out] normal normal to the face
 * @param[in] areaTolerance tolerance used in the geometric computations
 * @return area of the convex 3D polygon
 * @details if area < - areaTolerance, this function will throw an error,
 *          and if (- areaTolerance <= area <= areaTolerance), the area is set to zero
 */
template< typename CENTER_TYPE, typename NORMAL_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 centroid_3DPolygon( arraySlice1d< localIndex const > const pointsIndices,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & points,
                           CENTER_TYPE && center,
                           NORMAL_TYPE && normal,
                           real64 const areaTolerance = 0.0 )
{
  real64 area = 0.0;
  LvArray::tensorOps::fill< 3 >( center, 0 );
  LvArray::tensorOps::fill< 3 >( normal, 0 );

  localIndex const numberOfPoints = pointsIndices.size();

  GEOS_ERROR_IF_LT( numberOfPoints, 2 );

  real64 current[ 3 ], next[ 3 ], crossProduct[ 3 ];

  LvArray::tensorOps::copy< 3 >( next, points[ pointsIndices[ numberOfPoints - 1 ] ] );

  for( localIndex a=0; a<numberOfPoints; ++a )
  {
    LvArray::tensorOps::copy< 3 >( current, next );
    LvArray::tensorOps::copy< 3 >( next, points[ pointsIndices[ a ] ] );

    LvArray::tensorOps::crossProduct( crossProduct, current, next );

    LvArray::tensorOps::add< 3 >( normal, crossProduct );
    LvArray::tensorOps::add< 3 >( center, next );
  }

  area = LvArray::tensorOps::l2Norm< 3 >( normal );
  LvArray::tensorOps::scale< 3 >( center, 1.0 / numberOfPoints );

  if( area > areaTolerance )
  {
    LvArray::tensorOps::normalize< 3 >( normal );
    area *= 0.5;
  }
  else if( area < -areaTolerance )
  {
    for( localIndex a=0; a<numberOfPoints; ++a )
    {
      GEOS_LOG_RANK( "Points: " << points[ pointsIndices[ a ] ] << " " << pointsIndices[ a ] );
    }
    GEOS_ERROR( "Negative area found : " << area );
  }
  else
  {
    return 0.0;
  }

  return area;
}

/**
 * @brief Change the orientation of the input vector to be consistent in a global sense.
 * @tparam NORMAL_TYPE type of @p normal
 * @param[inout] normal normal to the face
 */
template< typename NORMAL_TYPE >
GEOS_HOST_DEVICE
void FixNormalOrientation_3D( NORMAL_TYPE && normal )
{
  real64 const orientationTolerance = 10 * machinePrecision;

  // Orient local normal in global sense.
  // First check: align with z direction
  if( normal[ 2 ] <= -orientationTolerance )
  {
    LvArray::tensorOps::scale< 3 >( normal, -1.0 );
  }
  else if( std::fabs( normal[ 2 ] ) < orientationTolerance )
  {
    // If needed, second check: align with y direction
    if( normal[ 1 ] <= -orientationTolerance )
    {
      LvArray::tensorOps::scale< 3 >( normal, -1.0 );
    }
    else if( fabs( normal[ 1 ] ) < orientationTolerance )
    {
      // If needed, third check: align with x direction
      if( normal[ 0 ] <= -orientationTolerance )
      {
        LvArray::tensorOps::scale< 3 >( normal, -1.0 );
      }
    }
  }
}

/**
 * @brief Calculate the rotation matrix for a face in the 3D space
 * @tparam NORMAL_TYPE type of @p normal
 * @tparam MATRIX_TYPE type of @p rotationMatrix
 * @param[in] normal normal to the face
 * @param[out] rotationMatrix rotation matrix for the face
 */
template< typename NORMAL_TYPE, typename MATRIX_TYPE >
GEOS_HOST_DEVICE
void RotationMatrix_3D( NORMAL_TYPE const & normal,
                        MATRIX_TYPE && rotationMatrix )
{
  real64 m1[ 3 ] = { normal[ 2 ], 0.0, -normal[ 0 ] };
  real64 m2[ 3 ] = { 0.0, normal[ 2 ], -normal[ 1 ] };
  real64 const norm_m1 = LvArray::tensorOps::l2Norm< 3 >( m1 );
  real64 const norm_m2 = LvArray::tensorOps::l2Norm< 3 >( m2 );

  // If present, looks for a vector with 0 norm
  // Fix the uncertain case of norm_m1 very close to norm_m2
  if( norm_m1+1.e+2*machinePrecision > norm_m2 )
  {
    LvArray::tensorOps::crossProduct( m2, normal, m1 );
    LvArray::tensorOps::normalize< 3 >( m2 );
    LvArray::tensorOps::normalize< 3 >( m1 );
  }
  else
  {
    LvArray::tensorOps::crossProduct( m1, normal, m2 );
    LvArray::tensorOps::scale< 3 >( m1, -1 );
    LvArray::tensorOps::normalize< 3 >( m1 );
    LvArray::tensorOps::normalize< 3 >( m2 );
  }

  // Save everything in the standard form (3x3 rotation matrix)
  rotationMatrix[ 0 ][ 0 ] = normal[ 0 ];
  rotationMatrix[ 1 ][ 0 ] = normal[ 1 ];
  rotationMatrix[ 2 ][ 0 ] = normal[ 2 ];
  rotationMatrix[ 0 ][ 1 ] = m1[ 0 ];
  rotationMatrix[ 1 ][ 1 ] = m1[ 1 ];
  rotationMatrix[ 2 ][ 1 ] = m1[ 2 ];
  rotationMatrix[ 0 ][ 2 ] = m2[ 0 ];
  rotationMatrix[ 1 ][ 2 ] = m2[ 1 ];
  rotationMatrix[ 2 ][ 2 ] = m2[ 2 ];

  GEOS_ERROR_IF( fabs( LvArray::tensorOps::determinant< 3 >( rotationMatrix ) - 1.0 ) > 1.e+1 * machinePrecision,
                 "Rotation matrix with determinant different from +1.0" );
}

/**
 * @brief Return the sign of a given value as an integer.
 * @tparam T type of value
 * @param val the value in question
 * @return -1, 0 or 1 depending on whether the value is negative, zero or positive
 */
template< typename T >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
int sign( T const val )
{
  return (T( 0 ) < val) - (val < T( 0 ));
}

/**
 * @brief Check if a point is inside a convex polyhedron (3D polygon)
 * @tparam POINT_TYPE type of @p point
 * @param[in] nodeCoordinates a global array of nodal coordinates
 * @param[in] faceIndices global indices of the faces of the cell
 * @param[in] facesToNodes map from face to nodes
 * @param[in] elemCenter coordinates of the element center
 * @param[in] point coordinates of the query point
 * @param[in] areaTolerance same as in centroid_3DPolygon
 * @return whether the point is inside
 *
 * @note For faces with n>3 nodes that are non-planar, average normal is used
 */
template< typename POINT_TYPE >
GEOS_HOST_DEVICE
bool isPointInsidePolyhedron( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodeCoordinates,
                              arraySlice1d< localIndex const > const & faceIndices,
                              ArrayOfArraysView< localIndex const > const & facesToNodes,
                              POINT_TYPE const & elemCenter,
                              POINT_TYPE const & point,
                              real64 const areaTolerance = 0.0 )
{
  localIndex const numFaces = faceIndices.size();
  R1Tensor faceCenter, faceNormal, cellToFaceVec;

  for( localIndex kf = 0; kf < numFaces; ++kf )
  {
    // compute the face normal at this face
    localIndex const faceIndex = faceIndices[kf];
    centroid_3DPolygon( facesToNodes[faceIndex], nodeCoordinates, faceCenter, faceNormal, areaTolerance );

    // make sure that the normal is outward pointing
    LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
    LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter );
    if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    // compute the vector face center to query point
    LvArray::tensorOps::subtract< 3 >( faceCenter, point );
    int const s = sign( LvArray::tensorOps::AiBi< 3 >( faceNormal, faceCenter ) );

    // all dot products should be non-negative (we enforce outward normals)
    if( s < 0 )
    {
      return false;
    }
  }
  return true;
}


/**
 * @brief Method to perform lexicographic comparison of two nodes based on coordinates.
 * @tparam COORD_TYPE type of coordinate
 * @tparam POINT_TYPE type of point
 * @param[in] ax x-coordinate of the first vertex
 * @param[in] ay y-coordinate of the first vertex
 * @param[in] az z-coordinate of the first vertex
 * @param[in] bx x-coordinate of the second vertex
 * @param[in] by y-coordinate of the second vertex
 * @param[in] bz z-coordinate of the second vertex
 * @return 0 if the vertices coincide, +1 if a > b and -1 if b > a
 */
template< typename COORD_TYPE, typename POINT_TYPE >
GEOS_HOST_DEVICE
int lexicographicalCompareVertex( POINT_TYPE const ax, POINT_TYPE const ay, POINT_TYPE const az,
                                  COORD_TYPE const bx, COORD_TYPE const by, COORD_TYPE const bz )
{
  if( ax < bx )
    return -1;
  else if( ax > bx )
    return 1;
  if( ay < by )
    return -1;
  else if( ay > by )
    return 1;
  if( az < bz )
    return -1;
  else if( az > bz )
    return 1;
  return 0;
}

/**
 * @brief Method to perform lexicographic comparison of a node and an edge based on coordinates.
 * @tparam COORD_TYPE type of coordinate
 * @tparam POINT_TYPE type of point
 * @param[in] ax x-coordinate of the first vertex
 * @param[in] ay y-coordinate of the first vertex
 * @param[in] az z-coordinate of the first vertex
 * @param[in] e1x x-coordinate of the first edge vertex
 * @param[in] e1y y-coordinate of the first edge vertex
 * @param[in] e1z z-coordinate of the first edge vertex
 * @param[in] e2x x-coordinate of the second edge vertex
 * @param[in] e2y y-coordinate of the second edge vertex
 * @param[in] e2z z-coordinate of the second edge vertex
 * @return 0 if the vertices lies on the edge, +1 if e > a and -1 if a > e
 */
template< typename COORD_TYPE, typename POINT_TYPE >
GEOS_HOST_DEVICE
int lexicographicalCompareEdge( POINT_TYPE const ax, POINT_TYPE const ay, POINT_TYPE const az,
                                COORD_TYPE const e1x, COORD_TYPE const e1y, COORD_TYPE const e1z,
                                COORD_TYPE const e2x, COORD_TYPE const e2y, COORD_TYPE const e2z )
{
  return lexicographicalCompareVertex( ( e1y - ay ) * ( e2x - ax ),
                                       ( e1z - az ) * ( e2x - ax ),
                                       ( e1z - az ) * ( e2y - ay ),
                                       ( e1x - ax ) * ( e2y - ay ),
                                       ( e1x - ax ) * ( e2z - az ),
                                       ( e1y - ay ) * ( e2z - az ) );
}

/**
 * @brief Method to perform lexicographic comparison of a node and a triangle based on coordinates.
 * @tparam COORD_TYPE type of coordinate
 * @tparam POINT_TYPE type of point
 * @param[in] ax x-coordinate of the first vertex
 * @param[in] ay y-coordinate of the first vertex
 * @param[in] az z-coordinate of the first vertex
 * @param[in] t1x x-coordinate of the first triangle vertex
 * @param[in] t1y y-coordinate of the first triangle vertex
 * @param[in] t1z z-coordinate of the first triangle vertex
 * @param[in] t2x x-coordinate of the second triangle vertex
 * @param[in] t2y y-coordinate of the second triangle vertex
 * @param[in] t2z z-coordinate of the second triangle vertex
 * @param[in] t3x x-coordinate of the third triangle vertex
 * @param[in] t3y y-coordinate of the third triangle vertex
 * @param[in] t3z z-coordinate of the third triangle vertex
 * @return 0 if the vertices lies on the triangle, +1 if t > a and -1 if t > e
 */
template< typename COORD_TYPE, typename POINT_TYPE >
GEOS_HOST_DEVICE
int lexicographicalCompareTriangle( POINT_TYPE const ax, POINT_TYPE const ay, POINT_TYPE const az,
                                    COORD_TYPE const t1x, COORD_TYPE const t1y, COORD_TYPE const t1z,
                                    COORD_TYPE const t2x, COORD_TYPE const t2y, COORD_TYPE const t2z,
                                    COORD_TYPE const t3x, COORD_TYPE const t3y, COORD_TYPE const t3z )
{
  COORD_TYPE v1x = t1x - ax;
  COORD_TYPE v1y = t1y - ay;
  COORD_TYPE v1z = t1z - az;
  COORD_TYPE v2x = t2x - ax;
  COORD_TYPE v2y = t2y - ay;
  COORD_TYPE v2z = t2z - az;
  COORD_TYPE v3x = t3x - ax;
  COORD_TYPE v3y = t3y - ay;
  COORD_TYPE v3z = t3z - az;
  COORD_TYPE sign = ( v1x * v2y - v1y * v2x ) * v3z +
                    ( v2x * v3y - v2y * v3x ) * v1z +
                    ( v3x * v1y - v3y * v1x ) * v2z;
  if( sign > 0 )
    return 1;
  else if( sign < 0 )
    return -1;
  return 0;
}
/**
 * @brief Method to find the reference element touching a vertex. The element with the lowest global ID is chosen
 * from the list.
 * @param[in] nodeElements the list of elements adjacent to the vertex
 * @param[in] elementGlobalIndex the global IDs for elements
 * @return element touching the vertex with the least global index
 */
template< typename ... LIST_TYPE >
GEOS_HOST_DEVICE
int findVertexRefElement( arraySlice1d< localIndex const > const & nodeElements,
                          arrayView1d< globalIndex const > const & elementGlobalIndex )
{
  localIndex minElement = -1;
  globalIndex minElementGID = LvArray::NumericLimits< globalIndex >::max;
  for( int i = 0; i < nodeElements.size(); i++ )
  {
    localIndex e = nodeElements( i );
    if( elementGlobalIndex[ e ] < minElementGID )
    {
      minElementGID = elementGlobalIndex[ e ];
      minElement = e;
    }
  }
  return minElement;
}

/**
 * @brief Method to find the reference element for an edge. The element with the lowest global ID is chosen
 * from the list.
 * @param[in] nodeElements1 the list of elements adjacent to the first node
 * @param[in] nodeElements2 the list of elements adjacent to the second node
 * @param[in] elementGlobalIndex the global IDs for elements
 * @return the element shared by the two nodes, with the minimal global index
 */
template< typename ... LIST_TYPE >
GEOS_HOST_DEVICE
int findEdgeRefElement( arraySlice1d< localIndex const > const & nodeElements1,
                        arraySlice1d< localIndex const > const & nodeElements2,
                        arrayView1d< globalIndex const > const & elementGlobalIndex )
{
  localIndex minElement = -1;
  globalIndex minElementGID = LvArray::NumericLimits< globalIndex >::max;
  for( int i = 0; i < nodeElements1.size(); i++ )
  {
    localIndex e1 = nodeElements1( i );
    for( int j = 0; j < nodeElements2.size(); j++ )
    {
      localIndex e2 = nodeElements2( j );
      if( e1 == e2 )
      {
        if( elementGlobalIndex[ e1 ] < minElementGID )
        {
          minElementGID = elementGlobalIndex[ e1 ];
          minElement = e1;
        }
      }
    }
  }
  return minElement;
}

/**
 * @brief Method to find the reference element for a triangle. The element with the lowest global ID is chosen
 * from the list.
 * @param[in] nodeElements1 the list of elements adjacent to the first node
 * @param[in] nodeElements2 the list of elements adjacent to the second node
 * @param[in] nodeElements3 the list of elements adjacent to the third node
 * @param[in] elementGlobalIndex the global IDs for elements
 * @return the element shared by the three nodes, with the minimal global index
 */
template< typename ... LIST_TYPE >
GEOS_HOST_DEVICE
int findTriangleRefElement( arraySlice1d< localIndex const > const & nodeElements1,
                            arraySlice1d< localIndex const > const & nodeElements2,
                            arraySlice1d< localIndex const > const & nodeElements3,
                            arrayView1d< globalIndex const > const & elementGlobalIndex )
{
  localIndex minElement = -1;
  globalIndex minElementGID = LvArray::NumericLimits< globalIndex >::max;
  for( int i = 0; i < nodeElements1.size(); i++ )
  {
    localIndex e1 = nodeElements1( i );
    for( int j = 0; j < nodeElements2.size(); j++ )
    {
      localIndex e2 = nodeElements2( j );
      for( int k = 0; k < nodeElements3.size(); k++ )
      {
        localIndex e3 = nodeElements3( k );
        if( e1 == e2 && e2 == e3 )
        {
          if( elementGlobalIndex[ e1 ] < minElementGID )
          {
            minElementGID = elementGlobalIndex[ e1 ];
            minElement = e1;
          }
        }
      }
    }
  }
  return minElement;
}

/**
 * @brief Computes the winding number of a point with respecto to a mesh element.
 * @tparam POINT_TYPE type of @p point
 * @param[in] element the element to be checked
 * @param[in] nodeCoordinates a global array of nodal coordinates
 * @param[in] elementsToFaces map from elements to faces
 * @param[in] facesToNodes map from faces to nodes
 * @param[in] nodesToElements map from nodes to elements
 * @param[in] nodeLocalToGlobal global indices of nodes
 * @param[in] elementLocalToGlobal global indices of elements
 * @param[in] elemCenter coordinates of the element centroid
 * @param[in] point coordinates of the query point
 * @return the signed winding number, which is positive if and only if the point is inside the mesh element.
 */
template< typename COORD_TYPE, typename POINT_TYPE >
GEOS_HOST_DEVICE
bool computeWindingNumber( localIndex element,
                           arrayView2d< COORD_TYPE const, nodes::REFERENCE_POSITION_USD > const & nodeCoordinates,
                           arrayView2d< localIndex const > const & elementsToFaces,
                           ArrayOfArraysView< localIndex const > const & facesToNodes,
                           ArrayOfArraysView< localIndex const > const & nodesToElements,
                           arrayView1d< globalIndex const > const & nodeLocalToGlobal,
                           arrayView1d< globalIndex const > const & elementLocalToGlobal,
                           POINT_TYPE const & elemCenter,
                           POINT_TYPE const & point )
{
  arraySlice1d< localIndex const > const & faceIndices = elementsToFaces[ element ];
  localIndex const numFaces = faceIndices.size();
  int omega = 0;
  for( localIndex kf = 0; kf < numFaces; ++kf )
  {
    // triangulate the face. The triangulation must be done in a consistent way across ranks.
    // This can be achieved by always picking the vertex with the lowest global index as root.
    localIndex const faceIndex = faceIndices[kf];
    globalIndex minGlobalId = LvArray::NumericLimits< globalIndex >::max;
    localIndex minVertex = -1;
    localIndex numFaceVertices = facesToNodes[faceIndex].size();
    for( localIndex v = 0; v < numFaceVertices; v++ )
    {
      localIndex vIndex = facesToNodes( faceIndex, v );
      globalIndex globalId = nodeLocalToGlobal[ vIndex ];
      if( globalId < minGlobalId )
      {
        minGlobalId = globalId;
        minVertex = vIndex;
      }
    }
    // triangulate the face using the minimum-id vertex as root
    localIndex vi[ 3 ] = { minVertex, -1, -1 };
    for( localIndex v = 0; v < numFaceVertices; v++ )
    {
      vi[ 1 ] = facesToNodes( faceIndex, v );
      vi[ 2 ] = facesToNodes( faceIndex, (v + 1) % numFaceVertices );
      if( vi[ 1 ] != minVertex && vi[ 2 ] != minVertex )
      {
        // To make the algorithm independent of rank, always take the two additional vertices in increasing global ID
        if( nodeLocalToGlobal[ vi[ 1 ] ] > nodeLocalToGlobal[ vi[ 2 ] ] )
        {
          localIndex temp = vi[ 1 ];
          vi[ 1 ] = vi[ 2 ];
          vi[ 2 ] = temp;
        }
        COORD_TYPE v1x = nodeCoordinates( vi[ 0 ], 0 );
        COORD_TYPE v1y = nodeCoordinates( vi[ 0 ], 1 );
        COORD_TYPE v1z = nodeCoordinates( vi[ 0 ], 2 );
        COORD_TYPE v2x = nodeCoordinates( vi[ 1 ], 0 );
        COORD_TYPE v2y = nodeCoordinates( vi[ 1 ], 1 );
        COORD_TYPE v2z = nodeCoordinates( vi[ 1 ], 2 );
        COORD_TYPE v3x = nodeCoordinates( vi[ 2 ], 0 );
        COORD_TYPE v3y = nodeCoordinates( vi[ 2 ], 1 );
        COORD_TYPE v3z = nodeCoordinates( vi[ 2 ], 2 );
        // check the orientation of this triangle
        R1Tensor vv1 = { v2x - v1x, v2y - v1y, v2z - v1z  };
        R1Tensor vv2 = { v3x - v1x, v3y - v1y, v3z - v1z  };
        R1Tensor dist = { elemCenter[ 0 ] - ( v1x + v2x + v3x )/3.0,
                          elemCenter[ 1 ] - ( v1y + v2y + v3y )/3.0,
                          elemCenter[ 2 ] - ( v1z + v2z + v3z )/3.0 };
        R1Tensor norm = { };
        LvArray::tensorOps::crossProduct( norm, vv1, vv2 );
        // check if face is oriented coherently, and change sign otherwise
        int sign = LvArray::tensorOps::AiBi< 3 >( norm, dist ) > 0 ? -1 : +1;
        // Compute the winding number contributed by this triangle
        int cmp1 = lexicographicalCompareVertex( point[ 0 ], point[ 1 ], point[ 2 ], v1x, v1y, v1z );
        if( cmp1 == 0 )
        {
          return findVertexRefElement( nodesToElements[ vi[ 0 ] ], elementLocalToGlobal ) == element;
        }
        int cmp2 = lexicographicalCompareVertex( point[ 0 ], point[ 1 ], point[ 2 ], v2x, v2y, v2z );
        if( cmp2 == 0 )
        {
          return findVertexRefElement( nodesToElements[ vi[ 1 ] ], elementLocalToGlobal ) == element;
        }
        int cmp3 = lexicographicalCompareVertex( point[ 0 ], point[ 1 ], point[ 2 ], v3x, v3y, v3z );
        if( cmp3 == 0 )
        {
          return findVertexRefElement( nodesToElements[ vi[ 2 ] ], elementLocalToGlobal ) == element;
        }
        int facecmp = 0;
        int edgecmp = 0;
        if( cmp1 != cmp2 )
        {
          edgecmp = lexicographicalCompareEdge( point[ 0 ], point[ 1 ], point[ 2 ],
                                                v1x, v1y, v1z,
                                                v2x, v2y, v2z );
          if( edgecmp == 0 )
          {
            return findEdgeRefElement( nodesToElements[ vi[ 0 ] ], nodesToElements[ vi[ 1 ] ], elementLocalToGlobal ) == element;
          }
          facecmp += sign * edgecmp;
        }
        if( cmp2 != cmp3 )
        {
          edgecmp = lexicographicalCompareEdge( point[ 0 ], point[ 1 ], point[ 2 ],
                                                v2x, v2y, v2z,
                                                v3x, v3y, v3z );
          if( edgecmp == 0 )
          {
            return findEdgeRefElement( nodesToElements[ vi[ 1 ] ], nodesToElements[ vi[ 2 ] ], elementLocalToGlobal ) == element;
          }
          facecmp += sign * edgecmp;
        }
        if( cmp3 != cmp1 )
        {
          edgecmp = lexicographicalCompareEdge( point[ 0 ], point[ 1 ], point[ 2 ],
                                                v3x, v3y, v3z,
                                                v1x, v1y, v1z );
          if( edgecmp == 0 )
          {
            return findEdgeRefElement( nodesToElements[ vi[ 0 ] ], nodesToElements[ vi[ 2 ] ], elementLocalToGlobal ) == element;
          }
          facecmp += sign * edgecmp;
        }
        // if all edges are on the same side, this triangle does not contribute to the winding number
        if( facecmp == 0 )
          continue;
        facecmp = lexicographicalCompareTriangle( point[ 0 ], point[ 1 ], point[ 2 ],
                                                  v1x, v1y, v1z,
                                                  v2x, v2y, v2z,
                                                  v3x, v3y, v3z );

        if( facecmp == 0 )
        {
          return findTriangleRefElement( nodesToElements[ vi[ 0 ] ], nodesToElements[ vi[ 1 ] ], nodesToElements[ vi[ 2 ] ], elementLocalToGlobal ) == element;
        }
        omega += sign * facecmp;
      }
    }
  }

  return omega;
}

/**
 * @brief Check if a point is inside a convex polyhedron (3D polygon), using a robust method
 *  to avoid ambiguity when the point lies on an interface.
 *  This method is based on the following method:
 *  - the winding number omega of the point with respect to the cell is used to determine if the point is inside (omega=1) or not (omega=0)
 *  - corner cases (point lying on a face, edge or vertex of the cell) are detected using a robust method based on lexicographical
 * comparisons
 *  - these comparisons are made consistent across MPI ranks by consistently arranging items based on global indices (GIDs). In particular:
 *    - Faces are triangulated using the vertex with the smallest GID as root;
 *    - Edges and faces are described by vertices in increasing GID order
 *  - Finally, if the point lies on a (vertex, edge, face), it is assigned to the first shared element with the least global index
 * @tparam POINT_TYPE type of @p point
 * @param[in] element the element to be checked
 * @param[in] nodeCoordinates a global array of nodal coordinates
 * @param[in] elementsToFaces map from elements to faces
 * @param[in] facesToNodes map from faces to nodes
 * @param[in] nodesToElements map from nodes to elements
 * @param[in] nodeLocalToGlobal global indices of nodes
 * @param[in] elementLocalToGlobal global indices of elements
 * @param[in] elemCenter coordinates of the element centroid
 * @param[in] point coordinates of the query point
 * @return whether the point is inside
 */
template< typename COORD_TYPE, typename POINT_TYPE >
GEOS_HOST_DEVICE
bool isPointInsideConvexPolyhedronRobust( localIndex element,
                                          arrayView2d< COORD_TYPE const, nodes::REFERENCE_POSITION_USD > const & nodeCoordinates,
                                          arrayView2d< localIndex const > const & elementsToFaces,
                                          ArrayOfArraysView< localIndex const > const & facesToNodes,
                                          ArrayOfArraysView< localIndex const > const & nodesToElements,
                                          arrayView1d< globalIndex const > const & nodeLocalToGlobal,
                                          arrayView1d< globalIndex const > const & elementLocalToGlobal,
                                          POINT_TYPE const & elemCenter,
                                          POINT_TYPE const & point )
{
  return computeWindingNumber( element, nodeCoordinates, elementsToFaces, facesToNodes, nodesToElements, nodeLocalToGlobal, elementLocalToGlobal, elemCenter, point ) > 0;
}

/**
 * @brief Compute the dimensions of the bounding box containing the element
 *   defined here by the coordinates of its vertices.
 * @tparam VEC_TYPE type of @p boxDims
 * @param[in] elemIndex index of the element in pointIndices.
 * @param[in] pointIndices the indices of the vertices in pointCoordinates.
 * @param[in] pointCoordinates the vertices coordinates.
 * @param[out] boxDims The dimensions of the bounding box.
 */
template< typename VEC_TYPE >
GEOS_HOST_DEVICE
void getBoundingBox( localIndex const elemIndex,
                     arrayView2d< localIndex const, cells::NODE_MAP_USD > const & pointIndices,
                     arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & pointCoordinates,
                     VEC_TYPE && boxDims )
{
  // This holds the min coordinates of the set in each direction
  R1Tensor minCoords = { LvArray::NumericLimits< real64 >::max,
                         LvArray::NumericLimits< real64 >::max,
                         LvArray::NumericLimits< real64 >::max };

  // boxDims is used to hold the max coordinates.
  LvArray::tensorOps::fill< 3 >( boxDims, LvArray::NumericLimits< real64 >::lowest );

  // loop over all the vertices of the element to get the min and max coords
  for( localIndex a = 0; a < pointIndices.size( 1 ); ++a )
  {
    localIndex const id = pointIndices( elemIndex, a );
    for( localIndex d = 0; d < 3; ++d )
    {
      minCoords[ d ] = fmin( minCoords[ d ], pointCoordinates( id, d ) );
      boxDims[ d ] = fmax( boxDims[ d ], pointCoordinates( id, d ) );
    }
  }

  LvArray::tensorOps::subtract< 3 >( boxDims, minCoords );
}

/**
 * @brief Compute the volume of an element (tetrahedron, pyramid, wedge, hexahedron)
 * @tparam FE_TYPE the type of finite element space
 * @param[in] X vertices of the element
 * @return the volume of the element
 */
template< typename FE_TYPE >
GEOS_HOST_DEVICE inline
real64 elementVolume( real64 const (&X)[FE_TYPE::numNodes][3] )
{
  real64 result{};
  for( localIndex q=0; q<FE_TYPE::numQuadraturePoints; ++q )
  {
    result = result + FE_TYPE::transformedQuadratureWeight( q, X );
  }
  return result;
}

/**
 * @brief Compute the volume of an hexahedron
 * @param[in] X vertices of the hexahedron
 * @return the volume of the hexahedron
 */
GEOS_HOST_DEVICE
inline
real64 hexahedronVolume( real64 const (&X)[8][3] )
{
  return elementVolume< finiteElement::H1_Hexahedron_Lagrange1_GaussLegendre2 >( X );
}

/**
 * @brief Compute the volume of an tetrahedron
 * @param[in] X vertices of the tetrahedron
 * @return the volume of the tetrahedron
 */
GEOS_HOST_DEVICE
inline
real64 tetrahedronVolume( real64 const (&X)[4][3] )
{
  return elementVolume< finiteElement::H1_Tetrahedron_Lagrange1_Gauss1 >( X );
}

/**
 * @brief Compute the volume of a wedge
 * @param[in] X vertices of the wedge
 * @return the volume of the wedge
 */
GEOS_HOST_DEVICE
inline
real64 wedgeVolume( real64 const (&X)[6][3] )
{
  return elementVolume< finiteElement::H1_Wedge_Lagrange1_Gauss6 >( X );
}

/**
 * @brief Compute the volume of a pyramid
 * @param[in] X vertices of the pyramid
 * @return the volume of the pyramid
 */
GEOS_HOST_DEVICE
inline
real64 pyramidVolume( real64 const (&X)[5][3] )
{
  return elementVolume< finiteElement::H1_Pyramid_Lagrange1_Gauss5 >( X );
}

/**
 * @brief Compute the volume of a prism with N-sided polygon base
 * @tparam N the number of sides in the polygon base
 * @param[in] X vertices of the prism
 * @return the volume of the prism
 *
 * @note The volume is computed splitting the prism into wedges. The function can be called only
 *       for N > 5. For N = 3 and N = 4 function wedgeVolume and hexahedronVolume, respectively,
 *       should be used.
 */
template< integer N >
GEOS_HOST_DEVICE
inline
real64 prismVolume( real64 const (&X)[2*N][3] )
{
  static_assert( N > 4,
                 "Function prismVolume can be called for a prism with N-sided polygon base where N > 5." );

  real64 result{};

  // Compute the barycenters of the prism bases
  real64 XGBot[3]{};
  real64 XGTop[3]{};
  for( integer a = 0; a < N; ++a )
  {
    LvArray::tensorOps::add< 3 >( XGBot, X[a] );
  }
  for( integer a = N; a < 2 * N; ++a )
  {
    LvArray::tensorOps::add< 3 >( XGTop, X[a] );
  }
  LvArray::tensorOps::scale< 3 >( XGBot, 1.0 / N );
  LvArray::tensorOps::scale< 3 >( XGTop, 1.0 / N );

  real64 XWedge[6][3];
  for( int a = 0; a < N - 1; ++a )
  {

    LvArray::tensorOps::copy< 3 >( XWedge[0], X[a] );
    LvArray::tensorOps::copy< 3 >( XWedge[1], X[a+N] );
    LvArray::tensorOps::copy< 3 >( XWedge[2], X[a+1] );
    LvArray::tensorOps::copy< 3 >( XWedge[3], X[a+1+N] );
    LvArray::tensorOps::copy< 3 >( XWedge[4], XGBot );
    LvArray::tensorOps::copy< 3 >( XWedge[5], XGTop );
    result = result + computationalGeometry::elementVolume< finiteElement::H1_Wedge_Lagrange1_Gauss6 >( XWedge );
  }
  LvArray::tensorOps::copy< 3 >( XWedge[0], X[N-1] );
  LvArray::tensorOps::copy< 3 >( XWedge[1], X[2*N-1] );
  LvArray::tensorOps::copy< 3 >( XWedge[2], X[0] );
  LvArray::tensorOps::copy< 3 >( XWedge[3], X[N] );
  LvArray::tensorOps::copy< 3 >( XWedge[4], XGBot );
  LvArray::tensorOps::copy< 3 >( XWedge[5], XGTop );
  result = result + computationalGeometry::elementVolume< finiteElement::H1_Wedge_Lagrange1_Gauss6 >( XWedge );
  return result;
}

} /* namespace computationalGeometry */
} /* namespace geos */

#endif /* GEOS_MESH_UTILITIES_COMPUTATIONALGEOMETRY_HPP_ */
