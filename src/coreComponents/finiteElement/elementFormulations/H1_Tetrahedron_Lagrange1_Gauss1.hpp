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
 * @file H1_Tetrahedron_Lagrange1_Gauss1.hpp
 */

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE1GAUSS1_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE1GAUSS1_HPP_

#include "FiniteElementBase.hpp"
#include "LagrangeBasis1.hpp"


namespace geos
{
namespace finiteElement
{

/**
 * This class contains the kernel accessible functions specific to the
 * H1-conforming nodal linear tetrahedron finite element with a 1-point
 * Gaussian quadrature rule.
 *
 *          3                                =====  ==  ==  ==
 *           +                               Node   r   s   t
 *           |\\_                            =====  ==  ==  ==
 *           ||  \_                          0      0   0   0
 *           | \   \_           t            1      1   0   0
 *           |  +__  \_         |   s        2      0   1   0
 *           | /2  \__ \_       |  /         3      0   0   1
 *           |/       \__\      | /          =====  ==  ==  ==
 *           +------------+     *------r
 *          0              1
 *
 */
class H1_Tetrahedron_Lagrange1_Gauss1 final : public FiniteElementBase
{
public:

  /// The type of basis used for this element
  using BASIS = LagrangeBasis1;

  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 4;

  /// The number of faces/support points per element.
  constexpr static localIndex numFaces = 4;

  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 14;
  //constexpr static localIndex numQuadraturePoints = 5;
  //constexpr static localIndex numQuadraturePoints = 4;
  //constexpr static localIndex numQuadraturePoints = 1;

  /// The number of sampling points per element.
  constexpr static int numSamplingPoints = numSamplingPointsPerDirection * numSamplingPointsPerDirection * numSamplingPointsPerDirection;

  virtual ~H1_Tetrahedron_Lagrange1_Gauss1() override
  {}

  GEOS_HOST_DEVICE
  virtual localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  /**
   * @brief Get the number of quadrature points.
   * @param stack Stack variables as filled by @ref setupStack.
   * @return The number of quadrature points.
   */
  GEOS_HOST_DEVICE
  static localIndex
  getNumQuadraturePoints( StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numQuadraturePoints;
  }

  GEOS_HOST_DEVICE
  virtual localIndex getNumSupportPoints() const override
  {
    return numNodes;
  }

  GEOS_HOST_DEVICE
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
  static localIndex getNumSupportPoints( StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numNodes;
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
    int const i0 = linearIndex % numSamplingPointsPerDirection;
    int const i1 = ( (linearIndex - i0)/numSamplingPointsPerDirection ) % numSamplingPointsPerDirection;
    int const i2 = ( (linearIndex - i0)/numSamplingPointsPerDirection - i1 ) / numSamplingPointsPerDirection;

    real64 const step = 1. / ( numSamplingPointsPerDirection - 1 );

    real64 const r = i0 * step;
    real64 const s = i1 * step;
    real64 const t = i2 * step;
    if( (r+s) <= t )
    {
      samplingPointCoord[0] = r;
      samplingPointCoord[1] = s;
      samplingPointCoord[2] = t;
    }
    else
    {
      // if outside of the triangle need to reproject it. Points will be doubled though.
      samplingPointCoord[0] = 1 - t - r;
      samplingPointCoord[1] = 1 - t - s;
      samplingPointCoord[2] = t;
    }
  }

  /**
   * @brief Calculate shape functions values for each support point at a
   *   given point in the parent space.
   * @param pointCoord coordinates of the given point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const (&pointCoord)[3],
                     real64 ( &N )[numNodes] );

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOS_HOST_DEVICE
  static void calcN( localIndex const q,
                     real64 ( &N )[numNodes] );

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
                     real64 ( &N )[numNodes] );

  /**
   * @brief Calculate face bubble functions values for each face at a
   *   given point in the parent space.
   * @param pointCoord coordinates of the given point.
   * @param N An array to pass back the shape function values for each support
   *   face.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcFaceBubbleN( real64 const (&pointCoord)[3],
                               real64 (& N)[numFaces] )
  {

    real64 const r = pointCoord[0];
    real64 const s = pointCoord[1];
    real64 const t = pointCoord[2];

    N[0] = (1 - r - s - t) * r * t;
    N[1] = (1 - r - s - t) * r * s;
    N[2] = (1 - r - s - t) * t * s;
    N[3] = r * s * t;
  }

  /**
   * @brief Calculate face bubble functions values for each face at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOS_HOST_DEVICE
  inline
  static void calcFaceBubbleN( localIndex const q,
                               real64 (& N)[numFaces] )
  {

    real64 const pointCoord[3] = {quadratureParentCoords0( q ),
                                  quadratureParentCoords1( q ),
                                  quadratureParentCoords2( q )};

    calcFaceBubbleN( pointCoord, N );
  }

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[numNodes][3],
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  inline
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[numNodes][3],
                           StackVariables const & stack,
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the bubble function derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param gradN Array to contain the shape bubble function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  static real64 calcGradFaceBubbleN( localIndex const q,
                                     real64 const (&X)[numNodes][3],
                                     real64 ( &gradN )[numFaces][3] );

  /**
   * @brief Calculate the integration weights for a quadrature point.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @return The product of the quadrature rule weight and the determinate of
   *   the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  static real64 transformedQuadratureWeight( localIndex const q,
                                             real64 const (&X)[numNodes][3] );

  /**
   * @brief Calculate the integration weights for a quadrature point.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @return The product of the quadrature rule weight and the determinate of
   *   the parent/physical transformation matrix.
   */
  GEOS_HOST_DEVICE
  static real64 transformedQuadratureWeight( localIndex const q,
                                             real64 const (&X)[numNodes][3],
                                             StackVariables const & stack )
  { GEOS_UNUSED_VAR( stack ); return transformedQuadratureWeight( q, X ); }

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param q The quadrature point index in 3d space.
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   * @return The determinant of the Jacobian transformation matrix.
   */
  GEOS_HOST_DEVICE
  static real64 invJacobianTransformation( int const q,
                                           real64 const (&X)[numNodes][3],
                                           real64 ( & J )[3][3] )
  {
    GEOS_UNUSED_VAR( q, X );
//    jacobianTransformation( q, X, J );
    return LvArray::tensorOps::invert< 3 >( J );
  }
private:
  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = 1.0 / 6.0;

  /// The weight of each quadrature point.
  constexpr static real64 weight = parentVolume;
  //constexpr static real64 weight = parentVolume / numQuadraturePoints;

  /**
   * @brief Calculate the weights
   * @param q
   * @return weight.
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 quadratureWeight( localIndex const q )
  {

    //real64 const w[numQuadraturePoints] = {-4.0/5.0, 9.0/20.0, 9.0/20.0, 9.0/20.0, 9.0/20.0 };
    //real64 const w[numQuadraturePoints] = {-4.0, 9.0/4.0, 9.0/4.0, 9.0/4.0, 9.0/4.0 };

    
    real64 const w[numQuadraturePoints] = { 0.073493043116361949544,
                                            0.073493043116361949544,
                                            0.073493043116361949544,
                                            0.073493043116361949544,
                                            0.11268792571801585080, 
                                            0.11268792571801585080, 
                                            0.11268792571801585080, 
                                            0.11268792571801585080, 
                                            0.042546020777081466438,
                                            0.042546020777081466438,
                                            0.042546020777081466438,
                                            0.042546020777081466438,
                                            0.042546020777081466438,
                                            0.042546020777081466438 };

    return w[q];

  }

  /**
   * @brief Calculate the parent coordinates for the r direction, given the
   *        linear index of a quadrature point.
   * @param a The linear index of quadrature point
   * @return parent coordinate in the r direction.
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 quadratureParentCoords0( localIndex const q )
  {

    //real64 const qCoords[numQuadraturePoints] = { 1.0/4.0, 1.0/2.0, 1.0/6.0, 1.0/6.0, 1.0/6.0 };

    real64 const qCoords[numQuadraturePoints] = { 0.72179424906732632079,
                                                  0.092735250310891226402,
                                                  0.092735250310891226402,
                                                  0.092735250310891226402,
                                                  0.067342242210098170608,
                                                  0.31088591926330060980,
                                                  0.31088591926330060980,
                                                  0.31088591926330060980,
                                                  0.045503704125649649492,
                                                  0.045503704125649649492,
                                                  0.045503704125649649492,
                                                  0.45449629587435035051,
                                                  0.45449629587435035051,
                                                  0.45449629587435035051 };


    return qCoords[q];
  }

  /**
   * @brief Calculate the parent coordinates for the s direction, given the
   *        linear index of a quadrature point.
   * @param q The linear index of quadrature point
   * @return parent coordinate in the s direction.
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 quadratureParentCoords1( localIndex const q )
  {

    //real64 const qCoords[numQuadraturePoints] = { 1.0/4.0, 1.0/6.0, 1.0/2.0,  1.0/6.0, 1.0/6.0 };

    real64 const qCoords[numQuadraturePoints] = { 0.092735250310891226402, 
                                                  0.72179424906732632079,
                                                  0.092735250310891226402, 
                                                  0.092735250310891226402, 
                                                  0.31088591926330060980,
                                                  0.067342242210098170608,
                                                  0.31088591926330060980,
                                                  0.31088591926330060980,
                                                  0.045503704125649649492, 
                                                  0.45449629587435035051, 
                                                  0.45449629587435035051, 
                                                  0.045503704125649649492, 
                                                  0.045503704125649649492, 
                                                  0.45449629587435035051 };

    return qCoords[q];
  }

  /**
   * @brief Calculate the parent coordinates for the t direction, given the
   *        linear index of a quadrature point.
   * @param q The linear index of quadrature point
   * @return parent coordinate in the t direction.
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 quadratureParentCoords2( localIndex const q )
  {

    //real64 const qCoords[numQuadraturePoints] = { 1.0/4.0, 1.0/6.0, 1.0/6.0, 1.0/2.0, 1.0/6.0 };

    real64 const qCoords[numQuadraturePoints] = { 0.092735250310891226402,
                                                  0.092735250310891226402,
                                                  0.72179424906732632079,
                                                  0.092735250310891226402,
                                                  0.31088591926330060980,
                                                  0.31088591926330060980,
                                                  0.067342242210098170608,
                                                  0.31088591926330060980,
                                                  0.45449629587435035051,
                                                  0.045503704125649649492,
                                                  0.45449629587435035051,
                                                  0.045503704125649649492, 
                                                  0.45449629587435035051,
                                                  0.045503704125649649492 }; 
    return qCoords[q];
  }

  /**
   * @brief Calculates the determinant of the Jacobian of the isoparametric
   *        mapping from the parent space to the physical space.
   * @param X Array containing the coordinates of the support points
   * @return determinant value
   */
  GEOS_HOST_DEVICE
  static real64 determinantJacobianTransformation( real64 const (&X)[numNodes][3] );

};

/// @cond Doxygen_Suppress

GEOS_HOST_DEVICE
inline
real64
H1_Tetrahedron_Lagrange1_Gauss1::
  determinantJacobianTransformation( real64 const (&X)[numNodes][3] )
{
  return ( X[1][0] - X[0][0] )*( ( X[2][1] - X[0][1] )*( X[3][2] - X[0][2] ) - ( X[3][1] - X[0][1] )*( X[2][2] - X[0][2] ) )
         + ( X[1][1] - X[0][1] )*( ( X[3][0] - X[0][0] )*( X[2][2] - X[0][2] ) - ( X[2][0] - X[0][0] )*( X[3][2] - X[0][2] ) )
         + ( X[1][2] - X[0][2] )*( ( X[2][0] - X[0][0] )*( X[3][1] - X[0][1] ) - ( X[3][0] - X[0][0] )*( X[2][1] - X[0][1] ) );
}

//*************************************************************************************************

GEOS_HOST_DEVICE
inline
void
H1_Tetrahedron_Lagrange1_Gauss1::
  calcN( real64 const (&coords)[3],
         real64 (& N)[numNodes] )
{
  real64 const r = coords[0];
  real64 const s = coords[1];
  real64 const t = coords[2];

  // single quadrature point (centroid), i.e.  r = s = t = 1/4
  N[0] = 1 - r - s - t;
  N[1] = r;
  N[2] = s;
  N[3] = t;
}


GEOS_HOST_DEVICE
inline
void
H1_Tetrahedron_Lagrange1_Gauss1::
  calcN( localIndex const q,
         real64 (& N)[numNodes] )
{
  real64 const pointCoord[3] = {quadratureParentCoords0( q ),
                                quadratureParentCoords1( q ),
                                quadratureParentCoords2( q )};

  calcN( pointCoord, N );
}

GEOS_HOST_DEVICE
inline
void H1_Tetrahedron_Lagrange1_Gauss1::
  calcN( localIndex const q,
         StackVariables const & GEOS_UNUSED_PARAM( stack ),
         real64 ( & N )[numNodes] )
{
  return calcN( q, N );
}

//*************************************************************************************************

GEOS_HOST_DEVICE
inline
real64
H1_Tetrahedron_Lagrange1_Gauss1::
  calcGradN( localIndex const q,
             real64 const (&X)[numNodes][3],
             real64 (& gradN)[numNodes][3] )
{

  gradN[0][0] =  X[1][1]*( X[3][2] - X[2][2] ) - X[2][1]*( X[3][2] - X[1][2] ) + X[3][1]*( X[2][2] - X[1][2] );
  gradN[0][1] = -X[1][0]*( X[3][2] - X[2][2] ) + X[2][0]*( X[3][2] - X[1][2] ) - X[3][0]*( X[2][2] - X[1][2] );
  gradN[0][2] =  X[1][0]*( X[3][1] - X[2][1] ) - X[2][0]*( X[3][1] - X[1][1] ) + X[3][0]*( X[2][1] - X[1][1] );

  gradN[1][0] = -X[0][1]*( X[3][2] - X[2][2] ) + X[2][1]*( X[3][2] - X[0][2] ) - X[3][1]*( X[2][2] - X[0][2] );
  gradN[1][1] =  X[0][0]*( X[3][2] - X[2][2] ) - X[2][0]*( X[3][2] - X[0][2] ) + X[3][0]*( X[2][2] - X[0][2] );
  gradN[1][2] = -X[0][0]*( X[3][1] - X[2][1] ) + X[2][0]*( X[3][1] - X[0][1] ) - X[3][0]*( X[2][1] - X[0][1] );

  gradN[2][0] =  X[0][1]*( X[3][2] - X[1][2] ) - X[1][1]*( X[3][2] - X[0][2] ) + X[3][1]*( X[1][2] - X[0][2] );
  gradN[2][1] = -X[0][0]*( X[3][2] - X[1][2] ) + X[1][0]*( X[3][2] - X[0][2] ) - X[3][0]*( X[1][2] - X[0][2] );
  gradN[2][2] =  X[0][0]*( X[3][1] - X[1][1] ) - X[1][0]*( X[3][1] - X[0][1] ) + X[3][0]*( X[1][1] - X[0][1] );

  gradN[3][0] = -X[0][1]*( X[2][2] - X[1][2] ) + X[1][1]*( X[2][2] - X[0][2] ) - X[2][1]*( X[1][2] - X[0][2] );
  gradN[3][1] =  X[0][0]*( X[2][2] - X[1][2] ) - X[1][0]*( X[2][2] - X[0][2] ) + X[2][0]*( X[1][2] - X[0][2] );
  gradN[3][2] = -X[0][0]*( X[2][1] - X[1][1] ) + X[1][0]*( X[2][1] - X[0][1] ) - X[2][0]*( X[1][1] - X[0][1] );

  real64 detJ = determinantJacobianTransformation( X );
  real64 factor = 1.0 / ( detJ );

  for( int i = 0; i < numNodes; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      gradN[i][j] *= factor;
    }
  }

  return detJ * weight * quadratureWeight(q);
}

GEOS_HOST_DEVICE
inline
real64 H1_Tetrahedron_Lagrange1_Gauss1::
  calcGradN( localIndex const q,
             real64 const (&X)[numNodes][3],
             StackVariables const & GEOS_UNUSED_PARAM( stack ),
             real64 ( & gradN )[numNodes][3] )
{
  return calcGradN( q, X, gradN );
}

GEOS_HOST_DEVICE
inline
real64
H1_Tetrahedron_Lagrange1_Gauss1::calcGradFaceBubbleN( localIndex const q,
                                                      real64 const (&X)[numNodes][3],
                                                      real64 (& gradN)[numFaces][3] )
{

  //real64 detJ = determinantJacobianTransformation( X );

  real64 J[3][3] = {{0}};

  J[0][0] = X[1][0] - X[0][0];
  J[0][1] = X[2][0] - X[0][0];
  J[0][2] = X[3][0] - X[0][0];

  J[1][0] = X[1][1] - X[0][1];
  J[1][1] = X[2][1] - X[0][1];
  J[1][2] = X[3][1] - X[0][1];

  J[2][0] = X[1][2] - X[0][2];
  J[2][1] = X[2][2] - X[0][2];
  J[2][2] = X[3][2] - X[0][2];

  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  real64 dNdXi[numFaces][3] = {{0}};

  real64 const r = quadratureParentCoords0( q );
  real64 const s = quadratureParentCoords1( q );
  real64 const t = quadratureParentCoords2( q );

  //std::cout << detJ1 << " " << detJ << " " << q << " r: " << r << " s: " << s << " t: " << t << std::endl; 

  dNdXi[0][0] = ( 1 - 2 * r - s - t ) * t; // dN0/dr
  dNdXi[0][1] = -r * t;                    // dN0/ds
  dNdXi[0][2] = ( 1 - r - s - 2 * t ) * r; // dN0/dt

  dNdXi[1][0] = ( 1 - 2 * r - s - t ) * s; // dN1/dr
  dNdXi[1][1] = ( 1 - r - 2 * s - t ) * r; // dN1/ds
  dNdXi[1][2] = -r * s;                    // dN1/dt

  dNdXi[2][0] = -t * s;                    // dN2/dr
  dNdXi[2][1] = ( 1 - r - 2 * s - t ) * t; // dN2/ds
  dNdXi[2][2] = ( 1 - r - s - 2 * t ) * s; // dN2/dt

  dNdXi[3][0] = t * s;                     // dN3/dr
  dNdXi[3][1] = r * t;                     // dN3/ds
  dNdXi[3][2] = r * s;                     // dN3/dt

  for( int fi=0; fi<numFaces; ++fi )
  {
    gradN[fi][0] = dNdXi[fi][0] * J[0][0] + dNdXi[fi][1] * J[1][0] + dNdXi[fi][2] * J[2][0];
    gradN[fi][1] = dNdXi[fi][0] * J[0][1] + dNdXi[fi][1] * J[1][1] + dNdXi[fi][2] * J[2][1];
    gradN[fi][2] = dNdXi[fi][0] * J[0][2] + dNdXi[fi][1] * J[1][2] + dNdXi[fi][2] * J[2][2];
  }

///////////////////////////////////////////////////////////////////////////////////////////////////
/*
  real64 dNdXi1[numNodes][3] = {{0}};

  dNdXi1[0][0] = -1; // dN0/dr
  dNdXi1[0][1] = -1; // dN0/ds
  dNdXi1[0][2] = -1; // dN0/dt

  dNdXi1[1][0] =  1; // dN1/dr
  dNdXi1[1][1] =  0; // dN1/ds
  dNdXi1[1][2] =  0; // dN1/dt

  dNdXi1[2][0] =  0; // dN2/dr
  dNdXi1[2][1] =  1; // dN2/ds
  dNdXi1[2][2] =  0; // dN2/dt

  dNdXi1[3][0] =  0; // dN3/dr
  dNdXi1[3][1] =  0; // dN3/ds
  dNdXi1[3][2] =  1; // dN3/dt

  real64 gradN1[numNodes][3] = {{0}};

  for( int fi=0; fi<numNodes; ++fi )
  {
    gradN1[fi][0] = dNdXi1[fi][0] * J[0][0] + dNdXi1[fi][1] * J[1][0] + dNdXi1[fi][2] * J[2][0];
    gradN1[fi][1] = dNdXi1[fi][0] * J[0][1] + dNdXi1[fi][1] * J[1][1] + dNdXi1[fi][2] * J[2][1];
    gradN1[fi][2] = dNdXi1[fi][0] * J[0][2] + dNdXi1[fi][1] * J[1][2] + dNdXi1[fi][2] * J[2][2];
  }

  real64 gradN2[numNodes][3] = {{0}};

  calcGradN( q, X, gradN2 );

  for( int fi=0; fi<numNodes; ++fi )
  {
    std::cout << gradN1[fi][0] << " " << gradN2[fi][0] << std::endl;
    std::cout << gradN1[fi][1] << " " << gradN2[fi][1] << std::endl;
    std::cout << gradN1[fi][2] << " " << gradN2[fi][2] << std::endl;
  }
  std::cout << std::endl;
*/
///////////////////////////////////////////////////////////////////////////////////////////////////

  return detJ * weight * quadratureWeight(q);
}

//*************************************************************************************************

GEOS_HOST_DEVICE
inline
real64
H1_Tetrahedron_Lagrange1_Gauss1::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{

  real64 detJ =  determinantJacobianTransformation( X );

  return detJ * weight * quadratureWeight(q);
}

/// @endcond

}
}

#endif //GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE1GAUSS1_HPP_
