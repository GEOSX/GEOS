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
 * @file H1_TriangleFace_Lagrange1_Gauss1.hpp
 */

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_H1TRIANGLEFACELAGRANGE1GAUSS1_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_H1TRIANGLEFACELAGRANGE1GAUSS1_HPP_

#include "FiniteElementBase.hpp"
#include "LagrangeBasis1.hpp"


namespace geos
{
namespace finiteElement
{

#define TRIANGLE_QUADRATURE_POINTS 4

/**
 * This class contains the kernel accessible functions specific to the
 * H1-conforming nodal linear triangular face finite element with a
 * 1-point Gaussian quadrature rule.
 *
 *
 *          3                                =====  ==  ==
 *           +                               Node   r   s
 *           |\_                             =====  ==  ==
 *           |  \_                           0      0   0
 *           |    \_           s             1      1   0
 *           |      \_         |             2      0   1
 *           |        \_       |             =====  ==  ==
 *           |          \      |
 *           +-----------+     *------r
 *          0              1
 *
 */
class H1_TriangleFace_Lagrange1_Gauss1 final : public FiniteElementBase
{
public:

  /// The type of basis used for this element
  using BASIS = LagrangeBasis1;

  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 3;
  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
#if TRIANGLE_QUADRATURE_POINTS == 6
  constexpr static localIndex numQuadraturePoints = 6;
#elif TRIANGLE_QUADRATURE_POINTS == 4
  constexpr static localIndex numQuadraturePoints = 4;
#elif TRIANGLE_QUADRATURE_POINTS == 1
  constexpr static localIndex numQuadraturePoints = 1;
#endif

  virtual ~H1_TriangleFace_Lagrange1_Gauss1() override
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
  static localIndex getNumQuadraturePoints( StackVariables const & stack )
  {
    GEOS_UNUSED_VAR( stack );
    return numQuadraturePoints;
  }

  GEOS_HOST_DEVICE
  virtual localIndex getNumSupportPoints() const override
  {
    return numNodes;
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

  GEOS_HOST_DEVICE
  virtual localIndex getMaxSupportPoints() const override
  {
    return maxSupportPoints;
  }

  /**
   * @brief Calculate shape functions values at a single point.
   * @param[in] coords The parent coordinates at which to evaluate the shape function value
   * @param[out] N The shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcN( real64 const (&coords)[2],
                     real64 ( &N )[numNodes] );


  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *          point.
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
   * @brief Calculate shape bubble functions values at a given point in the parent space.
   * @param pointCoord coordinates of the given point.
   * @param N An array to pass back the shape function values.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void calcBubbleN( real64 const (&pointCoord)[2],
                           real64 (& N)[1] )
  {

    real64 const r  = pointCoord[0];
    real64 const s  = pointCoord[1];

    N[0] = (1.0 - r - s) * r * s;

  }

  /**
   * @brief Calculate shape bubble functions values at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values.
   */
  GEOS_HOST_DEVICE
  inline
  static void calcBubbleN( localIndex const q,
                           real64 (& N)[1] )
  {

    real64 const qCoords[2] = {quadratureParentCoords0(q), quadratureParentCoords1(q) };

    calcBubbleN( qCoords, N );
  }

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
   * @brief Empty method, here for compatibility with methods that require a stabilization of the
   * grad-grad bilinear form.
   * @tparam NUMDOFSPERTRIALSUPPORTPOINT Number of degrees of freedom for each support point.
   * @tparam UPPER If true only the upper triangular part of @p matrix is modified.
   * @param stack Stack variables as filled by @ref setupStack.
   * @param matrix The matrix that needs to be stabilized.
   * @param scaleFactor Optional scaling of the stabilization matrix.
   */
  template< localIndex NUMDOFSPERTRIALSUPPORTPOINT, bool UPPER >
  GEOS_HOST_DEVICE
  inline
  static void addGradGradStabilization( StackVariables const & stack,
                                        real64 ( &matrix )
                                        [maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT]
                                        [maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT],
                                        real64 const & scaleFactor );

  /**
   * @brief Calculate the node permutation between the parent element and the geometric element.
   *   Note: The optimal location for this calculation is yet to be determined.
   * @param permutation An array to return the node permutation.
   */
  GEOS_HOST_DEVICE
  inline
  static void getPermutation( int (& permutation)[numNodes] )
  {
    permutation[0] = 0;
    permutation[1] = 1;
    permutation[2] = 2;
  }

private:
  /// The area of the element in the parent configuration.
  constexpr static real64 parentArea = 0.5;

  /// The weight of each quadrature point.
  //constexpr static real64 weight = parentArea / numQuadraturePoints;
  constexpr static real64 weight = parentArea;

  /**
   * @brief Calculate the weights
   * @param q
   * @return weight.
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 quadratureWeight( localIndex const q )
  {

#if TRIANGLE_QUADRATURE_POINTS == 6
    real64 const w[numQuadraturePoints] = { 1.0/6.0,
                                            1.0/6.0,
                                            1.0/6.0,
                                            1.0/6.0,
                                            1.0/6.0,
                                            1.0/6.0 };

#elif TRIANGLE_QUADRATURE_POINTS == 4

    real64 const w[numQuadraturePoints] = {-0.562500000000000,
                                            0.520833333333333,
                                            0.520833333333333,
                                            0.520833333333333 };

#elif TRIANGLE_QUADRATURE_POINTS == 1

    real64 const w[numQuadraturePoints] = { 1.0 };
#endif

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

#if TRIANGLE_QUADRATURE_POINTS == 6
    real64 const qCoords[numQuadraturePoints] = { 0.659027622374092,
                                                  0.109039009072877,
                                                  0.231933368553031,
                                                  0.659027622374092,
                                                  0.109039009072877,
                                                  0.231933368553031 };

#elif TRIANGLE_QUADRATURE_POINTS == 4
    real64 const qCoords[numQuadraturePoints] = { 0.333333333333333,
                                                  0.600000000000000,
                                                  0.200000000000000,
                                                  0.200000000000000 };

#elif TRIANGLE_QUADRATURE_POINTS == 1
    real64 const qCoords[numQuadraturePoints] = { 1/3 };

#endif


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

#if TRIANGLE_QUADRATURE_POINTS == 6
    real64 const qCoords[numQuadraturePoints] = { 0.231933368553031,
                                                  0.659027622374092,
                                                  0.109039009072877,
                                                  0.109039009072877,
                                                  0.231933368553031,
                                                  0.659027622374092 };

#elif TRIANGLE_QUADRATURE_POINTS == 4
    real64 const qCoords[numQuadraturePoints] = { 0.333333333333333,
                                                  0.200000000000000,
                                                  0.600000000000000,
                                                  0.200000000000000 };

#elif TRIANGLE_QUADRATURE_POINTS == 1
    real64 const qCoords[numQuadraturePoints] = { 1/3 };

#endif

    return qCoords[q];
  }

};

/// @cond Doxygen_Suppress

template< localIndex NUMDOFSPERTRIALSUPPORTPOINT, bool UPPER >
GEOS_HOST_DEVICE
inline
void H1_TriangleFace_Lagrange1_Gauss1::
  addGradGradStabilization( StackVariables const & stack,
                            real64 ( & matrix )
                            [maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT]
                            [maxSupportPoints * NUMDOFSPERTRIALSUPPORTPOINT],
                            real64 const & scaleFactor )
{
  GEOS_UNUSED_VAR( stack );
  GEOS_UNUSED_VAR( matrix );
  GEOS_UNUSED_VAR( scaleFactor );
}

GEOS_HOST_DEVICE
inline
void
H1_TriangleFace_Lagrange1_Gauss1::
  calcN( real64 const (&coords)[2],
         real64 ( & N )[numNodes] )
{
  real64 const r  = coords[0];
  real64 const s  = coords[1];

  N[0] = 1.0 - r - s;
  N[1] = r;
  N[2] = s;
}

GEOS_HOST_DEVICE
inline
void
H1_TriangleFace_Lagrange1_Gauss1::
  calcN( localIndex const q,
         real64 (& N)[numNodes] )
{

  real64 const qCoords[2] = {quadratureParentCoords0(q), quadratureParentCoords1(q) };

  calcN(qCoords, N);

}

GEOS_HOST_DEVICE
inline
void H1_TriangleFace_Lagrange1_Gauss1::
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
H1_TriangleFace_Lagrange1_Gauss1::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  //GEOS_UNUSED_VAR( q );
  real64 n[3] = { ( X[1][1] - X[0][1] ) * ( X[2][2] - X[0][2] ) - ( X[2][1] - X[0][1] ) * ( X[1][2] - X[0][2] ),
                  ( X[2][0] - X[0][0] ) * ( X[1][2] - X[0][2] ) - ( X[1][0] - X[0][0] ) * ( X[2][2] - X[0][2] ),
                  ( X[1][0] - X[0][0] ) * ( X[2][1] - X[0][1] ) - ( X[2][0] - X[0][0] ) * ( X[1][1] - X[0][1] )};

  return sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] ) * weight * quadratureWeight(q);
}

/// @endcond

#undef TRIANGLE_QUADRATURE_POINTS

}
}
#endif //GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_H1TRIANGLEFACELAGRANGE1GAUSS1_HPP_
