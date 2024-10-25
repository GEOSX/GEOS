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

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS2_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS2_HPP_
/**
 * @file LagrangeBasis2.hpp
 */

#include "common/DataTypes.hpp"

namespace geos
{
namespace finiteElement
{

/**
 * This class contains the implementation for a second order (quadratic) Lagrange
 * polynomial basis. The parent space is defined by:
 *
 *                 o-------------o-------------o  ---> xi
 *  Index:         0             1             2
 *  Coordinate:   -1             0             1
 *
 */
class LagrangeBasis2
{
public:
  /// The number of support points for the basis
  constexpr static localIndex numSupportPoints = 3;

  /**
   * @brief The value of the weight for the given support point
   * @param q The index of the support point
   * @return The value of the weight
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 weight( const int q )
  {
    switch( q )
    {
      case 0:
      case 2:
        return 1.0/3.0;
      default:
        return 4.0/3.0;
    }
  }

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *   linear index of a support point.
   * @param supportPointIndex The linear index of support point
   * @return parent coordinate in the xi0 direction.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 parentSupportCoord( const localIndex supportPointIndex )
  {
    switch( supportPointIndex )
    {
      case 0:
        return -1.0;
        break;
      case 2:
        return 1.0;
      case 1:
      default:
        return 0.0;
    }
  }

  /**
   * @brief The value of the basis function for a support point evaluated at a
   *   point along the axes.
   * @param index The index of the support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of basis function.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value( const int index,
                                 const real64 xi )
  {

    switch( index )
    {
      case 0:
        return value0( xi );
      case 2:
        return value2( xi );
      case 1:
      default:
        return value1( xi );
    }
  }

  /**
   * @brief The value of the basis function for support point 0.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value0( const real64 xi )
  {
    const real64 xi_div2 = 0.5 * xi;
    return -xi_div2 + xi_div2 * xi;
  }

  /**
   * @brief The value of the basis function for support point 1.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value1( const real64 xi )
  {
    return 1.0 - xi * xi;
  }

  /**
   * @brief The value of the basis function for support point 2.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value2( const real64 xi )
  {
    const real64 xi_div2 = 0.5 * xi;
    return xi_div2 + xi_div2 * xi;
  }

  /**
   * @brief The gradient of the basis function for support point 0 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient0( const real64 xi )
  {
    return -0.5 + xi;
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient1( const real64 xi )
  {
    return -2 * xi;
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient2( const real64 xi )
  {
    return 0.5 + xi;
  }

  /**
   * @brief The gradient of the basis function for a support point evaluated at a
   *   point along the axes.
   * @param index The index of the support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of basis function.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient( const int index,
                                    const real64 xi )
  {
    switch( index )
    {
      case 0:
        return gradient0( xi );
      case 2:
        return gradient2( xi );
      case 1:
      default:
        return gradient1( xi );
    }
  }

  /**
   * @brief The gradient of the basis function for a support point evaluated at
   *   a given support point. By symmetry, p is assumed to be in 0, ..., (N-1)/2
   * @param q The index of the basis function
   * @param p The index of the support point
   * @return The gradient of basis function.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradientAt( const int q,
                                      const int p )
  {
    switch( q )
    {
      case 0:
        return p == 0 ? -1.5 : -0.5;
      case 1:
        return p == 0 ? 2.0 : 0.0;
      case 2:
        return p == 0 ? -0.5 : 0.5;
      default:
        return 0;
    }
  }

  /**
   * @class TensorProduct2D
   *
   *         6               7               8
   *          o--------------o--------------o                          ________________
   *          |                             |                         |Node   xi0  xi1 |
   *          |                             |                         |=====  ===  === |
   *          |                             |                         |  0    -1   -1  |
   *          |                             |                         |  1     0   -1  |
   *          |                             |                         |  2     1   -1  |
   *          |                             |                         |  3    -1    0  |
   *        3 o              o 4            o 5                       |  4     0    0  |
   *          |                             |                         |  5     1    0  |
   *          |                             |                         |  6    -1    1  |
   *          |                             |                         |  7     0    1  |
   *          |                             |            xi1          |  8     1    1  |
   *          |                             |            |            |________________|
   *          |                             |            |
   *          o--------------o--------------o            |
   *         0               1               2           o----- xi0
   *
   *
   *
   */
  struct TensorProduct2D
  {

    /// The number of support points in the basis.
    constexpr static localIndex numSupportPoints = 9;

    /**
     * @brief Calculates the linear index for support/quadrature points from ij
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @return The linear index of the support/quadrature point (0-8)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j )
    {
      return i + 3 * j;
    }



    /**
     * @brief Calculate the Cartesian/TensorProduct index given the linear index
     *   of a support point.
     * @param linearIndex The linear index of support point
     * @param i0 The Cartesian index of the support point in the xi0 direction.
     * @param i1 The Cartesian index of the support point in the xi1 direction.
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static void multiIndex( const int linearIndex,
                                      int & i0,
                                      int & i1 )
    {

      i1 = ( ( linearIndex * 22 ) >> 6 );
      //i1 = a/3;

      i0 = linearIndex - i1 * 3;
    }

    /**
     * @brief The value of the basis function for a support point evaluated at a
     *   point along the axes.
     *
     * @param coords The coordinates (in the parent frame) at which to evaluate the basis
     * @param N Array to hold the value of the basis functions at each support point.
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    static void value( real64 const (&coords)[2],
                       real64 (& N)[numSupportPoints] )
    {
      for( int a=0; a<3; ++a )
      {
        for( int b=0; b<3; ++b )
        {
          const int lindex = LagrangeBasis2::TensorProduct2D::linearIndex( a, b );
          N[ lindex ] = LagrangeBasis2::value( a, coords[0] ) *
                        LagrangeBasis2::value( b, coords[1] );
        }
      }
    }
  };

  /**
   * @class TensorProduct3D
   *
   *                                                                  ____________________
   *                                                                 |Node   xi0  xi1  xi2|
   *                                                                 |=====  ===  ===  ===|
   *                                                                 |  0    -1   -1   -1 |
   *                                                                 |  1     0   -1   -1 |
   *                                                                 |  2     1   -1   -1 |
   *              24              25               26                |  3    -1    0   -1 |
   *                o--------------o--------------o                  |  4     0    0   -1 |
   *               /.                            /|                  |  5     1    0   -1 |
   *              / .                           / |                  |  6    -1    1   -1 |
   *          21 o  .           o 22        23 o  |                  |  7     0    1   -1 |
   *            /   .                         /   |                  |  8     1    1   -1 |
   *           /    .         19             /    |                  |  9    -1   -1    0 |
   *       18 o--------------o--------------o 20  |                  | 10     0   -1    0 |
   *          |     o              o        |     o                  | 11     1   -1    0 |
   *          |     .15             16      |     |17                | 12    -1    0    0 |
   *          |     .                       |     |                  | 13     0    0    0 |
   *          |  o  .           o           |  o  |                  | 14     1    0    0 |
   *          |   12.            13         |   14|                  | 15    -1    1    0 |
   *          |     .                       |     |                  | 16     0    1    0 |
   *        9 o     .        o 10           o 11  |                  | 17     1    1    0 |
   *          |     o..............o........|.....o                  | 18    -1   -1    1 |
   *          |    , 6              7       |    / 8                 | 19     0   -1    1 |
   *          |   ,                         |   /                    | 20     1   -1    1 |
   *          |  o              o           |  o         xi2         | 21    -1    0    1 |
   *          | , 3              4          | / 5        |           | 22     0    0    1 |
   *          |,                            |/           | / xi1     | 23     1    0    1 |
   *          o--------------o--------------o            |/          | 24    -1    1    1 |
   *         0                1              2           o----- xi0  | 25     0    1    1 |
   *                                                                 | 26     1    1    1 |
   *                                                                 |____________________|
   *
   */
  struct TensorProduct3D
  {

    /// The number of support points in the basis.
    constexpr static localIndex numSupportPoints = 27;

    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @param k The index in the xi2 direction (0,1)
     * @return The linear index of the support/quadrature point (0-26)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j,
                                      const int k )
    {
      return i + 3 * j + 9 * k;
    }



    /**
     * @brief Calculate the Cartesian/TensorProduct index given the linear index
     *   of a support point.
     * @param linearIndex The linear index of support point
     * @param i0 The Cartesian index of the support point in the xi0 direction.
     * @param i1 The Cartesian index of the support point in the xi1 direction.
     * @param i2 The Cartesian index of the support point in the xi2 direction.
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static void multiIndex( const int linearIndex,
                                      int & i0,
                                      int & i1,
                                      int & i2 )
    {
      i2 = ( linearIndex * 29 ) >> 8;
      //i2 = a/9;

      i1 = ( ( linearIndex * 22 ) >> 6 ) - i2 * 3;
      //i1 = a/3 - i2 * 3;

      i0 = linearIndex - i1 * 3 - i2 * 9;
    }

    /**
     * @brief The value of the basis function for a support point evaluated at a
     *   point along the axes.
     *
     * @param coords The coordinates (in the parent frame) at which to evaluate the basis
     * @param N Array to hold the value of the basis functions at each support point.
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    static void value( const real64 (& coords)[3],
                       real64 (& N)[numSupportPoints] )
    {
      for( int a=0; a<3; ++a )
      {
        for( int b=0; b<3; ++b )
        {
          for( int c=0; c<3; ++c )
          {
            const int lindex = LagrangeBasis2::TensorProduct3D::linearIndex( a, b, c );
            N[ lindex ] = LagrangeBasis2::value( a, coords[0] ) *
                          LagrangeBasis2::value( b, coords[1] ) *
                          LagrangeBasis2::value( c, coords[2] );
          }
        }
      }
    }
  };
};

}
}

#endif /* GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS2_HPP_ */
