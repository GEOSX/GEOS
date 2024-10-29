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

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_

/**
 * @file LagrangeBasis1.hpp
 */

#include "common/DataTypes.hpp"

namespace geos
{
namespace finiteElement
{

/**
 * This class contains the implementation for a first order (linear) Lagrange
 * polynomial basis. The parent space is defined by:
 *
 *                 o-------------o  ---> xi
 *  Index:         0             1
 *  Coordinate:   -1             1
 *
 */
class LagrangeBasis1
{
public:
  /// The number of support points for the basis
  constexpr static localIndex numSupportPoints = 2;

  /**
   * @brief The value of the weight for the given support point
   * @param q The index of the support point
   * @return The value of the weight
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 weight( const int q )
  {
    GEOS_UNUSED_VAR( q );
    return 1.0;
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
    return -1.0 + 2.0 * (supportPointIndex & 1);
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
    return 0.5 + 0.5 * xi * parentSupportCoord( index );
  }


  /**
   * @brief The value of the basis function for the 0 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value0( const real64 xi )
  {
    return 0.5 - 0.5 * xi;
  }

  /**
   * @brief The value of the basis function for the 1 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value1( const real64 xi )
  {
    return 0.5 + 0.5 * xi;
  }

  /**
   * @brief The value of the bubble basis function.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 valueBubble( const real64 xi )
  {
    return 1.0 - pow( xi, 2 );
  }


  /**
   * @brief The gradient of the basis function for a support point evaluated at
   *   a point along the axes.
   * @param index The index of the support point associated with the basis
   *   function.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient( const int index,
                                    const real64 xi )
  {
    GEOS_UNUSED_VAR( xi );
    return 0.5 * parentSupportCoord( index );
  }

  /**
   * @brief The gradient of the basis function for support point 0 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function (-0.5)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient0( const real64 xi )
  {
    GEOS_UNUSED_VAR( xi );
    return -0.5;
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function (0.5)
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient1( const real64 xi )
  {
    GEOS_UNUSED_VAR( xi );
    return 0.5;
  }

  /**
   * @brief The gradient of the bubble basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradientBubble( const real64 xi )
  {
    GEOS_UNUSED_VAR( xi );
    return -0.5*xi;
  }

  /**
   * @brief The gradient of the basis function for a support point evaluated at
   *   a given support point. By symmetry, p is assumed to be in 0, ..., (N-1)/2.
   *   in the case of the first-order basis, this value is independent of p.
   * @param q The index of the basis function
   * @return The gradient of basis function.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradientAt( const int q,
                                      const int )
  {
    return q == 0 ? -0.5 : 0.5;
  }

  /**
   * @struct TensorProduct2D
   *
   * A 2-dimensional basis formed from the tensor product of the 1d basis.
   *
   *               2                   3
   *                o-----------------o                           _______________
   *                |                 |                          |Node   xi0  xi1|
   *                |                 |                          |=====  ===  ===|
   *                |                 |                          | 0     -1   -1 |
   *                |                 |                          | 1      1   -1 |
   *                |                 |            xi1           | 2     -1    1 |
   *                |                 |            |             | 3      1    1 |
   *                |                 |            |             |_______________|
   *                o-----------------o            |
   *               0                   1           ------ xi0
   *
   */
  struct TensorProduct2D
  {
    /// The number of support points in the basis.
    constexpr static localIndex numSupportPoints = 4;

    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @return The linear index of the support/quadrature point (0-3)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j )
    {
      return i + 2 * j;
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
      i0 = ( linearIndex & 1 );
      i1 = ( linearIndex & 2 ) >> 1;
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
      for( int a=0; a<2; ++a )
      {
        for( int b=0; b<2; ++b )
        {
          const int lindex = LagrangeBasis1::TensorProduct2D::linearIndex( a, b );
          N[ lindex ] = LagrangeBasis1::value( a, coords[0] ) *
                        LagrangeBasis1::value( b, coords[1] );
        }
      }
    }

    /**
     * @brief The value of the bubble basis function evaluated at a
     *   point along the axes.
     *
     * @param coords The coordinates (in the parent frame) at which to evaluate the basis
     * @param N Array to hold the value of the basis functions.
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    static void valueBubble( real64 const (&coords)[2],
                             real64 (& N)[1] )
    {
      N[0] = LagrangeBasis1::valueBubble( coords[0] ) *
             LagrangeBasis1::valueBubble( coords[1] );
    }

    /**
     * @brief The parent coordinates for a support point in the xi0 direction.
     * @param linearIndex The linear index of the support point
     * @return
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static real64 parentCoords0( localIndex const linearIndex )
    {
      return -1.0 + 2.0 * (linearIndex & 1);
    }

    /**
     * @brief The parent coordinates for a support point in the xi1 direction.
     * @param linearIndex The linear index of the support point
     * @return
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static real64 parentCoords1( localIndex const linearIndex )
    {
      return -1.0 + ( linearIndex & 2 );
    }

  };

  /**
   * @struct TensorProduct3D
   *
   * A 3-dimensional basis formed from the tensor product of the 1d basis.
   *
   *                  6                   7                       ____________________
   *                   o-----------------o                       |Node   xi0  xi1  xi2|
   *                  /.                /|                       |=====  ===  ===  ===|
   *                 / .               / |                       | 0     -1   -1   -1 |
   *              4 o-----------------o 5|                       | 1      1   -1   -1 |
   *                |  .              |  |                       | 2     -1    1   -1 |
   *                |  .              |  |                       | 3      1    1   -1 |
   *                |  .              |  |                       | 4     -1   -1    1 |
   *                |  .              |  |                       | 5      1   -1    1 |
   *                |2 o..............|..o 3       xi2           | 6     -1    1    1 |
   *                | ,               | /          |             | 7      1    1    1 |
   *                |,                |/           | / xi1       |____________________|
   *                o-----------------o            |/
   *               0                   1           ------ xi0
   *
   */
  struct TensorProduct3D
  {
    /// The number of support points in the basis.
    constexpr static localIndex numSupportPoints = 8;

    /// The number of support faces in the basis.
    constexpr static localIndex numSupportFaces = 6;

    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @param k The index in the xi2 direction (0,1)
     * @return The linear index of the support/quadrature point (0-7)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j,
                                      const int k )
    {
      return i + 2 * j + 4 * k;
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
      i0 = ( linearIndex & 1 );
      i1 = ( linearIndex & 2 ) >> 1;
      i2 = ( linearIndex & 4 ) >> 2;
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
    static void value( real64 const (&coords)[3],
                       real64 (& N)[numSupportPoints] )
    {
      for( int a=0; a<2; ++a )
      {
        for( int b=0; b<2; ++b )
        {
          for( int c=0; c<2; ++c )
          {
            const int lindex = LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );
            N[ lindex ] = LagrangeBasis1::value( a, coords[0] ) *
                          LagrangeBasis1::value( b, coords[1] ) *
                          LagrangeBasis1::value( c, coords[2] );
          }
        }
      }
    }

    /**
     * @brief The value of the bubble basis function for a support face evaluated at a
     *   point along the axes.
     *
     * @param coords The coordinates (in the parent frame) at which to evaluate the basis
     * @param N Array to hold the value of the basis functions at each support face.
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    static void valueFaceBubble( real64 const (&coords)[3],
                                 real64 (& N)[numSupportFaces] )
    {
      for( int a=0; a<2; ++a )
      {
        N[ a*5 ]   = LagrangeBasis1::valueBubble( coords[0] ) *
                     LagrangeBasis1::valueBubble( coords[1] ) *
                     LagrangeBasis1::value( a, coords[2] );

        N[ a*3+1 ] = LagrangeBasis1::valueBubble( coords[0] ) *
                     LagrangeBasis1::value( a, coords[1] ) *
                     LagrangeBasis1::valueBubble( coords[2] );

        N[ a+2 ]   = LagrangeBasis1::value( a, coords[0] ) *
                     LagrangeBasis1::valueBubble( coords[1] ) *
                     LagrangeBasis1::valueBubble( coords[2] );
      }
    }

    /**
     * @brief The value of the bubble basis function derivatives for a support face evaluated at a
     *   point along the axes.
     *
     * @param coords The coordinates (in the parent frame) at which to evaluate the basis
     * @param dNdXi Array to hold the value of the basis function derivatives at each support face.
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    static void gradientFaceBubble( real64 const (&coords)[3],
                                    real64 (& dNdXi)[numSupportFaces][3] )
    {
      for( int a=0; a<2; ++a )
      {
        dNdXi[ a*5 ][0]   = LagrangeBasis1::gradientBubble( coords[0] ) *
                            LagrangeBasis1::valueBubble( coords[1] ) *
                            LagrangeBasis1::value( a, coords[2] );
        dNdXi[ a*5 ][1]   = LagrangeBasis1::valueBubble( coords[0] ) *
                            LagrangeBasis1::gradientBubble( coords[1] ) *
                            LagrangeBasis1::value( a, coords[2] );
        dNdXi[ a*5 ][2]   = LagrangeBasis1::valueBubble( coords[0] ) *
                            LagrangeBasis1::valueBubble( coords[1] ) *
                            LagrangeBasis1::gradient( a, coords[2] );

        dNdXi[ a*3+1 ][0] = LagrangeBasis1::gradientBubble( coords[0] ) *
                            LagrangeBasis1::value( a, coords[1] ) *
                            LagrangeBasis1::valueBubble( coords[2] );
        dNdXi[ a*3+1 ][1] = LagrangeBasis1::valueBubble( coords[0] ) *
                            LagrangeBasis1::gradient( a, coords[1] ) *
                            LagrangeBasis1::valueBubble( coords[2] );
        dNdXi[ a*3+1 ][2] = LagrangeBasis1::valueBubble( coords[0] ) *
                            LagrangeBasis1::value( a, coords[1] ) *
                            LagrangeBasis1::gradientBubble( coords[2] );

        dNdXi[ a+2 ][0]   = LagrangeBasis1::gradient( a, coords[0] ) *
                            LagrangeBasis1::valueBubble( coords[1] ) *
                            LagrangeBasis1::valueBubble( coords[2] );
        dNdXi[ a+2 ][1]   = LagrangeBasis1::value( a, coords[0] ) *
                            LagrangeBasis1::gradientBubble( coords[1] ) *
                            LagrangeBasis1::valueBubble( coords[2] );
        dNdXi[ a+2 ][2]   = LagrangeBasis1::value( a, coords[0] ) *
                            LagrangeBasis1::valueBubble( coords[1] ) *
                            LagrangeBasis1::gradientBubble( coords[2] );
      }
    }

    /**
     * @brief The parent coordinates for a support point in the xi0 direction.
     * @param linearIndex The linear index of the support point
     * @return
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static real64 parentCoords0( localIndex const linearIndex )
    {
      return -1.0 + 2.0 * (linearIndex & 1);
    }

    /**
     * @brief The parent coordinates for a support point in the xi1 direction.
     * @param linearIndex The linear index of the support point
     * @return
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static real64 parentCoords1( localIndex const linearIndex )
    {
      return -1.0 + ( linearIndex & 2 );
    }

    /**
     * @brief The parent coordinates for a support point in the xi2 direction.
     * @param linearIndex The linear index of the support point
     * @return
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static real64 parentCoords2( localIndex const linearIndex )
    {
      return -1.0 + 0.5 * ( linearIndex & 4 );
    }

  };

};

}
}


#endif /* GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_ */
