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

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS4GL_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS4GL_HPP_

/**
 * @file LagrangeBasis4GL.hpp
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
 *                 o------o------o------o------o  ---> xi
 *  Index:         0      1      2      3      4
 *  Coordinate:   -1 -sqrt(3/7)  0  sqrt(3/7)  1
 *
 */
class LagrangeBasis4GL
{
public:
  /// The number of support points for the basis
  constexpr static localIndex numSupportPoints = 5;

  /// sqrt(3/7)
  constexpr static real64 sqrt3_7 = 0.6546536707079771;

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
      case 4:
        return 1.0/10.0;
      case 1:
      case 3:
        return 49.0/90.0;
      default:
        return 32.0/45.0;
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
  // depending on the supportPointIndex value
  //Switch case
  constexpr static real64 parentSupportCoord( const localIndex supportPointIndex )
  {
    real64 result=0.0;

    switch( supportPointIndex )
    {
      case 0:
        result = -1.0;
        break;
      case 1:
        result = -sqrt3_7;
        break;
      case 2:
        result = 0.0;
        break;
      case 3:
        result = sqrt3_7;
        break;
      case 4:
        result = 1.0;
        break;
      default:
        break;
    }

    return result;
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
    real64 result=0.0;

    switch( index )
    {
      case 0:
        result = LagrangeBasis4GL::value0( xi );
        break;
      case 1:
        result = LagrangeBasis4GL::value1( xi );
        break;
      case 2:
        result = LagrangeBasis4GL::value2( xi );
        break;
      case 3:
        result = LagrangeBasis4GL::value3( xi );
        break;
      case 4:
        result = LagrangeBasis4GL::value4( xi );
        break;
      default:
        break;
    }

    return result;
  }


  /**
   * @brief The value of the basis function for support point 0.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  //MODFI4 : Implemented new base functions and their derivative for Q3
  constexpr static real64 value0( const real64 xi )
  {
    return (1.0/8.0)*(-1.0+xi)*xi*(-3.0+7.0*xi*xi);
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
    return (49.0/24.0)*(sqrt3_7-xi)*xi*(-1.0+xi*xi);
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
    return (1.0/3.0)*(3.0-10.0*xi*xi+7.0*xi*xi*xi*xi);
  }

  /**
   * @brief The value of the basis function for support point 3.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value3( const real64 xi )
  {
    return -(49.0/24.0)*(sqrt3_7+xi)*xi*(-1.0+xi*xi);
  }

  /**
   * @brief The value of the basis function for support point 4.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value4( const real64 xi )
  {
    return (1.0/8.0)*(1.0+xi)*xi*(-3.0+7.0*xi*xi);
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
    real64 result=0.0;

    switch( index )
    {
      case 0:
        result = LagrangeBasis4GL::gradient0( xi );
        break;
      case 1:
        result = LagrangeBasis4GL::gradient1( xi );
        break;
      case 2:
        result = LagrangeBasis4GL::gradient2( xi );
        break;
      case 3:
        result = LagrangeBasis4GL::gradient3( xi );
        break;
      case 4:
        result = LagrangeBasis4GL::gradient4( xi );
        break;
      default:
        break;
    }

    return result;
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
    return (1.0/8.0)*(3.0+xi*(-6.0+7.0*xi*(-3.0+4.0*xi)));
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
    return (49.0/24.0)*(-sqrt3_7+xi*(2.0+3.0*sqrt3_7*xi-4.0*xi*xi));
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
    return (4.0/3.0)*xi*(-5.0+7.0*xi*xi);
  }

  /**
   * @brief The gradient of the basis function for support point 3 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient3( const real64 xi )
  {
    return (49.0/24.0)*(sqrt3_7+xi*(2.0-3.0*sqrt3_7*xi-4.0*xi*xi));
  }

  /**
   * @brief The gradient of the basis function for support point 4 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient4( const real64 xi )
  {
    return (1.0/8.0)*(-3.0+xi*(-6.0+7.0*xi*(3.0+4.0*xi)));
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
        switch( p )
        {
          case 0: return -5.0000000000000000000;
          case 1: return -1.2409902530309828578;
          case 2: return 0.37500000000000000000;
        }
        break;
      case 1:
        switch( p )
        {
          case 0: return 6.7565024887242400038;
          case 1: return 0.0;
          case 2: return -1.3365845776954533353;
        }
        break;
      case 2:
        switch( p )
        {
          case 0: return -2.6666666666666666667;
          case 1: return 1.7457431218879390501;
          case 2: return 0.0;
        }
        break;
      case 3:
        switch( p )
        {
          case 0: return 1.4101641779424266628;
          case 1: return -0.7637626158259733344;
          case 2: return 1.3365845776954533353;
        }
        break;
      case 4:
        switch( p )
        {
          case 0: return -0.50000000000000000000;
          case 1: return 0.25900974696901714215;
          case 2: return -0.37500000000000000000;
        }
        break;
    }
    return 0;
  }

  /**
   * @class TensorProduct2D
   *
   *                                                                    _____________________________
   *        20       21      22         23     24                      |Node      xi0         xi1    |
   *          o-------o-------o-------o-------o                        |=====     ===         ===    |
   *          |                               |                        |  0       -1          -1     |
   *          |                               |                        |  1   -sqrt(3/7)      -1     |
   *       15 o    16 o    17 o       o 18    o 19                     |  2        0          -1     |
   *          |                               |                        |  3    sqrt(3/7)      -1     |
   *          |                               |                        |  4        1          -1     |
   *          |                               |                        |  5       -1      -sqrt(3/7) |
   *       10 o    11 o    12 o       o 13    o 14                     |  6   -sqrt(3/7)  -sqrt(3/7) |
   *          |                               |                        |  7        0      -sqrt(3/7) |
   *          |                               |                        |  8    sqrt(3/7)  -sqrt(3/7) |
   *          |                               |                        |  9        1      -sqrt(3/7) |
   *        5 o     6 o     7 o       o 8     o 9         xi1          | 10       -1           0     |
   *          |                               |            |           | 11   -sqrt(3/7)       0     |
   *          |                               |            |           | ........................... |
   *          o-------o-------o-------o-------o            |           | 13        0           1     |
   *         0        1       2       3       4            o----- xi0  | 14    sqrt(3/7)       1     |
   *                                                                   | 15        1           1     |
   *                                                                   |_____________________________|
   *
   */
  struct TensorProduct2D
  {
    /// The number of support points in the 2D tensor product
    constexpr static localIndex numSupportPoints = 25;

    /**
     * @brief Calculates the linear index for support/quadrature points from ij
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @return The linear index of the support/quadrature point (0-15)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j )
    {
      return i + 5 * j;
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
    constexpr static void multiIndex( int const linearIndex,
                                      int & i0,
                                      int & i1 )
    {

      i1 = linearIndex/5;

      i0 = linearIndex%5;

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
    static void value( const real64 (& coords)[2],
                       real64 (& N)[numSupportPoints] )
    {
      for( int a=0; a<5; ++a )
      {
        for( int b=0; b<5; ++b )
        {
          const int lindex = LagrangeBasis4GL::TensorProduct2D::linearIndex( a, b );
          N[ lindex ] = LagrangeBasis4GL::value( a, coords[0] ) *
                        LagrangeBasis4GL::value( b, coords[1] );
        }
      }
    }
  };

  /**
   * @class TensorProduct3D
   *
   *
   *
   *                                                                  _____________________________________
   *                120      121     122     123       124           |Node      xi0         xi1         xi2 |
   *                  o-------o-------o-------o-------o              |=====     ===         ===         === |
   *                 /.                              /|              |   0       -1          -1          -1 |
   *            115 o .  116o    117o    118o    119o |              |   1   -sqrt(3/7)      -1          -1 |
   *               /  o                            /  o              |   2        0          -1          -1 |
   *          110 o   .111o    112o    113o    114o   |              |   3    sqrt(3/7)      -1          -1 |
   *             /  o .                          /  o |              |   4        1          -1          -1 |
   *        105 o     . o106    o107    o108 109o     |              |   5       -1      -sqrt(3/7)      -1 |
   *           /  o   o      102     103    104/  o   o              |   6   -sqrt(3/7)  -sqrt(3/7)      -1 |
   *      100 o-------o-------o-------o-------o       |              |   7        0      -sqrt(3/7)      -1 |
   *          | o   o . 101                   | o   o |              |   8    sqrt(3/7)  -sqrt(3/7)      -1 |
   *          |       .                       |       |              |   9        1      -sqrt(3/7)      -1 |
   *          o   o   o       o       o       o   o   o              |  10       -1           0          -1 |
   *          |       .                       |       |              |  11   -sqrt(3/7)       0          -1 |
   *          | o   o .20     21      22    23| o   o |24            |  12        0           0          -1 |
   *          |       o.......o.......o.......|.......o              |  13    sqrt(3/7)       0          -1 |
   *          o   o  ,o       o       o       o   o  /               |  14        1           0          -1 |
   *          |     o       o       o       o |     o                |  ..       ..          ..          .. |
   *          | o  ,15      13      17      18| o  /19               |  ..       ..          ..          .. |
   *          |   o       o       o       o   |   o                  |  ..       ..          ..          .. |
   *          o  ,10  o   11  o   12  o   13  o  /14     xi2         | 121        -1          1           1 |
   *          | o       o       o       o     | o        |           | 122    -sqrt(3/7)      1           1 |
   *          |,5       6       7       8     |/9        | / xi1     | 123         0          1           1 |
   *          o-------o-------o-------o-------o          |/          | 124     sqrt(3/7)      1           1 |
   *         0        1       2       3        4         o----- xi0  | 125         1          1           1 |
   *                                                                 |______________________________________|
   *
   *
   */
  struct TensorProduct3D
  {
    /// The number of support points in the 3D tensor product
    constexpr static localIndex numSupportPoints = 125;

    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @param k The index in the xi2 direction (0,1)
     * @return The linear index of the support/quadrature point (0-124)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j,
                                      const int k )
    {
      return i + 5 * j + 25 * k;
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
    constexpr static void multiIndex( int const linearIndex,
                                      int & i0,
                                      int & i1,
                                      int & i2 )
    {

      i2 = linearIndex/25;

      i1 = (linearIndex%25)/5;

      i0 = (linearIndex%25)%5;

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
      for( int a=0; a<5; ++a )
      {
        for( int b=0; b<5; ++b )
        {
          for( int c=0; c<5; ++c )
          {
            const int lindex = LagrangeBasis4GL::TensorProduct3D::linearIndex( a, b, c );
            N[ lindex ] = LagrangeBasis4GL::value( a, coords[0] ) *
                          LagrangeBasis4GL::value( b, coords[1] ) *
                          LagrangeBasis4GL::value( c, coords[2] );
          }
        }
      }
    }
  };
};


}
}


#endif /* GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS4GL_HPP_ */
