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

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS3GL_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS3GL_HPP_

/**
 * @file LagrangeBasis3GL.hpp
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
 *                 o---------o--------o---------o  ---> xi
 *  Index:         0         1        2         3
 *  Coordinate:   -1    -1/sqrt(5) 1/sqrt(5)    1
 *
 */
class LagrangeBasis3GL
{
public:
  /// The number of support points for the basis
  constexpr static localIndex numSupportPoints = 4;

  /// sqrt(5)
  constexpr static real64 sqrt5 = 2.2360679774997897;

  /**
   * @brief The value of the weight for the given support point
   * @param q The index of the support point
   * @return The value of the weight
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 weight( const int q )
  {
    switch( q )
    {
      case 1:
      case 2:
        return 5.0/6.0;
      default:
        return 1.0/6.0;
    }
  }

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *   linear index of a support point.
   * @param supportPointIndex The linear index of support point
   * @return parent coordinate in the xi0 direction.
   */
  GEOS_HOST_DEVICE
  inline
  // MODIF1 : Harcoding the Gauss-Lobatto coordinates and return the right one
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
        result = -1.0/sqrt5;
        break;
      case 2:
        result = 1.0/sqrt5;
        break;
      case 3:
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
  inline
//MODIF3 : Change the  way to return the base function evaluated at the desired coord
  constexpr static real64 value( const int index,
                                 const real64 xi )
  {
    real64 result=0.0;

    switch( index )
    {
      case 0:
        result = LagrangeBasis3GL::value0( xi );
        break;
      case 1:
        result = LagrangeBasis3GL::value1( xi );
        break;
      case 2:
        result = LagrangeBasis3GL::value2( xi );
        break;
      case 3:
        result = LagrangeBasis3GL::value3( xi );
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
  inline
  //MODFI4 : Implemented new base functions and their derivative for Q3
  constexpr static real64 value0( const real64 xi )
  {
    return -(5.0/8.0)*(xi*xi*xi-xi*xi-(1.0/5.0)*xi+1.0/5.0);
  }

  /**
   * @brief The value of the basis function for support point 1.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 value1( const real64 xi )
  {
    return (5.0*sqrt5/8.0)*(xi*xi*xi-(1.0/sqrt5)*xi*xi-xi+1.0/sqrt5);
  }

  /**
   * @brief The value of the basis function for support point 2.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 value2( const real64 xi )
  {
    return -(5.0*sqrt5/8.0)*(xi*xi*xi+(1.0/sqrt5)*xi*xi-xi-1.0/sqrt5);
  }

  /**
   * @brief The value of the basis function for support point 3.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 value3( const real64 xi )
  {
    return (5.0/8.0)*(xi*xi*xi+xi*xi-(1.0/5.0)*xi-1.0/5.0);
  }


  /**
   * @brief The gradient of the basis function for a support point evaluated at a
   *   point along the axes.
   * @param index The index of the support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of basis function.
   */
  GEOS_HOST_DEVICE
  inline
  //MODIF5 : New function returning the derivated base function at desired coord
  constexpr static real64 gradient( const int index,
                                    const real64 xi )
  {
    real64 result=0.0;

    switch( index )
    {
      case 0:
        result = LagrangeBasis3GL::gradient0( xi );
        break;
      case 1:
        result = LagrangeBasis3GL::gradient1( xi );
        break;
      case 2:
        result = LagrangeBasis3GL::gradient2( xi );
        break;
      case 3:
        result = LagrangeBasis3GL::gradient3( xi );
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
  inline
  constexpr static real64 gradient0( const real64 xi )
  {
    return -(5.0/8.0)*(3.0*xi*xi-2.0*xi-(1.0/5.0));
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 gradient1( const real64 xi )
  {
    return (5.0*sqrt5/8.0)*(3.0*xi*xi-(2.0/sqrt5)*xi-1.0);
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 gradient2( const real64 xi )
  {
    return -(5.0*sqrt5/8.0)*(3.0*xi*xi+(2.0/sqrt5)*xi-1.0);
  }

  /**
   * @brief The gradient of the basis function for support point 3 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOS_HOST_DEVICE
  inline
  constexpr static real64 gradient3( const real64 xi )
  {
    return (5.0/8.0)*(3.0*xi*xi+2.0*xi-(1.0/5.0));;
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
        return p == 0 ? -3.0 : -0.80901699437494742410;
      case 1:
        return p == 0 ? 4.0450849718747371205 : 0.0;
      case 2:
        return p == 0 ? -1.5450849718747371205 : 1.1180339887498948482;
      case 3:
        return p == 0 ? 0.5 : -0.30901699437494742410;
      default:
        return 0;
    }
  }

  /**
   * @class TensorProduct2D
   *
   *                                                                  _____________________________
   *        12         13        14         15                       |Node      xi0         xi1    |
   *          o---------o---------o---------o                        |=====     ===         ===    |
   *          |                             |                        |  0       -1          -1     |
   *          |                             |                        |  1   -1/sqrt(5)      -1     |
   *          |                             |                        |  2    1/sqrt(5)      -1     |
   *          |                             |                        |  3        1          -1     |
   *        8 o       9 o         o 10      o 11                     |  4       -1      -1/sqrt(5) |
   *          |                             |                        |  5   -1/sqrt(5)  -1/sqrt(5) |
   *          |                             |                        |  6    1/sqrt(5)  -1/sqrt(5) |
   *          |                             |                        |  7        1      -1/sqrt(5) |
   *        4 o       5 o         o 6       o 7                      |  8       -1       1/sqrt(5) |
   *          |                             |                        |  9   -1/sqrt(5)   1/sqrt(5) |
   *          |                             |            xi1         | 10    1/sqrt(5)   1/sqrt(5) |
   *          |                             |            |           | 11        1       1/sqrt(5) |
   *          |                             |            |           | 12       -1           1     |
   *          o---------o---------o---------o            |           | 13   -1/sqrt(5)       1     |
   *         0          1         2          3           o----- xi0  | 14    1/sqrt(5)       1     |
   *                                                                 | 15        1           1     |
   *                                                                 |_____________________________|
   *
   */
  struct TensorProduct2D
  {
    /// The number of support points in the 2D tensor product
    constexpr static localIndex numSupportPoints = 16;

    /**
     * @brief Calculates the linear index for support/quadrature points from ij
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @return The linear index of the support/quadrature point (0-15)
     */
    GEOS_HOST_DEVICE
    inline
    constexpr static int linearIndex( const int i,
                                      const int j )
    {
      return i + 4 * j;
    }



    /**
     * @brief Calculate the Cartesian/TensorProduct index given the linear index
     *   of a support point.
     * @param linearIndex The linear index of support point
     * @param i0 The Cartesian index of the support point in the xi0 direction.
     * @param i1 The Cartesian index of the support point in the xi1 direction.
     */
    GEOS_HOST_DEVICE
    inline
    constexpr static void multiIndex( int const linearIndex,
                                      int & i0,
                                      int & i1 )
    {

      i1 = linearIndex/4;

      i0 = linearIndex%4;

    }


    /**
     * @brief The value of the basis function for a support point evaluated at a
     *   point along the axes.
     *
     * @param coords The coordinates (in the parent frame) at which to evaluate the basis
     * @param N Array to hold the value of the basis functions at each support point.
     */
    GEOS_HOST_DEVICE
    inline
    static void value( const real64 (& coords)[2],
                       real64 (& N)[numSupportPoints] )
    {
      for( int a=0; a<4; ++a )
      {
        for( int b=0; b<4; ++b )
        {
          const int lindex = LagrangeBasis3GL::TensorProduct2D::linearIndex( a, b );
          N[ lindex ] = LagrangeBasis3GL::value( a, coords[0] ) *
                        LagrangeBasis3GL::value( b, coords[1] );
        }
      }
    }
  };

  /**
   * @class TensorProduct3D
   *                                                                  _____________________________________
   *                                                                 |Node      xi0         xi1         xi2|
   *                                                                 |=====     ===         ===         ===|
   *                                                                 |  0       -1          -1          -1 |
   *                                                                 |  1   -1/sqrt(5)      -1          -1 |
   *                                                                 |  2    1/sqrt(5)      -1          -1 |
   *              60       61         62        63                   |  3        1          -1          -1 |
   *                o---------o---------o---------o                  |  4       -1      -1/sqrt(5)      -1 |
   *            56 /.     57        58        59 /|                  |  5   -1/sqrt(5)  -1/sqrt(5)      -1 |
   *              o .       o         o         o |                  |  6    1/sqrt(5)  -1/sqrt(5)      -1 |
   *          52 /  .   53        54        55 /  |                  |  7        1      -1/sqrt(5)      -1 |
   *            o   .     o         o         o   |                  |  8       -1       1/sqrt(5)      -1 |
   *        48 /    o 49      o 50      o 51 /    o                  |  9   -1/sqrt(5)   1/sqrt(5)      -1 |
   *          o---------o---------o---------o     |                  | 10    1/sqrt(5)   1/sqrt(5)      -1 |
   *          |   o .       o         o     |   o |                  | 11        1       1/sqrt(5)      -1 |
   *          |     .                       |     |                  | 12       -1           1          -1 |
   *          | o   o     o   o     o   o   | o   o                  | 13   -1/sqrt(5)       1          -1 |
   *          |     .                       |     |                  | 14    1/sqrt(5)       1          -1 |
   *          o   o .   o   o     o   o     o   o |                  | 15        1           1          -1 |
   *          |     .                       |     |                  | ..       ..          ..          .. |
   *          | o   .     o         o       | o   |                  | ..       ..          ..          .. |
   *          |     o.........o.........o...|.....o                  | 55        1      -1/sqrt(5)       1 |
   *          o    ,12  o     13  o     14  o    /15                 | 56       -1       1/sqrt(5)       1 |
   *          |   o         o         o     |   o                    | 57   -1/sqrt(5)   1/sqrt(5)       1 |
   *          |  ,8         9         10    |  /11       xi2         | 58    1/sqrt(5)   1/sqrt(5)       1 |
   *          | o         o         o       | o          |           | 59        1       1/sqrt(5)       1 |
   *          |,4         5         6       |/7          | / xi1     | 60       -1           1           1 |
   *          o---------o---------o---------o            |/          | 61   -1/sqrt(5)       1           1 |
   *         0         1         2         3             o----- xi0  | 62    1/sqrt(5)       1           1 |
   *                                                                 | 63        1           1           1 |
   *                                                                 |_____________________________________|
   *
   */
  struct TensorProduct3D
  {
    /// The number of support points in the 3D tensor product
    constexpr static localIndex numSupportPoints = 64;

    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @param k The index in the xi2 direction (0,1)
     * @return The linear index of the support/quadrature point (0-63)
     */
    GEOS_HOST_DEVICE
    inline
    //MODIF6 : Change linearIndex for 64 nodes
    constexpr static int linearIndex( const int i,
                                      const int j,
                                      const int k )
    {
      return i + 4 * j + 16 * k;
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
    inline
    //MODIF7 : Change calcul of multiIndex
    constexpr static void multiIndex( int const linearIndex,
                                      int & i0,
                                      int & i1,
                                      int & i2 )
    {

      i2 = linearIndex/16;

      i1 = (linearIndex%16)/4;

      i0 = (linearIndex%16)%4;

    }


    /**
     * @brief The value of the basis function for a support point evaluated at a
     *   point along the axes.
     *
     * @param coords The coordinates (in the parent frame) at which to evaluate the basis
     * @param N Array to hold the value of the basis functions at each support point.
     */
    GEOS_HOST_DEVICE
    inline
    static void value( const real64 (& coords)[3],
                       real64 (& N)[numSupportPoints] )
    {
      for( int a=0; a<4; ++a )
      {
        for( int b=0; b<4; ++b )
        {
          for( int c=0; c<4; ++c )
          {
            const int lindex = LagrangeBasis3GL::TensorProduct3D::linearIndex( a, b, c );
            N[ lindex ] = LagrangeBasis3GL::value( a, coords[0] ) *
                          LagrangeBasis3GL::value( b, coords[1] ) *
                          LagrangeBasis3GL::value( c, coords[2] );
          }
        }
      }
    }
  };
};


}
}


#endif /* GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS3GL_HPP_ */
