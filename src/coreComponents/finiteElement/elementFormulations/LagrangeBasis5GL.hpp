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

#ifndef GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS5GL_HPP_
#define GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS5GL_HPP_

/**
 * @file LagrangeBasis5GL.hpp
 */

#include "common/DataTypes.hpp"

namespace geos
{
namespace finiteElement
{

/**
 * This class contains the implementation for a 5th order Lagrange
 * polynomial basis. The parent space is defined by:
 *
 *                 o-------------o------------------------o------------------------o------------------------o-------------------o  ---> xi
 *  Index:         0             1                        2                        3                        4                   5
 *  Coordinate:   -1    -sqrt(1/21(7+2*sqrt(7))) -sqrt(1/21(7-2*sqrt(7)))  sqrt(1/21(7-2*sqrt(7)))  sqrt(1/21(7-2*sqrt(7)))     1
 *
 */
class LagrangeBasis5GL
{
public:

  /// The number of support points for the basis
  constexpr static localIndex numSupportPoints = 6;

  /// sqrt(7)
  static constexpr real64 sqrt_7_ = 2.64575131106459059;

  /// sqrt( 7 + 2 * sqrt(7) )
  static constexpr real64 sqrt__7_plus_2sqrt7__ = 3.50592393273573196;

  /// sqrt( 7 - 2 * sqrt(7) )
  static constexpr real64 sqrt__7_mins_2sqrt7__ = 1.30709501485960033;

  /// sqrt( 7 + 2 * sqrt(7) )
  static constexpr real64 sqrt__7_plus_sqrt7_div2__ = 2.884939454396278;

  /// sqrt( 7 - 2 * sqrt(7) )
  static constexpr real64 sqrt__7_mins_sqrt7_div2__ = 2.382671682055189;

  /// sqrt(1/21)
  static constexpr real64 sqrt_inv21 = 0.218217890235992381;

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
      case 1:
      case 4:
        return (1.0/30.0)*(14.0-sqrt_7_);
      case 2:
      case 3:
        return (1.0/30.0)*(14.0+sqrt_7_);
      default:
        return 1.0/15.0;
    }
  }

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *   linear index of a support point. Here we hardcode GL points at order 5.
   * @param supportPointIndex The linear index of support point
   * @return parent coordinate in the xi0 direction.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 parentSupportCoord( const localIndex supportPointIndex )
  {

    real64 result=0.0;

    switch( supportPointIndex )
    {
      case 0:
        result = -1.0;
        break;

      case 1:
        result = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
        break;

      case 2:
        result = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
        break;

      case 3:
        result = sqrt_inv21*sqrt__7_mins_2sqrt7__;
        break;

      case 4:
        result = sqrt_inv21*sqrt__7_plus_2sqrt7__;
        break;

      case 5:
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
        return result = LagrangeBasis5GL::value0( xi );
        break;

      case 1:
        return result = LagrangeBasis5GL::value1( xi );
        break;

      case 2:
        return result = LagrangeBasis5GL::value2( xi );
        break;

      case 3:
        return result = LagrangeBasis5GL::value3( xi );
        break;

      case 4:
        return result = LagrangeBasis5GL::value4( xi );
        break;

      case 5:
        return result = LagrangeBasis5GL::value5( xi );
        break;

      default:
        break;
    }

    return result;

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
    /* Define the two GL points needed to compute the basis function at point index 0. Here we need the points
       at index 3,4 called lambda3, lambda4. */

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return (-21.0/16.0)*(xi*xi*xi*xi*xi-xi*xi*xi*xi-(lambda3*lambda3+lambda4*lambda4)*xi*xi*xi+(lambda3*lambda3+lambda4*lambda4)*xi*xi+
                         lambda3*lambda3*lambda4*lambda4*xi-lambda3*lambda3*lambda4*lambda4);
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
    /* Define the two GL points needed to compute the basis function at point index 1. Here we need the points
       at index 3 and 4 called lambda3, lambda4. */

    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );
    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );

    return ((21.0/16.0)*sqrt__7_mins_sqrt7_div2__)*(xi*xi*xi*xi*xi-lambda4*xi*xi*xi*xi-(lambda3*lambda3+1)*xi*xi*xi+lambda4*(lambda3*lambda3+1)*xi*xi+
                                                    lambda3*lambda3*xi-lambda4*lambda3*lambda3);
  }

/**
 * @brief The value of the basis function for the 2 support point.
 * @param xi The coordinate at which to evaluate the basis.
 * @return The value of the basis.
 */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value2( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 2. Here we need the points
       at index 3 and 4 called lambda3, lambda4. */

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return ((-21.0/16.0)*sqrt__7_plus_sqrt7_div2__)*(xi*xi*xi*xi*xi-lambda3*xi*xi*xi*xi-(lambda4*lambda4+1)*xi*xi*xi+lambda3*(lambda4*lambda4+1)*xi*xi+
                                                     lambda4*lambda4*xi-lambda3*lambda4*lambda4);
  }

/**
 * @brief The value of the basis function for the 3 support point.
 * @param xi The coordinate at which to evaluate the basis.
 * @return The value of the basis.
 */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value3( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 3. Here we need the points
       at index 3 and 4 called lambda1, lambda2. */

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return ((21.0/16.0)*sqrt__7_plus_sqrt7_div2__)*(xi*xi*xi*xi*xi+lambda3*xi*xi*xi*xi-(lambda4*lambda4+1)*xi*xi*xi-lambda3*(lambda4*lambda4+1)*xi*xi+
                                                    lambda4*lambda4*xi+lambda3*lambda4*lambda4);
  }

/**
 * @brief The value of the basis function for the 4 support point.
 * @param xi The coordinate at which to evaluate the basis.
 * @return The value of the basis.
 */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value4( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 4. Here we need the points
       at index 3 and 4 called lambda3, lambda4. */

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return ((-21.0/16.0)*sqrt__7_mins_sqrt7_div2__)*(xi*xi*xi*xi*xi+lambda4*xi*xi*xi*xi-(lambda3*lambda3+1)*xi*xi*xi-lambda4*(lambda3*lambda3+1)*xi*xi+
                                                     lambda3*lambda3*xi+lambda4*lambda3*lambda3);
  }

/**
 * @brief The value of the basis function for the 5 support point.
 * @param xi The coordinate at which to evaluate the basis.
 * @return The value of the basis.
 */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 value5( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 5. Here we need the points
       at index 3 and 4 called lambda3, lambda4. */

    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );
    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );

    return (21.0/16.0)*(xi*xi*xi*xi*xi+xi*xi*xi*xi-(lambda4*lambda4+lambda3*lambda3)*xi*xi*xi-(lambda3*lambda3+lambda4*lambda4)*xi*xi+
                        lambda3*lambda3*lambda4*lambda4*xi+lambda4*lambda4*lambda3*lambda3);
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

    real64 result=0.0;

    switch( index )
    {
      case 0:
        result = LagrangeBasis5GL::gradient0( xi );
        break;

      case 1:
        result = LagrangeBasis5GL::gradient1( xi );
        break;

      case 2:
        result = LagrangeBasis5GL::gradient2( xi );
        break;

      case 3:
        result = LagrangeBasis5GL::gradient3( xi );
        break;

      case 4:
        result = LagrangeBasis5GL::gradient4( xi );
        break;

      case 5:
        result = LagrangeBasis5GL::gradient5( xi );
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
   * @return The gradient of basis function at point 0.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient0( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return (-21.0/16.0)*(5.0*xi*xi*xi*xi-4.0*xi*xi*xi-3.0*(lambda3*lambda3+lambda4*lambda4)*xi*xi+2.0*(lambda3*lambda3+lambda4*lambda4)*xi+lambda3*lambda3*lambda4*lambda4);

  }

/**
 * @brief The gradient of the basis function for support point 1 evaluated at
 *   a point along the axes.
 * @param xi The coordinate at which to evaluate the gradient.
 * @return The gradient of basis function at point 1.
 */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient1( const real64 xi )
  {

    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );
    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );

    return (21.0/16.0)*sqrt__7_mins_sqrt7_div2__*(5.0*xi*xi*xi*xi-4.0*lambda4*xi*xi*xi-3.0*(lambda3*lambda3+1.0)*xi*xi+2.0*lambda4*(lambda3*lambda3+1.0)*xi+lambda3*lambda3);

  }

  /**
   * @brief The gradient of the basis function for support point 2 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 2.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient2( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return (-21.0/16.0)*sqrt__7_plus_sqrt7_div2__*(5.0*xi*xi*xi*xi-4.0*lambda3*xi*xi*xi-3.0*(lambda4*lambda4+1.0)*xi*xi+2.0*lambda3*(lambda4*lambda4+1.0)*xi+lambda4*lambda4);

  }

  /**
   * @brief The gradient of the basis function for support point 3 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 3.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient3( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return (21.0/16.0)*sqrt__7_plus_sqrt7_div2__*(5.0*xi*xi*xi*xi+4.0*lambda3*xi*xi*xi-3.0*(lambda4*lambda4+1.0)*xi*xi-2*lambda3*(lambda4*lambda4+1.0)*xi+lambda4*lambda4);

  }

  /**
   * @brief The gradient of the basis function for support point 4 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 0.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient4( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return (-21.0/16.0)*sqrt__7_mins_sqrt7_div2__*(5.0*xi*xi*xi*xi+4.0*lambda4*xi*xi*xi-3.0*(lambda3*lambda3+1.0)*xi*xi-2.0*lambda4*(lambda3*lambda3+1.0)*xi+lambda3*lambda3);
  }

  /**
   * @brief The gradient of the basis function for support point 5 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 5.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr static real64 gradient5( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord( 4 );
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord( 3 );

    return (21.0/16.0)*(5.0*xi*xi*xi*xi+4.0*xi*xi*xi-3.0*(lambda3*lambda3+lambda4*lambda4)*xi*xi-2.0*(lambda3*lambda3+lambda4*lambda4)*xi+lambda3*lambda3*lambda4*lambda4);

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
          case 0: return -7.5000000000000000000;
          case 1: return -1.7863649483390948939;
          case 2: return 0.48495104785356916930;
        }
        break;
      case 1:
        switch( p )
        {
          case 0: return 10.14141593631966928023;
          case 1: return 0.0;
          case 2: return -1.72125695283023338321;
        }
        break;
      case 2:
        switch( p )
        {
          case 0: return -4.03618727030534800527;
          case 1: return 2.5234267774294554319088;
          case 2: return 0.0;
        }
        break;
      case 3:
        switch( p )
        {
          case 0: return 2.2446846481761668242712;
          case 1: return -1.1528281585359293413318;
          case 2: return 1.7529619663678659788775;
        }
        break;
      case 4:
        switch( p )
        {
          case 0: return -1.3499133141904880992312;
          case 1: return 0.6535475074298001672007;
          case 2: return -0.7863566722232407374395;
        }
        break;
      case 5:
        switch( p )
        {
          case 0: return 0.500000000000000000000;
          case 1: return -0.2377811779842313638052;
          case 2: return 0.2697006108320389724720;
        }
        break;
    }
    return 0;
  }

  /* UNCRUSTIFY-OFF */

  /**
   * @struct TensorProduct2D
   *
   * A 2-dimensional basis formed from the tensor product of the 1d basis.
   *
   *  30        31       32       33       34         35                              ____________________________________________________
   *    o--------o--------o--------o--------o--------o                               |Node      xi0                        xi1            |
   *    |                                            |                               |=====     ===                        ===            |
   *    |                                            |                               |  0       -1                         -1             |
   *    |                                            |                               |  1   -sqrt(1/21(7+/sqrt(7))         -1             |
   *    |                                            |                               |  2    -sqrt(1/21(7-/sqrt(7))        -1             |
   * 24 o     25 o     26 o        o 27     o 28     o 29                            |  3    sqrt(1/21(7-/sqrt(7))         -1             |
   *    |                                            |                               |  4    sqrt(1/21(7+/sqrt(7))         -1             |
   *    |                                            |                               |  5        1                         -1             |
   *    |                                            |                               |  6       -1                 -sqrt(1/21(7+/sqrt(7)) |
   * 18 o     19 o     20 o        o 21     o 22     o 23                            |  7   -sqrt(1/21(7+/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) |
   *    |                                            |                               |  8   -sqrt(1/21(7-/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) |
   *    |                                            |                               |  9    sqrt(1/21(7-/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) |
   *    |                                            |                               | 10    sqrt(1/21(7+/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) |
   * 12 o     13 o     14 o        o 15     o 16     o 17                            | 11        1                 -sqrt(1/21(7+/sqrt(7)) |
   *    |                                            |                               | ..       ..                         ..             |
   *    |                                            |                               | ..       ..                         ..             |
   *    |                                            |                               | 30       -1                          1             |
   * 6  o     7  o     8  o        o 9      o 10     o 11                            | 31   -sqrt(1/21(7+/sqrt(7))          1             |
   *    |                                            |                               | 32    -sqrt(1/21(7-/sqrt(7))         1             |
   *    |                                            |                  xi1          | 33    sqrt(1/21(7-/sqrt(7))          1             |
   *    |                                            |                   |           | 34    sqrt(1/21(7+/sqrt(7))          1             |
   *    |                                            |                   |           | 35        1                          1             |
   *    o--------o--------o--------o--------o--------o                   |           |____________________________________________________|
   *   0         1        2        3        4         5                  o----- xi0
   */

  /* UNCRUSTIFY-ON */

  struct TensorProduct2D
  {
    /// The number of support points in the basis.
    constexpr static localIndex numSupportPoints = 36;

    /**
     * @brief Calculates the linear index for support/quadrature points from ij
     *   coordinates.
     * @param i The index in the xi0 direction (0,5)
     * @param j The index in the xi1 direction (0,5)
     * @return The linear index of the support/quadrature point (0-35)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j )
    {
      return i + 6 * j;
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

      i1 = linearIndex/6;

      i0 = linearIndex%6;

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
      for( int a=0; a<6; ++a )
      {
        for( int b=0; b<6; ++b )
        {
          const int lindex = LagrangeBasis5GL::TensorProduct2D::linearIndex( a, b );
          N[ lindex ] = LagrangeBasis5GL::value( a, coords[0] ) *
                        LagrangeBasis5GL::value( b, coords[1] );
        }
      }
    }

  };

  /**
   * @struct TensorProduct3D
   *
   * A 3-dimensional basis formed from the tensor product of the 1d basis.
   *
   *                210      211      212      213      214      215    _______________________________________________________
   *                   o--------o--------o--------o--------o--------o  |Node      xi0                        xi1            xi2|
   *                  /.                                           /|  |=====     ===                        ===            ===|
   *            204  / .  205      206      207      208      209 / |  |  0       -1                         -1             -1 |
   *                o  .     o        o        o        o        o  |  |  1   -sqrt(1/21(7+/sqrt(7))         -1             -1 |
   *               /   o                                        /   |  |  2    -sqrt(1/21(7-/sqrt(7))        -1             -1 |
   *         198  /    .199     200      201      202      203 /    o  |  3    sqrt(1/21(7-/sqrt(7))         -1             -1 |
   *             o     .  o        o        o        o        o     |  |  4    sqrt(1/21(7+/sqrt(7))         -1             -1 |
   *            /      .                                     /      |  |  5        1                         -1             -1 |
   *      192  /   193 o     194      195      196      197 /    o  |  |  6       -1                 -sqrt(1/21(7+/sqrt(7)) -1 |
   *          o        o        o        o        o        o        o  |  7   -sqrt(1/21(7+/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) -1 |
   *         /         .                                  /         |  |  8   -sqrt(1/21(7-/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) -1 |
   *    186 /    187   .  188      189      190      191 /    o     |  |  9    sqrt(1/21(7-/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) -1 |
   *       o        o  o     o        o        o        o        o  |  | 10    sqrt(1/21(7+/sqrt(7)) -sqrt(1/21(7+/sqrt(7)) -1 |
   *      /            .                               /            o  | 11        1                 -sqrt(1/21(7+/sqrt(7)) -1 |
   * 180 /    181      . 182    183      184      185 /    o        |  | ..       ..                         ..             .. |
   *    o--------o--------o--------o--------o--------o        o     |  | ..       ..                         ..             .. |
   *    |           o  .                             |           o  |  | 204      -1                  sqrt(1/21(7+/sqrt(7)) 1  |
   *    |  o           o        o        o        o  |  o  o        o  | 205  -sqrt(1/21(7+/sqrt(7))  sqrt(1/21(7+/sqrt(7)) 1  |
   *    |     o        .                             |     o        |  | 206  -sqrt(1/21(7-/sqrt(7))  sqrt(1/21(7-/sqrt(7)) 1  |
   *    |        o     .                             |        o     |  | 207  sqrt(1/21(7+/sqrt(7))   sqrt(1/21(7-/sqrt(7)) 1  |
   *    o           o  .                             o           o  |  | 208  sqrt(1/21(7-/sqrt(7))   sqrt(1/21(7+/sqrt(7)) 1  |
   *    |  o           .                             |  o           |  | 209       1                  sqrt(1/21(7+/sqrt( *  1  |
   *    |     o        o--------o--------o--------o--|-----o--------o  | 210      -1                          1             1  |
   *    |        o    ,30       31      32        33 |     34 o    /35 | 211  -sqrt(1/21(7+/sqrt(7))          1             1  |
   *    o            ,                               o            /    | 212  -sqrt(1/21(7-/sqrt(7))          1             1  |
   *    |  o        o        o         o       o     |  o        o     | 213   sqrt(1/21(7-/sqrt(7))          1             1  |
   *    |     o    ,24       25        26      27    |  28 o    /29    | 214   sqrt(1/21(7+/sqrt(7))          1             1  |
   *    |         ,                                  |         /       | 215       1                          1             1  |
   *    o        o        o         o       o     22 o        o *      |_______________________________________________________|
   *    |  o    ,18       19        20      21       |  o    /23
   *    |      ,                                     |      /
   *    |     o        o         o       o        o  |     o
   *    o    ,12       13        14      15       16 o    /17
   *    |   ,                                        |   /
   *    |  o        o        o        o        o     |  o               xi2
   *    | ,6        7        8        9        10    | /11               |
   *    |,                                           |/                  | / xi1
   *    o--------o--------o--------o--------o--------o                   |/
   *    0        1        2        3        4        5                   o----- xi0
   */

  /* UNCRUSTIFY-ON */

  struct TensorProduct3D
  {
    /// The number of support points in the basis.
    constexpr static localIndex numSupportPoints = 216;

    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,5)
     * @param j The index in the xi1 direction (0,5)
     * @param k The index in the xi2 direction (0,5)
     * @return The linear index of the support/quadrature point (0-215)
     */
    GEOS_HOST_DEVICE
    GEOS_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j,
                                      const int k )
    {
      return i + 6 * j + 36 * k;
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

      i2 = linearIndex/36;

      i1 = (linearIndex%36)/6;

      i0 = (linearIndex%36)%6;

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
      for( int a=0; a<6; ++a )
      {
        for( int b=0; b<6; ++b )
        {
          for( int c=0; c<6; ++c )
          {
            const int lindex = LagrangeBasis5GL::TensorProduct3D::linearIndex( a, b, c );
            N[ lindex ] = LagrangeBasis5GL::value( a, coords[0] ) *
                          LagrangeBasis5GL::value( b, coords[1] ) *
                          LagrangeBasis5GL::value( c, coords[2] );
          }
        }
      }
    }

  };

};

}
}


#endif /* GEOS_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS5GL_HPP_ */
