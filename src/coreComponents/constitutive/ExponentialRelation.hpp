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
 * @file ExponentialRelations.hpp
 */

#ifndef GEOS_CONSITUTIVE_EXPONENTIALRELATION_HPP_
#define GEOS_CONSITUTIVE_EXPONENTIALRELATION_HPP_

#include "common/DataTypes.hpp"
#include "common/format/EnumStrings.hpp"

#include <cmath>

namespace geos
{

namespace constitutive
{

/**
 * @enum Enumeration describing available approximations of the exponent
 */
enum class ExponentApproximationType : integer
{
  Full,
  Linear,
  Quadratic
};

ENUM_STRINGS( ExponentApproximationType,
              "exponential",
              "linear",
              "quadratic" );

namespace detail
{

// bring into scope for ADL
using std::exp;
using std::log;
using std::sqrt;

/**
 * @brief Specializations of this struct provide compute methods for different exponent approximations
 * @tparam T the data type
 * @tparam EAT the type/order of exponent approximation (linear, quadratic, full)
 * @tparam INVERSE whether to compute the inverse of the exponent
 */
template< typename T, ExponentApproximationType EAT, bool INVERSE = false >
struct ExponentialCompute
{};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Full, false >
{
  // - Two variables
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    y = y0 * exp( alpha * (x - x0));
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
  {
    y = y0 * exp( alpha * (x - x0));
    dy_dx = alpha * y;
  }

  // - Three variables
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & w0, const T & y0, const T & alpha, const T & beta, const T & x, const T & w, T & y )
  {
    y = y0 * exp( alpha * (x - x0) + beta * (w - w0) );
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & w0, const T & y0, const T & alpha, const T & beta, const T & x, const T & w, T & y, T & dy_dx, T & dy_dw )
  {
    y = y0 * exp( alpha * (x - x0) + beta * (w - w0) );
    dy_dx = alpha * y;
    dy_dw = beta * y;
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Full, true >
{
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + log( y / y0 ) / alpha;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * log( y / y0 );
    dx_dy = alpha_inv / y;
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Quadratic, false >
{
  // - Two variables
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    const T z = T( 1.0 ) + alpha * (x - x0);
    y = y0 / T( 2.0 ) * (T( 1.0 ) + z * z);
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
  {
    const T z = T( 1.0 ) + alpha * (x - x0);
    y = y0 / T( 2.0 ) * (T( 1.0 ) + z * z);
    dy_dx = alpha * y0 * z;
  }

  // - Three variables
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & w0, const T & y0, const T & alpha, const T & beta, const T & x, const T & w, T & y )
  {
    const T z1 = T( 1.0 ) + alpha * (x - x0);
    const T z2 = T( 1.0 ) + beta * (w - w0);
    y = y0 / T( 4.0 ) * (T( 1.0 ) + z1 * z1) * (T( 1.0 ) + z2 * z2);
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & w0, const T & y0, const T & alpha, const T & beta, const T & x, const T & w, T & y, T & dy_dx, T & dy_dw )
  {
    const T z1 = T( 1.0 ) + alpha * (x - x0);
    const T z2 = T( 1.0 ) + beta * (w - w0);
    y = y0 / T( 4.0 ) * (T( 1.0 ) + z1 * z1) * (T( 1.0 ) + z2 * z2);

    dy_dx = y0 / T( 2.0 ) * (T( 1.0 ) + z2 * z2) * alpha * z1;
    dy_dw = y0 / T( 2.0 ) * (T( 1.0 ) + z1 * z1) * beta * z2;
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Quadratic, true >
{
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    const T z = sqrt( T( 2.0 ) * y / y0 - T( 1.0 ));
    x = x0 + (z - 1) / alpha;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    const T z = sqrt( T( 2.0 ) * y / y0 - T( 1.0 ));
    x = x0 + alpha_inv * (z - 1);
    dx_dy = alpha_inv / (y0 * z);
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Linear, false >
{
  // - Two variables
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    y = y0 * (T( 1.0 ) + alpha * (x - x0));
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y, T & dy_dx )
  {
    y = y0 * (T( 1.0 ) + alpha * (x - x0));
    dy_dx = alpha * y0;
  }

  // - Three variables
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & w0, const T & y0, const T & alpha, const T & beta, const T & x, const T & w, T & y )
  {
    y = y0 * (T( 1.0 ) + alpha * (x - x0)) * (T( 1.0 ) + beta * (w - w0));
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & w0, const T & y0, const T & alpha, const T & beta, const T & x, const T & w, T & y, T & dy_dx, T & dy_dw )
  {
    y = y0 * (T( 1.0 ) + alpha * (x - x0)) * (T( 1.0 ) + beta * (w - w0));

    dy_dx = y0 * alpha * (T( 1.0 ) + beta * (w - w0));
    dy_dw = y0 * beta * (T( 1.0 ) + alpha * (x - x0));
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Linear, true >
{
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + (y / y0 - 1) / alpha;
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x, T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * (y / y0 - 1);
    dx_dy = alpha_inv / y0;
  }
};

}

/**
 * @class ExponentialRelation
 *
 * Describes a simple exponential relationship with three or four parameters.
 * The form of the relationship with three parameters: \f$ y = y0 * exp(a * (x - x0)) \f$.
 * The form of the relationship with four parameters: \f$ y = y0 * exp(a * (x - x0) + b * (w - w0)) \f$.
 * Optional linear and quadratic Taylor series approximation of the exponent.
 * Implements point and batch direct (y(x) and y(x, w)) and inverse (x(y) only for three parameters) compute functions and provides
 * derivatives.
 *
 * @tparam T scalar real-valued type used in computation
 * @tparam EAT the type/order of exponent approximation (linear, quadratic, full, linearThreeVar, quadraticThreeVar, fullThreeVar)
 * @tparam VAR integer no. of variables (only allows 2 and 3 now), default value is 2
 */
template< typename T, ExponentApproximationType EAT, integer VAR = 2 >
class ExponentialRelation
{
public:

  static constexpr ExponentApproximationType type = EAT;

  // *** constructors ***

  /**
   * @brief Default constructor. Sets \f$ x0 = 0, y0 = 1, alpha = 1 \f$
   */
  GEOS_HOST_DEVICE
  ExponentialRelation()
    : ExponentialRelation( T( 0 ), T( 1 ), T( 1 ) )
  {}

  /**
   * @brief Constructor with arguments
   * @param x0 variable shift
   * @param y0 scaling coefficient
   * @param alpha exponential coefficient
   */
  GEOS_HOST_DEVICE
  ExponentialRelation( T x0, T y0, T alpha )
  {
    GEOS_ERROR_IF( VAR != 2, GEOS_FMT( "The constructor is inconsistent with the number of variables {}", VAR ) );

    setCoefficients( x0, y0, alpha );
  }

  /**
   * @brief Constructor with arguments
   * @param x0 variable shift
   * @param w0 variable shift
   * @param y0 scaling coefficient
   * @param alpha exponential coefficient
   * @param beta  exponential coefficient
   */
  GEOS_HOST_DEVICE
  ExponentialRelation( T x0, T w0, T y0, T alpha, T beta )
  {
    GEOS_ERROR_IF( VAR != 3, GEOS_FMT( "The constructor is inconsistent with the number of variables {}", VAR ) );

    setCoefficients( x0, w0, y0, alpha, beta );
  }

  // *** setters ***

  /**
   * @brief Sets the relation coefficients
   * @param x0 variable shift
   * @param y0 scaling coefficient
   * @param alpha exponential coefficient
   */
  GEOS_HOST_DEVICE
  void setCoefficients( T x0, T y0, T alpha )
  {
    m_x0 = x0;
    m_y0 = y0;
    m_alpha = alpha;
  }

  /**
   * @brief Constructor with arguments
   * @param x0 variable shift
   * @param w0 variable shift
   * @param y0 scaling coefficient
   * @param alpha exponential coefficient
   * @param beta  exponential coefficient
   */
  GEOS_HOST_DEVICE
  void setCoefficients( T x0, T w0, T y0, T alpha, T beta )
  {
    m_x0 = x0;
    m_w0 = w0;
    m_y0 = y0;
    m_alpha = alpha;
    m_beta = beta;
  }

  // *** no-derivative computes ***

  /**
   * @brief Compute value
   * @param x
   * @param y
   */
  GEOS_HOST_DEVICE
  void compute( const T & x, T & y ) const
  {
    detail::ExponentialCompute< T, EAT >::compute( m_x0, m_y0, m_alpha, x, y );
  }

  /**
   * @brief Compute value
   * @param x
   * @param w
   * @param y
   */
  GEOS_HOST_DEVICE
  void compute( const T & x, const T & w, T & y ) const
  {
    detail::ExponentialCompute< T, EAT >::compute( m_x0, m_w0, m_y0, m_alpha, m_beta, x, w, y );
  }

  /**
   * @brief Compute inverse value
   * @param y
   * @param x
   */
  GEOS_HOST_DEVICE
  void inverse( const T & y, T & x ) const
  {
    detail::ExponentialCompute< T, EAT, true >::compute( m_x0, m_y0, m_alpha, y, x );
  }

  // *** derivative computes ***

  /**
   * @brief Compute value and derivative
   * @param x
   * @param y
   * @param dy_dx
   */
  GEOS_HOST_DEVICE
  void compute( const T & x, T & y, T & dy_dx ) const
  {
    detail::ExponentialCompute< T, EAT >::compute( m_x0, m_y0, m_alpha, x, y, dy_dx );
  }

  /**
   * @brief Compute value and derivative
   * @param x
   * @param w
   * @param y
   * @param dy_dx
   * @param dy_dw
   */
  GEOS_HOST_DEVICE
  void compute( const T & x, const T & w, T & y, T & dy_dx, T & dy_dw ) const
  {
    detail::ExponentialCompute< T, EAT >::compute( m_x0, m_w0, m_y0, m_alpha, m_beta, x, w, y, dy_dx, dy_dw );
  }

  /**
   * @brief Compute inverse value and derivative
   * @param y
   * @param x
   * @param dx_dy
   */
  GEOS_HOST_DEVICE
  void inverse( const T & y, T & x, T & dx_dy ) const
  {
    detail::ExponentialCompute< T, EAT, true >::compute( m_x0, m_y0, m_alpha, y, x, dx_dy );
  }

private:

  /// variable shift
  T m_x0;

  /// variable shift
  T m_w0;

  /// scaling coefficient
  T m_y0;

  /// exponential coefficient
  T m_alpha;

  /// exponential coefficient
  T m_beta;

};



template< ExponentApproximationType EAT >
struct ExponentApproximationTypeWrapper
{
  static constexpr ExponentApproximationType value = EAT;
};

/**
 * @brief Function that implements the switchyard for ExponentApproximationType
 * @tparam LAMBDA the type of the user-provided generic lambda
 * @param type type/order of the approximation (linear, quadratic, full)
 * @param lambda user-provided generic lambda to be called with an instance of EAT wrapper type
 */
template< typename T, typename LAMBDA >
void ExponentApproximationTypeSwitchBlock( ExponentApproximationType const type, T const & x0, T const & y0, T const & alpha, LAMBDA lambda )
{
  switch( type )
  {
    case ExponentApproximationType::Full:
    {
      return lambda( ExponentialRelation< T, ExponentApproximationType::Full >( x0, y0, alpha ) );
    }
    case ExponentApproximationType::Quadratic:
    {
      return lambda( ExponentialRelation< T, ExponentApproximationType::Quadratic >( x0, y0, alpha ) );
    }
    case ExponentApproximationType::Linear:
    {
      return lambda( ExponentialRelation< T, ExponentApproximationType::Linear >( x0, y0, alpha ) );
    }
    default:
    {
      GEOS_ERROR( "ExponentApproximationTypeSwitchBlock() ExponentApproximationType is invalid!" );
    }
  }
}

template< typename T, typename LAMBDA >
void ExponentApproximationTypeSwitchBlock( ExponentApproximationType const type, T const & x0, T const & w0, T const & y0, T const & alpha, T const & beta, LAMBDA lambda )
{
  switch( type )
  {
    case ExponentApproximationType::Full:
    {
      return lambda( ExponentialRelation< T, ExponentApproximationType::Full, 3 >( x0, w0, y0, alpha, beta ) );
    }
    case ExponentApproximationType::Quadratic:
    {
      return lambda( ExponentialRelation< T, ExponentApproximationType::Quadratic, 3 >( x0, w0, y0, alpha, beta ) );
    }
    case ExponentApproximationType::Linear:
    {
      return lambda( ExponentialRelation< T, ExponentApproximationType::Linear, 3 >( x0, w0, y0, alpha, beta ) );
    }
    default:
    {
      GEOS_ERROR( "ExponentApproximationTypeSwitchBlock() ExponentApproximationType is invalid!" );
    }
  }
}

/**
 * @brief Same as above, but with non-default constructor arguments (relation coefficients)
 * @tparam T scalar real-valued type
 * @tparam LAMBDA type of the user-provided generic lambda
 * @param type runtime choice of ExponentApproximationType
 * @param x0 variable shift
 * @param y0 scaling coefficient
 * @param alpha exponential coefficient
 * @param lambda user-provided generic lambda that will be passed an instance of ExponentialRelation<T,?>
 */
// For two variables
template< typename T, typename LAMBDA >
void makeExponentialRelation( ExponentApproximationType type,
                              T const & x0, T const & y0, T const & alpha,
                              LAMBDA && lambda )
{
  ExponentApproximationTypeSwitchBlock( type, x0, y0, alpha, std::forward< LAMBDA >( lambda ) );
}

// For three variables
template< typename T, typename LAMBDA >
void makeExponentialRelation( ExponentApproximationType type,
                              T const & x0, T const & w0, T const & y0, T const & alpha, T const & beta,
                              LAMBDA && lambda )
{
  ExponentApproximationTypeSwitchBlock( type, x0, w0, y0, alpha, beta, std::forward< LAMBDA >( lambda ) );
}

} // namespace constitutive

} // namespace geos

#endif // GEOS_CONSITUTIVE_EXPONENTIALRELATION_HPP_
