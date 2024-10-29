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

/**
 * @file common.hpp
 */
#ifndef GEOS_LINEARALGEBRA_INTERFACES_COMMON_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_COMMON_HPP_

#include "common/DataTypes.hpp"

/**
 * Whether to check preconditions at runtime in LAI functions
 */
#define GEOS_LAI_RUNTIME_ASSERT 1
#if GEOS_LAI_RUNTIME_ASSERT

/**
 * Assert expression is true
 * @param expr expression
 */
#define GEOS_LAI_ASSERT( expr ) GEOS_ERROR_IF( !(expr), "" )

/**
 * Assert expression and output message if false
 * @param expr expression
 * @param msg message
 */
#define GEOS_LAI_ASSERT_MSG( expr, msg ) GEOS_ERROR_IF( !(expr), msg )

/**
 * Assert lhs equals rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
#define GEOS_LAI_ASSERT_EQ( lhs, rhs ) GEOS_ERROR_IF_NE( lhs, rhs )

/**
 * Assert lhs not equal to rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
#define GEOS_LAI_ASSERT_NE( lhs, rhs ) GEOS_ERROR_IF_EQ( lhs, rhs )

/**
 * Assert lhs greater than rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
#define GEOS_LAI_ASSERT_GT( lhs, rhs ) GEOS_ERROR_IF_GE( rhs, lhs )

/**
 * Assert lhs greater than or equal to rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
#define GEOS_LAI_ASSERT_GE( lhs, rhs ) GEOS_ERROR_IF_GT( rhs, lhs )
#else

/**
 * Assert expression is true
 * @param expr expression
 */
#define GEOS_LAI_ASSERT( expr ) GEOS_ASSERT( expr )

/**
 * Assert expression and output message if false
 * @param expr expression
 * @param msg message
 */
#define GEOS_LAI_ASSERT_MSG( expr, msg ) GEOS_ASSERT_MSG( expr, msg )

/**
 * Assert lhs equals rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
#define GEOS_LAI_ASSERT_EQ( lhs, rhs ) GEOS_ASSERT_EQ( lhs, rhs )

/**
 * Assert lhs not equal to rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
#define GEOS_LAI_ASSERT_NE( lhs, rhs ) GEOS_ASSERT_NE( lhs, rhs )

/**
 * Assert lhs greater rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
#define GEOS_LAI_ASSERT_GT( lhs, rhs ) GEOS_ASSERT_GT( lhs, rhs )

/**
 * Assert lhs greater than or equal to rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
#define GEOS_LAI_ASSERT_GE( lhs, rhs ) GEOS_ASSERT_GE( lhs, rhs )
#endif

/**
 * Macro for checking and reporting error codes from TPL packages
 * @param call call to check function
 */
#define GEOS_LAI_CHECK_ERROR( call ) \
  do { \
    auto const ierr = call; \
    GEOS_ERROR_IF_NE_MSG( ierr, 0, "Error in call to " << #call ); \
  } while( false )

/**
 * Macro for checking and reporting non-negative error codes from TPL packages
 * @param call call to check function
 */
#define GEOS_LAI_CHECK_ERROR_NNEG( call ) \
  do { \
    auto const ierr = call; \
    GEOS_ERROR_IF_GT_MSG( 0, ierr, "Error in call to " << #call ); \
  } while( false )

namespace geos
{

/**
 *  Enumeration of available output formats for LAI objects
 */
enum class LAIOutputFormat
{
  NATIVE_ASCII,
  NATIVE_BINARY,
  MATLAB_ASCII,
  MATLAB_BINARY,
  MATRIX_MARKET
};

}

#endif //GEOS_LINEARALGEBRA_INTERFACES_COMMON_HPP_
