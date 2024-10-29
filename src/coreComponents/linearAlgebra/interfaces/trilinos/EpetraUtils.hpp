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
 * @file EpetraUtils.hpp
 */
#ifndef GEOS_LINEARALGEBRA_INTERFACES_EPETRAUTILS_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_EPETRAUTILS_HPP_

#include "common/DataTypes.hpp"

#ifdef GEOS_USE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

namespace geos
{

namespace trilinos
{

// Check matching requirements on index/value types between GEOSX and PETSc

static_assert( sizeof( long long ) == sizeof( geos::globalIndex ),
               "long long and geos::globalIndex must have the same size" );

static_assert( std::is_signed< long long >::value == std::is_signed< geos::globalIndex >::value,
               "long long and geos::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, geos::real64 >::value,
               "double and geos::real64 must be the same type" );

#ifdef GEOS_USE_MPI
/// Alias for Epetra communicator
using EpetraComm = Epetra_MpiComm;
#else
/// Serial fallback alias for Epetra communicator
using EpetraComm = Epetra_SerialComm;
#endif

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to Epetra long long
 * @param[in] index the input array
 * @return the converted array
 */
inline long long * toEpetraLongLong( globalIndex * const index )
{
  return reinterpret_cast< long long * >(index);
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to Epetra long long
 * @param[in] index the input array
 * @return the converted array
 */
inline long long const * toEpetraLongLong( globalIndex const * const index )
{
  return reinterpret_cast< long long const * >(index);
}

} // namespace trilinos

} // namespace geos

#endif //GEOS_LINEARALGEBRA_INTERFACES_EPETRAUTILS_HPP_
