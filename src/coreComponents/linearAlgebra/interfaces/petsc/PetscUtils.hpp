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
 * @file PetscUtils.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_PETSCUTILS_HPP
#define GEOS_LINEARALGEBRA_INTERFACES_PETSCUTILS_HPP

#include "common/DataTypes.hpp"

#include <petscsys.h>

namespace geos
{

namespace petsc
{

// Check matching requirements on index/value types between GEOSX and PETSc

static_assert( sizeof( PetscInt ) == sizeof( geos::globalIndex ),
               "PetscInt and geos::globalIndex must have the same size" );

static_assert( std::is_signed< PetscInt >::value == std::is_signed< geos::globalIndex >::value,
               "PetscInt and geoex::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< PetscScalar, geos::real64 >::value,
               "PetscScalar and geos::real64 must be the same type" );

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to PetscInt
 * @param[in] index the input array
 * @return the converted array
 */
inline PetscInt * toPetscInt( globalIndex * const index )
{
  return reinterpret_cast< PetscInt * >(index);
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to PetscInt
 * @param[in] index the input array
 * @return the converted array
 */
inline PetscInt const * toPetscInt( globalIndex const * const index )
{
  return reinterpret_cast< PetscInt const * >(index);
}

} // namespace petsc

} // namespace geos

#endif //GEOS_LINEARALGEBRA_INTERFACES_PETSCUTILS_HPP
