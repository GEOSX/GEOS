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
 * @file InterfaceTypes.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_INTERFACETYPES_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_INTERFACETYPES_HPP_

#include "common/GeosxConfig.hpp"

#ifdef GEOS_USE_TRILINOS
#include "linearAlgebra/interfaces/trilinos/TrilinosInterface.hpp"
#endif

#ifdef GEOS_USE_HYPRE
#include "linearAlgebra/interfaces/hypre/HypreInterface.hpp"
#endif

#ifdef GEOS_USE_PETSC
#include "linearAlgebra/interfaces/petsc/PetscInterface.hpp"
#endif

namespace geos
{

/// Alias for current interface
using LAInterface = GEOS_CONCAT( GEOS_LA_INTERFACE, Interface );

/// Alias for ParallelMatrix
using ParallelMatrix = LAInterface::ParallelMatrix;

/// Alias for ParallelVector
using ParallelVector = LAInterface::ParallelVector;

/**
 * @brief High level interface to call the proper LAI initialize function.
 */
inline void setupLAI()
{
#ifdef GEOS_USE_TRILINOS
  TrilinosInterface::initialize();
#endif
#ifdef GEOS_USE_HYPRE
  HypreInterface::initialize();
#endif
#ifdef GEOS_USE_PETSC
  PetscInterface::initialize();
#endif
}

/**
 * @brief High level interface to call the proper LAI finalize function.
 */
inline void finalizeLAI()
{
#ifdef GEOS_USE_TRILINOS
  TrilinosInterface::finalize();
#endif
#ifdef GEOS_USE_HYPRE
  HypreInterface::finalize();
#endif
#ifdef GEOS_USE_PETSC
  PetscInterface::finalize();
#endif
}

}



#endif /*GEOS_LINEARALGEBRA_INTERFACES_INTERFACETYPES_HPP_*/
