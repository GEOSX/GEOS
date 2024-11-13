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
 * @file HypreInterface.cpp
 */

#include "HypreInterface.hpp"

#include "linearAlgebra/interfaces/direct/SuiteSparse.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/hypre/HyprePreconditioner.hpp"
#include "linearAlgebra/interfaces/hypre/HypreSolver.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"

#if defined(GEOS_USE_SUPERLU_DIST)
#include "linearAlgebra/interfaces/direct/SuperLUDist.hpp"
#endif

#include "HYPRE_utilities.h"
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
#include "_hypre_utilities.h"
#include "_hypre_utilities.hpp"
#endif

namespace geos
{

void HypreInterface::initialize()
{
#ifdef GEOS_USE_OPENMP
  GEOS_LOG_RANK_0_IF( omp_get_max_threads()>1,
                      "OMP_NUM_THREADS > 1 may not be optimal for certain hypre preconditioning options. " );
#endif

  HYPRE_Init();
#if GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CUDA || GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_HIP
  HYPRE_SetExecutionPolicy( HYPRE_EXEC_DEVICE );
  HYPRE_SetSpGemmUseVendor( 0 );
  HYPRE_DeviceInitialize();
#endif
  HYPRE_SetMemoryLocation( hypre::memoryLocation );
  HYPRE_SetPrintErrorMode( 1 );
}

void HypreInterface::finalize()
{
  HYPRE_Finalize();
}

std::unique_ptr< LinearSolverBase< HypreInterface > >
HypreInterface::createSolver( LinearSolverParameters params )
{
  if( params.solverType == LinearSolverParameters::SolverType::direct )
  {
    if( params.direct.parallel )
    {
#if defined(GEOS_USE_SUPERLU_DIST)
      return std::make_unique< SuperLUDist< HypreInterface > >( std::move( params ) );
#else
      GEOS_ERROR( "GEOSX is configured without support for SuperLU_dist." );
      return std::unique_ptr< LinearSolverBase< HypreInterface > >( nullptr );
#endif
    }
    else
    {
      return std::make_unique< SuiteSparse< HypreInterface > >( std::move( params ) );
    }
  }
  else
  {
    return std::make_unique< HypreSolver >( std::move( params ) );
  }
}

std::unique_ptr< PreconditionerBase< HypreInterface > >
geos::HypreInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< HyprePreconditioner >( std::move( params ) );
}

std::unique_ptr< PreconditionerBase< HypreInterface > >
geos::HypreInterface::createPreconditioner( LinearSolverParameters params,
                                            array1d< HypreVector > const & nearNullKernel )
{
  return std::make_unique< HyprePreconditioner >( std::move( params ), nearNullKernel );
}

}
