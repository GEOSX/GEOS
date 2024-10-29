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
 * @file HypreUtils.cpp
 */

#include "HypreUtils.hpp"

#include "linearAlgebra/interfaces/hypre/HypreVector.hpp"

#include <_hypre_parcsr_mv.h>
#include <_hypre_parcsr_ls.h>

namespace geos
{

namespace hypre
{

HYPRE_Vector parVectorToVectorAll( HYPRE_ParVector const vec )
{
  if( hypre_ParVectorMemoryLocation( vec ) == HYPRE_MEMORY_HOST )
  {
    return (HYPRE_Vector)hypre_ParVectorToVectorAll( vec );
  }
  else
  {
    // Explicitly copy data to host before gathering, since hypre_ParVectorToVectorAll relies on UM.
    hypre_ParVector * const hostVec = hypre_ParVectorCloneDeep_v2( vec, HYPRE_MEMORY_HOST );
    hypre_Vector * const fullVec = hypre_ParVectorToVectorAll( hostVec );
    hypre_ParVectorDestroy( hostVec );
    return (HYPRE_Vector)fullVec;
  }
}

HYPRE_Int dummySetup( HYPRE_Solver,
                      HYPRE_ParCSRMatrix,
                      HYPRE_ParVector,
                      HYPRE_ParVector )
{
  return 0;
}

HYPRE_Int SuperLUDistSolve( HYPRE_Solver solver,
                            HYPRE_ParCSRMatrix A,
                            HYPRE_ParVector b,
                            HYPRE_ParVector x )
{
  GEOS_UNUSED_VAR( A );
#if defined(GEOS_USE_SUPERLU_DIST)
  return hypre_SLUDistSolve( solver, b, x );
#else
  GEOS_UNUSED_VAR( solver );
  GEOS_UNUSED_VAR( b );
  GEOS_UNUSED_VAR( x );
  GEOS_ERROR( "GEOSX is configured without support for SuperLU_dist." );
  return -1;
#endif
}

HYPRE_Int SuperLUDistDestroy( HYPRE_Solver solver )
{
#if defined(GEOS_USE_SUPERLU_DIST)
  return hypre_SLUDistDestroy( solver );
#else
  GEOS_UNUSED_VAR( solver );
  GEOS_ERROR( "GEOSX is configured without support for SuperLU_dist." );
  return -1;
#endif
}

/**
 * @brief Holds temporary data needed by hypre's relaxation methods.
 *
 * Many implementations dispatched to by hypre_BoomerAMGRelax crash unless Vtemp,
 * Ztemp (temporary vectors) or l1Norms (row norm vector) parameters are provided.
 * We want these temporary vectors to persist across multiple applications of the relaxation method.
 * In order to make it conform to hypre's solver interface, we allocate these vectors during setup
 * and disguise them behind HYPRE_Solver, hypre's opaque solver struct pointer type.
 * Hence multiple reinterpret_casts (replacing C-style casts used in hypre) are required
 * to convert back and forth between HYPRE_Solver and pointer to actual data struct.
 * We also have to manage memory with raw new/delete here, because a single pointer
 * that is maintained by hypre is all we got to manage the struct with.
 */
struct RelaxationData
{
  HYPRE_Int type = -1;
  HypreVector vtemp;
  HypreVector ztemp;
  HYPRE_Real * norms{};
};

HYPRE_Int relaxationCreate( HYPRE_Solver & solver,
                            HYPRE_Int const type )
{
  RelaxationData * const data = new RelaxationData;
  data->type = type;
  solver = reinterpret_cast< HYPRE_Solver >( data );
  return 0;
}

HYPRE_Int getL1NormType( HYPRE_Int const relaxType )
{
  // Extracted from the depths of hypre_BoomerAMGSetup()
  switch( relaxType )
  {
    case 7: return 5;
    case 8:
    case 13:
    case 14: return 4;
    case 18: return 1;
    default: return -1;
  }
}

HYPRE_Int relaxationSetup( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x )
{
  GEOS_UNUSED_VAR( b, x );

  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  data->vtemp.create( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );
  data->ztemp.create( hypre_ParCSRMatrixNumRows( A ), hypre_ParCSRMatrixComm( A ) );

  // L1-smoothers need row L1-norms precomputed
  HYPRE_Int const normType = getL1NormType( data->type );
  if( normType > 0 )
  {
    if( data->norms )
    {
      hypre_TFree( data->norms, hypre::memoryLocation );
    }
    hypre_ParCSRComputeL1Norms( A, normType, nullptr, &data->norms );
  }
  return 0;
}

HYPRE_Int relaxationSolve( HYPRE_Solver solver,
                           HYPRE_ParCSRMatrix A,
                           HYPRE_ParVector b,
                           HYPRE_ParVector x )
{
  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  return hypre_BoomerAMGRelax( A, b, nullptr, data->type, 0, 1.0, 1.0, data->norms, x, data->vtemp.unwrapped(), data->ztemp.unwrapped() );
}

HYPRE_Int relaxationDestroy( HYPRE_Solver solver )
{
  // Refer to RelaxationData doxygen above for explanation of reinterpret_cast
  RelaxationData * const data = reinterpret_cast< RelaxationData * >( solver );
  if( data->norms )
  {
    hypre_TFree( data->norms, hypre::memoryLocation );
  }
  delete data;
  return 0;
}

} // namespace hypre

} // namespace geos
