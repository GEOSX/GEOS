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

#include "common/GeosxConfig.hpp"
#include "common/logger/Logger.hpp"
#include "mainInterface/GeosxVersion.hpp"

namespace geos
{

std::string getVersion()
{
#if defined(GEOS_GIT_BRANCH) && defined(GEOS_GIT_HASH)
  return GEOS_VERSION_FULL " (" GEOS_GIT_BRANCH ", sha1: " GEOS_GIT_HASH ")";
#else
  return GEOS_VERSION_FULL;
#endif
}

static std::string getCppCompilerIdString()
{
  std::ostringstream oss;

  #if defined(__clang__)
  oss << "clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__;
#if defined(__apple_build_version__)
  oss << " (apple version " << __apple_build_version__ << ")";
#endif
#if defined(__ibmxl_vrm__)
  oss << "IBMXL "
      << __ibmxl_version__ << "."
      << __ibmxl_release__ << "."
      << __ibmxl_modification__ << "."
      << __ibmxl_ptf_fix_level__;
#endif
#elif defined(__GNUC__)
  oss << "gcc "
      << __GNUC__ << "."
      << __GNUC_MINOR__ << "."
      << __GNUC_PATCHLEVEL__;
#endif
  return oss.str();
}

static std::string getGpuCompilerIdString()
{
  std::ostringstream oss;

#if defined( GEOS_USE_CUDA )
  oss << "  - CUDA compiler version: " << CUDA_VERSION/10/100 << "." << CUDA_VERSION/10%100;
#endif
#if defined( GEOS_USE_HIP )
  oss << "  - ROCm compiler version: " << ROCM_VERSION/100/100 << "." << ROCM_VERSION/100%100;
#endif
  return oss.str();
}

void outputVersionInfo()
{

  GEOS_LOG_RANK_0( "GEOS version: " << getVersion() );

  GEOS_LOG_RANK_0( "  - c++ compiler: " << getCppCompilerIdString() );

  std::string const gpuCompilerIdString = getGpuCompilerIdString();
  GEOS_LOG_RANK_0_IF( !gpuCompilerIdString.empty(), gpuCompilerIdString );

#if defined(_OPENMP)
  GEOS_LOG_RANK_0( "  - openmp version: " << _OPENMP );
#endif

#if defined(GEOS_USE_MPI)
  {
    char version[MPI_MAX_LIBRARY_VERSION_STRING];
    int len;
    MPI_Get_library_version( version, &len );
    GEOS_LOG_RANK_0( "  - MPI version: " << version );
  }
#endif

#if defined(HDF5_VERSION)
  GEOS_LOG_RANK_0( "  - HDF5 version: " << STRINGIZE( HDF5_VERSION ) );
#endif

#if defined(Conduit_VERSION)
  GEOS_LOG_RANK_0( "  - Conduit version: " << STRINGIZE( Conduit_VERSION ) );
#endif

#if defined(VTK_VERSION)
  GEOS_LOG_RANK_0( "  - VTK version: " << STRINGIZE( VTK_VERSION ) );
#endif

#if defined(RAJA_VERSION)
  GEOS_LOG_RANK_0( "  - RAJA version: " << STRINGIZE( RAJA_VERSION ) );
#endif

#if defined(umpire_VERSION)
  GEOS_LOG_RANK_0( "  - umpire version: " << STRINGIZE( umpire_VERSION ) );
#endif

#if defined(chai_VERSION)
  GEOS_LOG_RANK_0( "  - chai version: " << STRINGIZE( chai_VERSION ) );
#endif

#if defined(adiak_VERSION)
  GEOS_LOG_RANK_0( "  - adiak version: " << STRINGIZE( adiak_VERSION ) );
#endif

#if defined(caliper_VERSION)
  GEOS_LOG_RANK_0( "  - caliper version: " << STRINGIZE( caliper_VERSION ) );
#endif

#if defined(metis_VERSION)
  GEOS_LOG_RANK_0( "  - METIS version: " << STRINGIZE( metis_VERSION ) );
#endif

#if defined(parmetis_VERSION)
  GEOS_LOG_RANK_0( "  - PARAMETIS version: " << STRINGIZE( parmetis_VERSION ) );
#endif

#if defined(scotch_VERSION)
  GEOS_LOG_RANK_0( "  - scotch version: " << STRINGIZE( scotch_VERSION ) );
#endif

#if defined(superlu_dist_VERSION)
  GEOS_LOG_RANK_0( "  - superlu_dist version: " << STRINGIZE( superlu_dist_VERSION ) );
#endif

#if defined(suitesparse_VERSION)
  GEOS_LOG_RANK_0( "  - suitesparse version: " << STRINGIZE( suitesparse_VERSION ) );
#endif

#if defined(hypre_VERSION)
  GEOS_LOG_RANK_0( "  - hypre version: " << STRINGIZE( hypre_VERSION ) );
#endif

#if defined(trilinos_VERSION)
  GEOS_LOG_RANK_0( "  - trilinos version: " << STRINGIZE( trilinos_VERSION ) );
#endif

#if defined(petsc_VERSION)
  GEOS_LOG_RANK_0( "  - petsc version: " << STRINGIZE( petsc_VERSION ) );
#endif

#if defined(Python3_VERSION)
  GEOS_LOG_RANK_0( "  - Python3 version: " << STRINGIZE( Python3_VERSION ) );
#endif

#if defined(CUDAToolkit_VERSION)
  GEOS_LOG_RANK_0( "  - CUDAToolkit version: " << STRINGIZE( CUDAToolkit_VERSION ) );
#endif

#if \
  defined(GEOS_USE_DEVICE) && \
  defined(GEOS_USE_HYPRE) && \
  ( GEOS_USE_HYPRE_DEVICE == GEOS_USE_HYPRE_CPU )
  GEOS_LOG_RANK_0( "" );
  GEOS_LOG_RANK_0( "**************************************************" );
  GEOS_LOG_RANK_0( "*                   WARNING!!!                   *" );
  GEOS_LOG_RANK_0( "*                                                *" );
  GEOS_LOG_RANK_0( "*  GEOS has GPU support enabled, but not HYPRE!  *" );
  GEOS_LOG_RANK_0( "**************************************************" );
  GEOS_LOG_RANK_0( "" );
#endif

}

}
