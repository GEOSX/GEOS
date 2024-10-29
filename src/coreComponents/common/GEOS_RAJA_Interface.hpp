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

#ifndef GEOS_RAJAINTERFACE_RAJAINTERFACE_HPP
#define GEOS_RAJAINTERFACE_RAJAINTERFACE_HPP

// Source includes
#include "common/GeosxConfig.hpp"
#include "common/DataTypes.hpp"

// TPL includes
#include <RAJA/RAJA.hpp>

#include <chrono>
#include <thread>

#define GEOS_ASYNC_WAIT( UPPER, NANOSLEEP, TEST ) while( !TEST ) { }

namespace geos
{

auto const hostMemorySpace = LvArray::MemorySpace::host;

using serialPolicy = RAJA::seq_exec;
using serialAtomic = RAJA::seq_atomic;
using serialReduce = RAJA::seq_reduce;

using serialStream = RAJA::resources::Host;
using serialEvent = RAJA::resources::HostEvent;

#if defined( GEOS_USE_OPENMP )

auto const parallelHostMemorySpace = hostMemorySpace;

using parallelHostPolicy = RAJA::omp_parallel_for_exec;
using parallelHostReduce = RAJA::omp_reduce;
using parallelHostAtomic = RAJA::builtin_atomic;

// issues with Raja::resources::Omp on lassen
using parallelHostStream = serialStream;
using parallelHostEvent = serialEvent;

void RAJA_INLINE parallelHostSync() { RAJA::synchronize< RAJA::omp_synchronize >(); }

#else

auto const parallelHostMemorySpace = hostMemorySpace;

using parallelHostPolicy = serialPolicy;
using parallelHostReduce = serialReduce;
using parallelHostAtomic = serialAtomic;
using parallelHostStream = serialStream;
using parallelHostEvent = serialEvent;

void RAJA_INLINE parallelHostSync() { }

#endif

#if defined( GEOS_USE_CUDA )
auto const parallelDeviceMemorySpace = LvArray::MemorySpace::cuda;

template< size_t BLOCK_SIZE = GEOS_BLOCK_SIZE >
using parallelDevicePolicy = RAJA::cuda_exec< BLOCK_SIZE >;

template< size_t BLOCK_SIZE = GEOS_BLOCK_SIZE >
using parallelDeviceAsyncPolicy = RAJA::cuda_exec_async< BLOCK_SIZE >;

using parallelDeviceStream = RAJA::resources::Cuda;
using parallelDeviceEvent = RAJA::resources::Event;

using parallelDeviceReduce = RAJA::cuda_reduce;
using parallelDeviceAtomic = RAJA::cuda_atomic;

void RAJA_INLINE parallelDeviceSync() { RAJA::synchronize< RAJA::cuda_synchronize >(); }

template< typename POLICY, typename RESOURCE, typename LAMBDA >
RAJA_INLINE parallelDeviceEvent forAll( RESOURCE && stream, const localIndex end, LAMBDA && body )
{
  return RAJA::forall< POLICY >( std::forward< RESOURCE >( stream ),
                                 RAJA::TypedRangeSegment< localIndex >( 0, end ),
                                 std::forward< LAMBDA >( body ) );
}

#elif defined( GEOS_USE_HIP )

auto const parallelDeviceMemorySpace = LvArray::MemorySpace::hip;

template< size_t BLOCK_SIZE = GEOS_BLOCK_SIZE >
using parallelDevicePolicy = RAJA::hip_exec< BLOCK_SIZE >;


using parallelDeviceStream = RAJA::resources::Hip;
using parallelDeviceEvent = RAJA::resources::Event;

using parallelDeviceReduce = RAJA::hip_reduce;
using parallelDeviceAtomic = RAJA::hip_atomic;

void RAJA_INLINE parallelDeviceSync() { RAJA::synchronize< RAJA::hip_synchronize >( ); }

// the async dispatch policy caused runtime issues as of rocm@4.5.2, hasn't been checked in rocm@5:
template< size_t BLOCK_SIZE = GEOS_BLOCK_SIZE >
using parallelDeviceAsyncPolicy = parallelDevicePolicy< BLOCK_SIZE >; // RAJA::hip_exec_async< BLOCK_SIZE >;

template< typename POLICY, typename RESOURCE, typename LAMBDA >
RAJA_INLINE parallelDeviceEvent forAll( RESOURCE && GEOS_UNUSED_PARAM( stream ), const localIndex end, LAMBDA && body )
{
  RAJA::forall< POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, end ), std::forward< LAMBDA >( body ) );
  return parallelDeviceEvent();
}
#else

auto const parallelDeviceMemorySpace = parallelHostMemorySpace;

template< size_t BLOCK_SIZE = 0 >
using parallelDevicePolicy = parallelHostPolicy;

template< size_t BLOCK_SIZE = 0 >
using parallelDeviceAsyncPolicy = parallelHostPolicy;

using parallelDeviceStream = parallelHostStream;
using parallelDeviceEvent = parallelHostEvent;

using parallelDeviceReduce = parallelHostReduce;
using parallelDeviceAtomic = parallelHostAtomic;

void RAJA_INLINE parallelDeviceSync() { parallelHostSync( ); }

template< typename POLICY, typename RESOURCE, typename LAMBDA >
RAJA_INLINE parallelDeviceEvent forAll( RESOURCE && GEOS_UNUSED_PARAM( stream ), const localIndex end, LAMBDA && body )
{
  RAJA::forall< POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, end ), std::forward< LAMBDA >( body ) );
  return parallelDeviceEvent();
}

#endif

using parallelDeviceEvents = std::vector< parallelDeviceEvent >;

namespace internalRajaInterface
{
template< typename >
struct PolicyMap
{};

template<>
struct PolicyMap< serialPolicy >
{
  using atomic = serialAtomic;
  using reduce = serialReduce;
};

#if defined(GEOS_USE_OPENMP)
template<>
struct PolicyMap< RAJA::omp_parallel_for_exec >
{
  using atomic = RAJA::builtin_atomic;
  using reduce = RAJA::omp_reduce;
};
#endif

#if defined(GEOS_USE_CUDA)
template< typename X, typename Y, typename C, size_t BLOCK_SIZE, bool ASYNC >
struct PolicyMap< RAJA::policy::cuda::cuda_exec_explicit< X, Y, C, BLOCK_SIZE, ASYNC > >
{
  using atomic = RAJA::cuda_atomic;
  using reduce = RAJA::cuda_reduce;
};
#endif

#if defined(GEOS_USE_HIP)
template< size_t BLOCK_SIZE, bool ASYNC >
struct PolicyMap< RAJA::hip_exec< BLOCK_SIZE, ASYNC > >
{
  using atomic = RAJA::hip_atomic;
  using reduce = RAJA::hip_reduce;
};
#endif
}


template< typename POLICY >
using ReducePolicy = typename internalRajaInterface::PolicyMap< POLICY >::reduce;

template< typename POLICY >
using AtomicPolicy = typename internalRajaInterface::PolicyMap< POLICY >::atomic;


RAJA_INLINE bool testAllDeviceEvents( parallelDeviceEvents & events )
{
  bool allDone = true;
  for( auto & event : events )
  {
    if( !event.check() )
    {
      allDone = false;
      break;
    }
  }
  return allDone;
}

RAJA_INLINE void waitAllDeviceEvents( parallelDeviceEvents & events )
{
  // poll device events for completion then wait 10 nanoseconds 6,000,000,000 times (60 sec timeout)
  // 10 nsecs ~= 30 clock cycles @ 3Ghz
  GEOS_ASYNC_WAIT( 6000000000, 10, testAllDeviceEvents( events ) );
}

template< typename POLICY, typename INDEX, typename LAMBDA >
RAJA_INLINE void forAll( INDEX const end, LAMBDA && body )
{
  RAJA::forall< POLICY >( RAJA::TypedRangeSegment< INDEX >( 0, end ), std::forward< LAMBDA >( body ) );
}

template< typename POLICY, typename INDEX, typename LAMBDA >
RAJA_INLINE void forRange( INDEX const begin, INDEX const end, LAMBDA && body )
{
  RAJA::forall< POLICY >( RAJA::TypedRangeSegment< INDEX >( begin, end ), std::forward< LAMBDA >( body ) );
}

} // namespace geos

#endif // GEOS_RAJAINTERFACE_RAJAINTERFACE_HPP
