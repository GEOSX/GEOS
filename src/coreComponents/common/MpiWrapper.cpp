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
 * @file MpiWrapper.cpp
 */

#include "MpiWrapper.hpp"
#include <unistd.h>

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif


namespace geos
{

void MpiWrapper::barrier( MPI_Comm const & MPI_PARAM( comm ) )
{
#ifdef GEOS_USE_MPI
  MPI_Barrier( comm );
#endif
}

int MpiWrapper::cartCoords( MPI_Comm comm, int rank, int maxdims, int coords[] )
{
#ifdef GEOS_USE_MPI
  return MPI_Cart_coords( comm, rank, maxdims, coords );
#else
  return 0;
#endif
}

int MpiWrapper::cartCreate( MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                            int reorder, MPI_Comm * comm_cart )
{
#ifdef GEOS_USE_MPI
  return MPI_Cart_create( comm_old, ndims, dims, periods, reorder, comm_cart );
#else
  return 0;
#endif
}

int MpiWrapper::cartRank( MPI_Comm comm, const int coords[] )
{
  int rank = 0;
#ifdef GEOS_USE_MPI
  MPI_Cart_rank( comm, coords, &rank );
#endif
  return rank;
}

void MpiWrapper::commFree( MPI_Comm & comm )
{
#ifdef GEOS_USE_MPI
  MPI_CHECK_ERROR( MPI_Comm_free( &comm ) );
#else
//  comm = MPI_COMM_NULL;
#endif
}

int MpiWrapper::commRank( MPI_Comm const & MPI_PARAM( comm ) )
{
  int rank = 0;
#ifdef GEOS_USE_MPI
  MPI_Comm_rank( comm, &rank );
#endif
  return rank;
}

int MpiWrapper::commSize( MPI_Comm const & MPI_PARAM( comm ) )
{
  int size = 1;
#ifdef GEOS_USE_MPI
  MPI_Comm_size( comm, &size );
#endif
  return size;
}

bool MpiWrapper::commCompare( MPI_Comm const & comm1, MPI_Comm const & comm2 )
{
#ifdef GEOS_USE_MPI
  int result;
  MPI_Comm_compare( comm1, comm2, &result );
  return result == MPI_IDENT || result == MPI_CONGRUENT;
#else
  return comm1 == comm2;
#endif
}

bool MpiWrapper::initialized()
{
#ifdef GEOS_USE_MPI
  int ret = false;
  MPI_CHECK_ERROR( MPI_Initialized( &ret ) );
  return ret;
#else
  return false;
#endif
}

int MpiWrapper::init( int * argc, char * * * argv )
{
#ifdef GEOS_USE_MPI
  return MPI_Init( argc, argv );
#else
  return 0;
#endif
}

void MpiWrapper::finalize()
{
#ifdef GEOS_USE_MPI
  MPI_CHECK_ERROR( MPI_Finalize() );
#endif
}


MPI_Comm MpiWrapper::commDup( MPI_Comm const comm )
{
#ifdef GEOS_USE_MPI
  MPI_Comm duplicate;
  MPI_CHECK_ERROR( MPI_Comm_dup( comm, &duplicate ) );
  return duplicate;
#else
  return comm;
#endif
}

MPI_Comm MpiWrapper::commSplit( MPI_Comm const comm, int color, int key )
{
#ifdef GEOS_USE_MPI
  MPI_Comm scomm;
  MPI_CHECK_ERROR( MPI_Comm_split( comm, color, key, &scomm ) );
  return scomm;
#else
  return comm;
#endif
}

int MpiWrapper::test( MPI_Request * request, int * flag, MPI_Status * status )
{
#ifdef GEOS_USE_MPI
  return MPI_Test( request, flag, status );
#else
  *flag = 0;
  return 0;
#endif
}

int MpiWrapper::testAny( int count, MPI_Request array_of_requests[], int * idx, int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOS_USE_MPI
  return MPI_Testany( count, array_of_requests, idx, flag, array_of_statuses );
#else
  *flag = 0;
  return 0;
#endif
}

int MpiWrapper::testSome( int count, MPI_Request array_of_requests[], int * outcount, int array_of_indices[], MPI_Status array_of_statuses[] )
{
#ifdef GEOS_USE_MPI
  return MPI_Testsome( count, array_of_requests, outcount, array_of_indices, array_of_statuses );
#else
  *outcount = 0;
  return 0;
#endif
}

int MpiWrapper::testAll( int count, MPI_Request array_of_requests[], int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOS_USE_MPI
  return MPI_Testall( count, array_of_requests, flag, array_of_statuses );
#else
  *flag = 0;
  return 0;
#endif
}

int MpiWrapper::check( MPI_Request * request, int * flag, MPI_Status * status )
{
#ifdef GEOS_USE_MPI
  return MPI_Request_get_status( *request, flag, status );
#else
  *flag = 0;
  return 0;
#endif
}

int MpiWrapper::checkAny( int count, MPI_Request array_of_requests[], int * idx, int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOS_USE_MPI
  bool found = false;
  int flagCache = -1;
  int rval = MPI_SUCCESS;
  std::vector< int > rvals( count );
  for( int jdx = 0; jdx < count; ++jdx )
  {
    *flag = 0;
    rvals[ jdx ] = MPI_Request_get_status( array_of_requests[ jdx ], flag, &array_of_statuses[ jdx ] );
    if( *flag && !found )
    {
      *idx = jdx;
      flagCache = *flag;
    }
    if( rvals[ jdx ] != MPI_SUCCESS )
    {
      rval = rvals[ jdx ];
    }
  }
  if( found )
  {
    *flag = flagCache;
  }
  return rval;
#else
  *flag = 0;
  return 0;
#endif
}

int MpiWrapper::checkAll( int count, MPI_Request array_of_requests[], int * flag, MPI_Status array_of_statuses[] )
{
#ifdef GEOS_USE_MPI
  // assume all passing, any that don't pass set the flag to false
  *flag = 1;
  int rval = MPI_SUCCESS;
  std::vector< int > rvals( count );
  int iFlag = 0;
  for( int idx = 0; idx < count; ++idx )
  {
    rvals[ idx ] = MPI_Request_get_status( array_of_requests[ idx ], &iFlag, &array_of_statuses[ idx ] );
    if( !iFlag )
    {
      *flag = iFlag;
    }
    if( rvals[ idx ] != MPI_SUCCESS )
    {
      rval = rvals[ idx ];
    }
  }
  return rval;
#else
  *flag = 0;
  return 0;
#endif
}

int MpiWrapper::wait( MPI_Request * request, MPI_Status * status )
{
#ifdef GEOS_USE_MPI
  return MPI_Wait( request, status );
#else
  return 0;
#endif
}

int MpiWrapper::waitAny( int count, MPI_Request array_of_requests[], int * indx, MPI_Status array_of_statuses[] )
{
#ifdef GEOS_USE_MPI
  return MPI_Waitany( count, array_of_requests, indx, array_of_statuses );
#else
  return 0;
#endif
}

int MpiWrapper::waitSome( int count, MPI_Request array_of_requests[], int * outcount, int array_of_indices[], MPI_Status array_of_statuses[] )
{
#ifdef GEOS_USE_MPI
  return MPI_Waitsome( count, array_of_requests, outcount, array_of_indices, array_of_statuses );
#else
  // *outcount = 0;
  return 0;
#endif
}

int MpiWrapper::waitAll( int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[] )
{
#ifdef GEOS_USE_MPI
  return MPI_Waitall( count, array_of_requests, array_of_statuses );
#else
  return 0;
#endif
}

double MpiWrapper::wtime( void )
{
#ifdef GEOS_USE_MPI
  return MPI_Wtime( );
#else
  return 0;
#endif

}

int MpiWrapper::activeWaitAny( const int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[], std::function< MPI_Request ( int ) > func )
{
  int cmp = 0;
  while( cmp < count )
  {
    int idx = 0;
    int err = waitAny( count, array_of_requests, &idx, array_of_statuses );
    if( err != MPI_SUCCESS )
      return err;
    if( idx != MPI_UNDEFINED )   // only if all(requests == MPI_REQUEST_NULL)
    {
      func( idx );
    }
    cmp++;
  }
  return MPI_SUCCESS;
}

int MpiWrapper::activeWaitSome( const int count,
                                MPI_Request array_of_requests[],
                                MPI_Status array_of_statuses[],
                                std::function< MPI_Request ( int ) > func )
{
  int cmp = 0;
  while( cmp < count )
  {
    int rcvd = 0;
    std::vector< int > indices( count, -1 );
    int err = waitSome( count, array_of_requests, &rcvd, &indices[0], array_of_statuses );
    if( err != MPI_SUCCESS )
      return err;
    if( rcvd > 0 )
    {
      for( int ii = 0; ii < rcvd; ++ii )
      {
        if( indices[ii] != MPI_UNDEFINED )
        {
          func( indices[ii] );
        }
      }
    }
    cmp += rcvd;
  }
  return MPI_SUCCESS;
}


int MpiWrapper::activeWaitSomeCompletePhase( const int participants,
                                             std::vector< std::tuple< MPI_Request *, MPI_Status *, std::function< MPI_Request ( int ) > > > const & phases )
{
  const int num_phases = phases.size();
  int err = 0;
  for( int phase = 0; phase < num_phases; ++phase )
  {
    MPI_Request * const requests = std::get< 0 >( phases[phase] );
    MPI_Status * const statuses = std::get< 1 >( phases[phase] );
    std::function< MPI_Request ( int ) > func = std::get< 2 >( phases[phase] );
    if( requests!=nullptr )
    {
      err = activeWaitSome( participants,
                            requests,
                            statuses,
                            func );
    }
    else
    {
      for( int idx = 0; idx < participants; ++idx )
      {
        func( idx );
      }
    }
    if( err != MPI_SUCCESS )
      break;
  }
  return err;
}

int MpiWrapper::activeWaitOrderedCompletePhase( const int participants,
                                                std::vector< std::tuple< MPI_Request *, MPI_Status *, std::function< MPI_Request ( int ) > > > const & phases )
{
  const int num_phases = phases.size();
  for( int phase = 0; phase < num_phases; ++phase )
  {
    MPI_Request * const requests = std::get< 0 >( phases[phase] );
    MPI_Status * const statuses = std::get< 1 >( phases[phase] );
    std::function< MPI_Request ( int ) > func = std::get< 2 >( phases[phase] );

    for( int idx = 0; idx < participants; ++idx )
    {
      if( requests!=nullptr )
      {
        wait( &requests[idx], &statuses[idx] );
      }
      func( idx );
    }
  }
  return MPI_SUCCESS;
}

int MpiWrapper::nodeCommSize()
{
  // if not initialized then we guess there is no MPI.
  if( !initialized() )
    return 1;

  int len;
  std::array< char, MPI_MAX_PROCESSOR_NAME + 1 > hostname;
  MPI_Get_processor_name( hostname.data(), &len );
  hostname[len] = '\0';
  int color = (int)std::hash< string >{} (hostname.data());
  if( color < 0 )
    color *= -1;

  /**
   * Create intra-node communicator
   */
  MPI_Comm nodeComm;
  int nodeCommSize;
  MPI_Comm_split( MPI_COMM_WORLD, color, -1, &nodeComm );
  MPI_Comm_size( nodeComm, &nodeCommSize );
  return nodeCommSize;
}
} /* namespace geos */

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
