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
 * @file HypreVector.cpp
 */

#include "HypreVector.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"

#include <_hypre_IJ_mv.h>

#include <iomanip>

namespace geos
{

HypreVector::HypreVector()
  : VectorBase(),
  m_vec{}
{}

// Copy constructor
HypreVector::HypreVector( HypreVector const & src )
  : HypreVector()
{
  *this = src;
}

// Move constructor
HypreVector::HypreVector( HypreVector && src ) noexcept
  : HypreVector()
{
  *this = std::move( src );
}

HypreVector & HypreVector::operator=( HypreVector const & src )
{
  if( &src != this )
  {
    reset();
    if( src.created() )
    {
      create( src.localSize(), src.comm() );
      copy( src );
    }
  }
  return *this;
}

HypreVector & HypreVector::operator=( HypreVector && src ) noexcept
{
  if( &src != this )
  {
    m_vec = src.m_vec;
    src.m_vec = nullptr;
    VectorBase::operator=( std::move( src ) );
  }
  return *this;
}

void HypreVector::reset()
{
  VectorBase::reset();
  if( m_vec )
  {
    hypre_ParVectorDestroy( m_vec );
    m_vec = nullptr;
  }
}

HypreVector::~HypreVector()
{
  reset();
}

void HypreVector::create( localIndex const localSize,
                          MPI_Comm const & comm )
{
  VectorBase::create( localSize, comm );

  // Compute partitioning information
  HYPRE_BigInt partitioning[2];
  partitioning[0] = MpiWrapper::prefixSum< HYPRE_BigInt >( localSize, comm );
  partitioning[1] = partitioning[0] + localSize;
  HYPRE_BigInt globalSize = partitioning[1];
  MpiWrapper::broadcast( globalSize, MpiWrapper::commSize( comm ) - 1, comm );

  // Set up the parallel and local vector data structures
  m_vec = hypre_ParVectorCreate( comm, globalSize, partitioning );
  hypre_ParVectorOwnsData( m_vec ) = false;

  hypre_Vector * const localVector = hypre_ParVectorLocalVector( m_vec );
  hypre_VectorOwnsData( localVector ) = false;

  // Inject the memory managed by m_values in the correct space into hypre vector
  m_values.move( hypre::memorySpace, false );
  hypre_VectorData( localVector ) = m_values.data();

  // Complete the initialization (vector will not allocate if data is already set)
  GEOS_LAI_CHECK_ERROR( hypre_ParVectorInitialize_v2( m_vec, hypre::memoryLocation ) );
  GEOS_LAI_CHECK_ERROR( hypre_ParVectorSetConstantValues( m_vec, 0.0 ) );
}

bool HypreVector::created() const
{
  return m_vec != nullptr;
}

void HypreVector::set( real64 value )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParVectorSetConstantValues( m_vec, value ) );
  touch();
}

void HypreVector::rand( unsigned const seed )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( HYPRE_ParVectorSetRandomValues( m_vec, seed ) );
  touch();
}

void HypreVector::close()
{
  GEOS_LAI_ASSERT( !closed() );
  m_values.move( hypre::memorySpace, false );
  m_closed = true;
}

void HypreVector::touch()
{
  GEOS_LAI_ASSERT( ready() );
  m_values.registerTouch( hypre::memorySpace );
}

void HypreVector::scale( real64 const scalingFactor )
{
  GEOS_LAI_ASSERT( ready() );
  if( !isEqual( scalingFactor, 1.0 ) )
  {
    GEOS_LAI_CHECK_ERROR( HYPRE_ParVectorScale( scalingFactor, m_vec ) );
    touch();
  }
}

void HypreVector::reciprocal()
{
  GEOS_LAI_ASSERT( ready() );
  arrayView1d< real64 > values = m_values.toView();
  forAll< hypre::execPolicy >( localSize(), [values] GEOS_HYPRE_DEVICE ( localIndex const i )
  {
    values[i] = 1.0 / values[i];
  } );
}

real64 HypreVector::dot( HypreVector const & vec ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( vec.ready() );
  GEOS_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  HYPRE_Real result;
  GEOS_LAI_CHECK_ERROR( HYPRE_ParVectorInnerProd( m_vec, vec.m_vec, &result ) );
  return result;
}

void HypreVector::copy( HypreVector const & x )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( x.ready() );
  GEOS_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOS_LAI_CHECK_ERROR( HYPRE_ParVectorCopy( x.m_vec, m_vec ) );
  touch();
}

void HypreVector::axpy( real64 const alpha,
                        HypreVector const & x )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( x.ready() );
  GEOS_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  if( !isEqual( alpha, 0.0 ) )
  {
    if( &x != this )
    {
      GEOS_LAI_CHECK_ERROR( HYPRE_ParVectorAxpy( alpha, x.m_vec, m_vec ) );
      touch();
    }
    else
    {
      scale( 1.0 + alpha );
    }
  }
}

void HypreVector::axpby( real64 const alpha,
                         HypreVector const & x,
                         real64 const beta )
{
  if( &x != this )
  {
    scale( beta );
    axpy( alpha, x );
  }
  else
  {
    scale( alpha + beta );
  }
}

void HypreVector::pointwiseProduct( HypreVector const & x,
                                    HypreVector & y ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( x.ready() );
  GEOS_LAI_ASSERT( y.ready() );
  GEOS_LAI_ASSERT_EQ( localSize(), x.localSize() );
  GEOS_LAI_ASSERT_EQ( localSize(), y.localSize() );

  arrayView1d< real64 const > const my_values = m_values.toViewConst();
  arrayView1d< real64 const > const x_values = x.m_values.toViewConst();
  arrayView1d< real64 > const y_values = y.m_values.toView();
  forAll< hypre::execPolicy >( localSize(), [y_values, my_values, x_values] GEOS_HYPRE_DEVICE ( localIndex const i )
  {
    y_values[i] = my_values[i] * x_values[i];
  } );
}

real64 HypreVector::norm1() const
{
  GEOS_LAI_ASSERT( ready() );

  arrayView1d< real64 const > values = m_values.toViewConst();
  RAJA::ReduceSum< ReducePolicy< hypre::execPolicy >, real64 > localNorm( 0.0 );
  forAll< hypre::execPolicy >( localSize(), [localNorm, values] GEOS_HYPRE_DEVICE ( localIndex const i )
  {
    localNorm += LvArray::math::abs( values[i] );
  } );
  return MpiWrapper::sum( localNorm.get(), comm() );
}

real64 HypreVector::norm2() const
{
  GEOS_LAI_ASSERT( ready() );
  return std::sqrt( dot( *this ) );
}

real64 HypreVector::normInf() const
{
  GEOS_LAI_ASSERT( ready() );

  arrayView1d< real64 const > values = m_values.toViewConst();
  RAJA::ReduceMax< ReducePolicy< hypre::execPolicy >, real64 > localNorm( 0.0 );
  forAll< hypre::execPolicy >( localSize(), [localNorm, values] GEOS_HYPRE_DEVICE ( localIndex const i )
  {
    localNorm.max( LvArray::math::abs( values[i] ) );
  } );
  return MpiWrapper::max( localNorm.get(), comm() );
}

globalIndex HypreVector::globalSize() const
{
  GEOS_LAI_ASSERT( created() );
  return hypre_ParVectorGlobalSize( m_vec );
}

localIndex HypreVector::localSize() const
{
  GEOS_LAI_ASSERT( created() );
  return hypre_ParVectorActualLocalSize( m_vec );
}

globalIndex HypreVector::ilower() const
{
  GEOS_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( hypre_ParVectorFirstIndex( m_vec ) );
}

globalIndex HypreVector::iupper() const
{
  GEOS_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >( hypre_ParVectorLastIndex( m_vec ) ) + 1;
}

void HypreVector::print( std::ostream & os ) const
{
  GEOS_LAI_ASSERT( ready() );

  int const myRank = MpiWrapper::commRank( comm() );
  int const numProcs = MpiWrapper::commSize( comm() );
  char str[77];

  constexpr char const lineFormat[] = "{:>11}{:>18}{:>28.16e}\n";
  constexpr char const headFormat[] = "{:>11}{:>18}{:>28}\n";

  if( myRank == 0 )
  {
    GEOS_FMT_TO( str, sizeof( str ), headFormat, "MPI_Process", "GlobalRowID", "Value" );
    os << str;
  }

  for( int rank = 0; rank < numProcs; ++rank )
  {
    MpiWrapper::barrier( comm() );
    if( rank == myRank )
    {
      arrayView1d< real64 const > const data = values();
      globalIndex const firstRowID = ilower();
      forAll< serialPolicy >( localSize(), [&, data]( localIndex const i )
      {
        GEOS_FMT_TO( str, sizeof( str ), lineFormat,
                     rank,
                     firstRowID + i,
                     data[i] );
        os << str;
      } );
    }
  }
}

void HypreVector::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOS_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      GEOS_LAI_CHECK_ERROR( hypre_ParVectorPrint( m_vec, filename.c_str() ) );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      int const rank = MpiWrapper::commRank( comm() );

      // Write MatrixMarket header
      if( rank == 0 )
      {
        std::ofstream os( filename );
        GEOS_ERROR_IF( !os, GEOS_FMT( "Unable to open file for writing: {}", filename ) );
        os << "%%MatrixMarket matrix array real general\n";
        os << GEOS_FMT( "{} {}\n", globalSize(), 1 );
      }

      if( globalSize() > 0 )
      {
        // Copy distributed parVector in a local vector on every process with at least one component
        // Warning: works for a parVector that is smaller than 2^31-1
        hypre_Vector * const fullVector = (hypre_Vector *)hypre::parVectorToVectorAll( m_vec );

        // Identify the smallest process where vector exists
        int const printRank = MpiWrapper::min( fullVector ? rank : MpiWrapper::commSize( comm() ), comm() );

        // Write to file vector
        if( MpiWrapper::commRank( comm() ) == printRank )
        {
          std::ofstream os( filename, std::ios_base::app );
          GEOS_ERROR_IF( !os, GEOS_FMT( "Unable to open file for writing on rank {}: {}", rank, filename ) );
          char str[32];

          HYPRE_Real const * const data = hypre_VectorData( fullVector );
          HYPRE_Int const size = hypre_VectorSize( fullVector );

          for( HYPRE_Int i = 0; i < size; i++ )
          {
            GEOS_FMT_TO( str, sizeof( str ), "{:>28.16e}\n", data[i] );
            os << str;
          }
        }

        // Destroy temporary vector
        GEOS_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( fullVector ) );
      }
      break;
    }
    default:
    {
      GEOS_ERROR( "Unsupported vector output format" );
    }
  }
}

HYPRE_ParVector const & HypreVector::unwrapped() const
{
  return m_vec;
}

MPI_Comm HypreVector::comm() const
{
  GEOS_LAI_ASSERT( created() );
  return hypre_ParVectorComm( m_vec );
}

} // end namespace geos
