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
 * @file EpetraVector.cpp
 */

#include "EpetraVector.hpp"

#include "codingUtilities/RTTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraUtils.hpp"

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <EpetraExt_MultiVectorOut.h>

namespace geos
{

EpetraVector::EpetraVector()
  : VectorBase(),
  m_vec{}
{}

EpetraVector::EpetraVector( EpetraVector const & src )
  : EpetraVector()
{
  *this = src;
}

EpetraVector::EpetraVector( EpetraVector && src ) noexcept
  : EpetraVector()
{
  *this = std::move( src );
}

EpetraVector & EpetraVector::operator=( EpetraVector const & src )
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

EpetraVector & EpetraVector::operator=( EpetraVector && src ) noexcept
{
  if( &src != this )
  {
    m_vec = std::move( src.m_vec );
    src.m_vec = nullptr;
    VectorBase::operator=( std::move( src ) );
  }
  return *this;
}

EpetraVector::~EpetraVector() = default;

bool EpetraVector::created() const
{
  return bool( m_vec);
}

void EpetraVector::create( localIndex const localSize,
                           MPI_Comm const & MPI_PARAM( comm ) )
{
  VectorBase::create( localSize, comm );
  Epetra_Map const map( LvArray::integerConversion< long long >( -1 ),
                        LvArray::integerConversion< int >( localSize ),
                        0,
                        trilinos::EpetraComm( MPI_PARAM( comm ) ) );
  m_values.move( hostMemorySpace, false );
  m_vec = std::make_unique< Epetra_Vector >( View, map, m_values.data() );
}

void EpetraVector::set( real64 value )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( m_vec->PutScalar( value ) );
  touch();
}

void EpetraVector::rand( unsigned const seed )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( m_vec->SetSeed( seed ) );
  GEOS_LAI_CHECK_ERROR( m_vec->Random() );
  touch();
}

void EpetraVector::close()
{
  GEOS_LAI_ASSERT( !closed() );
  m_values.move( hostMemorySpace, false );
  m_closed = true;
}

void EpetraVector::touch()
{
  GEOS_LAI_ASSERT( ready() );
  m_values.registerTouch( hostMemorySpace );
}

void EpetraVector::reset()
{
  VectorBase::reset();
  m_vec.reset();
}

void EpetraVector::scale( real64 const scalingFactor )
{
  GEOS_LAI_ASSERT( ready() );

  if( !isEqual( scalingFactor, 1.0 ) )
  {
    GEOS_LAI_CHECK_ERROR( m_vec->Scale( scalingFactor ) );
    touch();
  }
}

void EpetraVector::reciprocal()
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_CHECK_ERROR( m_vec->Reciprocal( *m_vec ) );
  touch();
}

real64 EpetraVector::dot( EpetraVector const & vec ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( vec.ready() );
  GEOS_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  real64 tmp;
  GEOS_LAI_CHECK_ERROR( m_vec->Dot( vec.unwrapped(), &tmp ) );
  return tmp;
}

void EpetraVector::copy( EpetraVector const & x )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( x.ready() );
  GEOS_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOS_LAI_CHECK_ERROR( m_vec->Update( 1., x.unwrapped(), 0. ) );
  touch();
}

void EpetraVector::axpy( real64 const alpha,
                         EpetraVector const & x )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( x.ready() );
  GEOS_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOS_LAI_CHECK_ERROR( m_vec->Update( alpha, x.unwrapped(), 1. ) );
  touch();
}

void EpetraVector::axpby( real64 const alpha,
                          EpetraVector const & x,
                          real64 const beta )
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( x.ready() );
  GEOS_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOS_LAI_CHECK_ERROR( m_vec->Update( alpha, x.unwrapped(), beta ) );
  touch();
}

void EpetraVector::pointwiseProduct( EpetraVector const & x,
                                     EpetraVector & y ) const
{
  GEOS_LAI_ASSERT( ready() );
  GEOS_LAI_ASSERT( x.ready() );
  GEOS_LAI_ASSERT( y.ready() );
  GEOS_LAI_ASSERT_EQ( globalSize(), x.globalSize() );
  GEOS_LAI_ASSERT_EQ( globalSize(), y.globalSize() );

  GEOS_LAI_CHECK_ERROR( ( y.unwrapped() ).Multiply( 1.0, unwrapped(), x.unwrapped(), 0.0 ) );
  y.touch();
}

real64 EpetraVector::norm1() const
{
  GEOS_LAI_ASSERT( ready() );

  real64 tmp;
  m_vec->Norm1( &tmp );
  return tmp;
}

real64 EpetraVector::norm2() const
{
  GEOS_LAI_ASSERT( ready() );
  real64 tmp;
  m_vec->Norm2( &tmp );
  return tmp;
}

real64 EpetraVector::normInf() const
{
  GEOS_LAI_ASSERT( ready() );
  real64 tmp;
  m_vec->NormInf( &tmp );
  return tmp;
}

void EpetraVector::print( std::ostream & os ) const
{
  GEOS_LAI_ASSERT( ready() );
  m_vec->Print( os );
}

void EpetraVector::write( string const & filename,
                          LAIOutputFormat const format ) const
{
  GEOS_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::MATLAB_ASCII:
    {
      GEOS_LAI_CHECK_ERROR( EpetraExt::MultiVectorToMatlabFile( filename.c_str(), *m_vec ) );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      GEOS_LAI_CHECK_ERROR( EpetraExt::MultiVectorToMatrixMarketFile( filename.c_str(), *m_vec ) );
      break;
    }
    default:
      GEOS_ERROR( "Unsupported vector output format" );
  }
}

Epetra_Vector const & EpetraVector::unwrapped() const
{
  GEOS_LAI_ASSERT( created() );
  return *m_vec;
}

Epetra_Vector & EpetraVector::unwrapped()
{
  GEOS_LAI_ASSERT( created() );
  return *m_vec;
}

globalIndex EpetraVector::globalSize() const
{
  GEOS_LAI_ASSERT( created() );
  return m_vec->GlobalLength64();
}

localIndex EpetraVector::localSize() const
{
  GEOS_LAI_ASSERT( created() );
  return m_vec->MyLength();
}

globalIndex EpetraVector::ilower() const
{
  GEOS_LAI_ASSERT( created() );
  return m_vec->Map().MinMyGID64();
}

globalIndex EpetraVector::iupper() const
{
  GEOS_LAI_ASSERT( created() );
  return m_vec->Map().MaxMyGID64() + 1;
}

MPI_Comm EpetraVector::comm() const
{
  GEOS_LAI_ASSERT( created() );
#ifdef GEOS_USE_MPI
  return dynamicCast< Epetra_MpiComm const & >( m_vec->Map().Comm() ).Comm();
#else
  return MPI_COMM_GEOS;
#endif
}

} // end geos
