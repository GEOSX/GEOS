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
 * @file EpetraExport.cpp
 */

#include "EpetraExport.hpp"

#include "common/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraMatrix.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraUtils.hpp"

#include <Epetra_Map.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>

namespace geos
{

EpetraExport::EpetraExport() = default;

EpetraExport::EpetraExport( EpetraMatrix const & mat,
                            integer targetRank )
  : m_targetRank( targetRank )
{
  globalIndex const numGlobalRows = mat.numGlobalRows();
  localIndex const numLocalRows = ( m_targetRank == MpiWrapper::commRank( mat.comm() ) ) ? numGlobalRows : 0;
  m_serialMap = std::make_unique< Epetra_Map >( numGlobalRows, numLocalRows, 0, mat.unwrapped().Comm() );
  m_serialImport = std::make_unique< Epetra_Import >( *m_serialMap, mat.unwrapped().RowMap() );
}

EpetraExport::~EpetraExport() = default;

template< typename OFFSET_TYPE, typename COLUMN_TYPE >
void EpetraExport::exportCRS( EpetraMatrix const & mat,
                              arrayView1d< OFFSET_TYPE > const & rowOffsets,
                              arrayView1d< COLUMN_TYPE > const & colIndices,
                              arrayView1d< real64 > const & values ) const
{

  int const rank = MpiWrapper::commRank( mat.comm() );
  Epetra_CrsMatrix const * localMatrix = &mat.unwrapped();

  rowOffsets.move( hostMemorySpace, false );
  colIndices.move( hostMemorySpace, false );
  values.move( hostMemorySpace, false );

  // import on target rank if needed
  std::unique_ptr< Epetra_CrsMatrix > serialMatrix;
  if( m_targetRank >= 0 )
  {
    serialMatrix = std::make_unique< Epetra_CrsMatrix >( mat.unwrapped(), *m_serialImport, m_serialMap.get() );
    localMatrix = serialMatrix.get();
  }

  // export the raw CRS data
  if( m_targetRank < 0 || m_targetRank == rank )
  {
    int * ia;
    int * ja;
    real64 * va;
    GEOS_LAI_CHECK_ERROR( localMatrix->ExtractCrsDataPointers( ia, ja, va ) );

    // contains the global ID of local columns
    globalIndex const * const globalColumns = localMatrix->ColMap().MyGlobalElements64();
    std::transform( ia, ia + localMatrix->NumMyRows() + 1, rowOffsets.data(), LvArray::integerConversion< OFFSET_TYPE, int > );
    std::transform( ja, ja + localMatrix->NumMyNonzeros(), colIndices.data(),
                    [globalColumns]( int const i ){ return LvArray::integerConversion< COLUMN_TYPE >( globalColumns[i] ); } );
    std::copy( va, va + localMatrix->NumMyNonzeros(), values.data() );
  }
}

void EpetraExport::exportVector( EpetraVector const & vec,
                                 arrayView1d< real64 > const & values ) const
{
  values.move( hostMemorySpace, false );
  if( m_targetRank >= 0 )
  {
    // Create a local vector that directly wraps the user-provided buffer and gather
    Epetra_MultiVector localVec( View, *m_serialMap, values.data(), 0, 1 );
    localVec.Import( vec.unwrapped(), *m_serialImport, Insert );
  }
  else
  {
    arrayView1d< real64 const > const data = vec.values();
    data.move( hostMemorySpace, false );
    std::copy( data.begin(), data.end(), values.data() );
  }
}

void EpetraExport::importVector( arrayView1d< const real64 > const & values,
                                 EpetraVector & vec ) const
{
  values.move( hostMemorySpace, false );
  if( m_targetRank >= 0 )
  {
    // HACK: const_cast required in order to create an Epetra vector that wraps user data;
    //       we promise the values won't be changed since the only thing we do is scatter.
    Epetra_MultiVector localVector( View, *m_serialMap, const_cast< real64 * >( values.data() ), LvArray::integerConversion< int >( vec.globalSize() ), 1 );
    vec.unwrapped().Export( localVector, *m_serialImport, Insert );
  }
  else
  {
    arrayView1d< real64 > const data = vec.open();
    std::copy( values.data(), values.data() + vec.localSize(), data.begin() );
    vec.close();
  }
}

/**
 * Explicit template instantiation macro for EpetraExport::exportCRS.
 * We need to explicitly instantiate this function template because:
 * - we want CRS consumers to specify their own destination buffer types;
 * - we're "hiding" Epetra headers from being included by consumer code.
 */
#define INST_EPETRA_EXPORT_CRS( OFFSET_TYPE, COLUMN_TYPE ) \
  template void \
  EpetraExport::exportCRS< OFFSET_TYPE, COLUMN_TYPE >( EpetraMatrix const &, \
                                                       arrayView1d< OFFSET_TYPE > const &, \
                                                       arrayView1d< COLUMN_TYPE > const &, \
                                                       arrayView1d< real64 > const & ) const

// Add other instantiations as needed (only use built-in types)
INST_EPETRA_EXPORT_CRS( int, int );
INST_EPETRA_EXPORT_CRS( long, long );
INST_EPETRA_EXPORT_CRS( long long, long long );

#undef INST_EPETRA_EXPORT_CRS

} // namespace geos
