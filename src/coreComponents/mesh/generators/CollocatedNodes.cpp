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

#include "CollocatedNodes.hpp"

#include "common/MpiWrapper.hpp"


#include <vtkPointData.h>

namespace geos::vtk
{

CollocatedNodes::CollocatedNodes( string const & faceBlockName,
                                  vtkSmartPointer< vtkDataSet > faceMesh,
                                  bool isParallel )
{
  // The vtk field to the collocated nodes for fractures.
  string const COLLOCATED_NODES = "collocated_nodes";

  vtkIdTypeArray const * collocatedNodes = vtkIdTypeArray::FastDownCast( faceMesh->GetPointData()->GetArray( COLLOCATED_NODES.c_str() ) );
  if( isParallel )
  {
    // Depending on the parallel split, the vtk face mesh may be empty on a rank.
    // In that case, vtk will not provide any field for the emtpy mesh.
    // Therefore, not finding the duplicated nodes field on a rank cannot be interpreted as a globally missing field.
    // Converting the address into an integer and exchanging it through the MPI ranks let us find out
    // if the field is globally missing or not.
    std::uintptr_t const address = MpiWrapper::max( reinterpret_cast< std::uintptr_t >(collocatedNodes) );
    if( address == 0 )
    {
      GEOS_ERROR_IF( collocatedNodes == nullptr, "Could not find valid field \"" << COLLOCATED_NODES << "\" for fracture \"" << faceBlockName << "\"." );
    }
  }

  if( collocatedNodes )
  {
    init( collocatedNodes );
  }
}

void CollocatedNodes::init( vtkIdTypeArray const * collocatedNodes )
{
  vtkIdType const numTuples = collocatedNodes->GetNumberOfTuples();
  int const numComponents = collocatedNodes->GetNumberOfComponents();
  m_collocatedNodes.resize( numTuples );
  for( vtkIdType i = 0; i < numTuples; ++i )
  {
    m_collocatedNodes[i].reserve( numComponents );
  }

  for( vtkIdType i = 0; i < numTuples; ++i )
  {
    for( int j = 0; j < numComponents; ++j )
    {
      vtkIdType const tmp = collocatedNodes->GetTypedComponent( i, j );
      if( tmp > -1 )
      {
        m_collocatedNodes[i].emplace_back( tmp );
      }
    }
  }
}


} // geos::vtk
