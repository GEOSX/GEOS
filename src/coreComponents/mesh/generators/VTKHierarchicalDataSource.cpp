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
 * @file VTKHierarchicalDataSource.cpp
 */

#include "mesh/generators/VTKHierarchicalDataSource.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#include <vtkXMLPartitionedDataSetCollectionReader.h>

namespace geos
{
using namespace dataRepository;

VTKHierarchicalDataSource::VTKHierarchicalDataSource( string const & name,
                                                      Group * const parent )
  : ExternalDataSourceBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setApplyDefaultValue( "attribute" ).
    setDescription( "Path to the mesh file" );
}

void VTKHierarchicalDataSource::open()
{
  string const extension = m_filePath.extension();
  GEOS_ERROR_IF( extension != "vtpc", "Unsupported vtk extension. File must be a vtpc file" );

  vtkNew< vtkXMLPartitionedDataSetCollectionReader > reader;
  reader->SetFileName( m_filePath.c_str());
  reader->Update();

  m_collection = vtkSmartPointer< vtkPartitionedDataSetCollection >( vtkPartitionedDataSetCollection::SafeDownCast( reader->GetOutput()));
  m_dataAssembly = vtkSmartPointer< vtkDataAssembly >( m_collection->GetDataAssembly() );

  GEOS_ERROR_IF( m_dataAssembly == nullptr, "No data Assembly attached to this collection" );
}

vtkSmartPointer< vtkPartitionedDataSet >
VTKHierarchicalDataSource::search( string const & path )
{
  int node = m_dataAssembly->GetFirstNodeByPath( path.c_str());
  GEOS_ERROR_IF( node == -1, "Node doesn't exist" );
  GEOS_ERROR_IF( m_dataAssembly->GetNumberOfChildren( node ) > 0, "Only leaf nodes can be queried." );

  std::vector< unsigned int > indices = m_dataAssembly->GetDataSetIndices( node, false );

  GEOS_ERROR_IF( indices.size() == 0, "Queried node has no dataset attached." );
  GEOS_ERROR_IF( indices.size() > 1, "Current constraint each tree node has only one dataset." );

  return m_collection->GetPartitionedDataSet( indices[0] );
}

REGISTER_CATALOG_ENTRY( ExternalDataSourceBase, VTKHierarchicalDataSource, string const &, Group * const )


}
