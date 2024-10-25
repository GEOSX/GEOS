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

#include "VTKPVDWriter.hpp"

#include "common/MpiWrapper.hpp"

namespace geos
{
namespace vtk
{
VTKPVDWriter::VTKPVDWriter( string fileName ):
  m_fileName( std::move( fileName ) )
{
  // Declaration of XML version
  auto declarationNode = m_pvdFile.appendChild( pugi::node_declaration );
  declarationNode.append_attribute( "version" ) = "1.0";

  // Declaration of the node VTKFile
  auto vtkFileNode = m_pvdFile.appendChild( "VTKFile" );
  vtkFileNode.append_attribute( "type" ) = "Collection";
  vtkFileNode.append_attribute( "version" ) = "0.1";

  vtkFileNode.append_child( "Collection" );
}

void VTKPVDWriter::setFileName( string fileName )
{
  m_fileName = std::move( fileName );
}

void VTKPVDWriter::save() const
{
  m_pvdFile.saveFile( m_fileName );
}

void VTKPVDWriter::addData( real64 time, string const & filePath ) const
{
  auto collectionNode = m_pvdFile.getChild( "VTKFile" ).child( "Collection" );
  auto dataSetNode = collectionNode.append_child( "DataSet" );
  dataSetNode.append_attribute( "timestep" ) = time;
  dataSetNode.append_attribute( "file" ) = filePath.c_str();
}

void VTKPVDWriter::reinitData() const
{
  auto collectionNode = m_pvdFile.getChild( "VTKFile" ).child( "Collection" );

  while( collectionNode.remove_child( "DataSet" ) )
  {}
}
}
}
