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
 * @file OutputBase.cpp
 */

#include "OutputBase.hpp"
#include "common/MpiWrapper.hpp"
#include "functions/FunctionBase.hpp"


namespace geos
{
using namespace dataRepository;

OutputBase::OutputBase( string const & name,
                        Group * const parent ):
  ExecutableGroup( name, parent ),
  m_childDirectory(),
  m_parallelThreads( 1 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeysStruct::childDirectoryString, &m_childDirectory ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Child directory path" );

  registerWrapper( viewKeysStruct::parallelThreadsString, &m_parallelThreads ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Number of plot files." );

}

OutputBase::~OutputBase()
{}

OutputBase::CatalogInterface::CatalogType & OutputBase::getCatalog()
{
  static OutputBase::CatalogInterface::CatalogType catalog;
  return catalog;
}



void OutputBase::initializePreSubGroups()
{
  // This command doesn't seem to work anymore
  // SetupDirectoryStructure();
}



string const & OutputBase::getOutputDirectory()
{
  static string m_outputDirectory;
  return m_outputDirectory;
}

void OutputBase::setOutputDirectory( string const & outputDir )
{
  string & outputDirectory = const_cast< string & >( getOutputDirectory() );
  outputDirectory = outputDir;
  FunctionBase::setOutputDirectory( outputDirectory );
}



string const & OutputBase::getFileNameRoot()
{
  static string m_fileNameRoot;
  return m_fileNameRoot;
}

void OutputBase::setFileNameRoot( string const & root )
{
  string & fileRootName = const_cast< string & >( getFileNameRoot() );
  fileRootName = root;
}


void OutputBase::setupDirectoryStructure()
{
  string childDirectory = m_childDirectory;

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  if( rank == 0 )
  {
    if( !childDirectory.empty())
    {
      makeDirsForPath( childDirectory );
    }
  }
}


} /* namespace geos */
