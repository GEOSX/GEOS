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
 * @file DataContext.cpp
 */

#include "DataContext.hpp"

namespace geos
{
namespace dataRepository
{


DataContext::DataContext( string const & targetName ):
  m_targetName( targetName )
{}

std::ostream & operator<<( std::ostream & os, DataContext const & ctx )
{
  os << ctx.toString();
  return os;
}

DataContext::ToStringInfo::ToStringInfo( string const & targetName, string const & filePath, size_t line ):
  m_targetName( targetName ),
  m_filePath( filePath ),
  m_line( line )
{}
DataContext::ToStringInfo::ToStringInfo( string const & targetName ):
  m_targetName( targetName )
{}


/**
 * @return the node 'name' attribute if it exists, return the node tag name otherwise.
 * @param node the target node.
 */
string getNodeName( xmlWrapper::xmlNode const & node )
{
  xmlWrapper::xmlAttribute const nameAtt = node.attribute( "name" );
  if( !nameAtt.empty() )
  {
    return string( node.attribute( "name" ).value() );
  }
  else
  {
    return string( node.name() );
  }
}

DataFileContext::DataFileContext( xmlWrapper::xmlNode const & targetNode,
                                  xmlWrapper::xmlNodePos const & nodePos ):
  DataContext( getNodeName( targetNode ) ),
  m_typeName( targetNode.name() ),
  m_filePath( nodePos.filePath ),
  m_line( nodePos.line ),
  m_offsetInLine( nodePos.offsetInLine ),
  m_offset( nodePos.offset )
{}

DataFileContext::DataFileContext( xmlWrapper::xmlNode const & targetNode,
                                  xmlWrapper::xmlAttribute const & att,
                                  xmlWrapper::xmlAttributePos const & attPos ):
  DataContext( getNodeName( targetNode ) + '/' + att.name() ),
  m_typeName( att.name() ),
  m_filePath( attPos.filePath ),
  m_line( attPos.line ),
  m_offsetInLine( attPos.offsetInLine ),
  m_offset( attPos.offset )
{}

string DataFileContext::toString() const
{
  if( m_line != xmlWrapper::xmlDocument::npos )
  {
    return GEOS_FMT( "{} ({}, l.{})", m_targetName, splitPath( m_filePath ).second, m_line );
  }
  else if( m_offset != xmlWrapper::xmlDocument::npos )
  {
    return GEOS_FMT( "{} ({}, offset {})", m_targetName, splitPath( m_filePath ).second, m_offset );
  }
  else
  {
    return GEOS_FMT( "{} (Source file not found)", m_targetName );
  }
}

DataContext::ToStringInfo DataFileContext::getToStringInfo() const
{ return ToStringInfo( m_targetName, m_filePath, m_line ); }



} /* namespace dataRepository */
} /* namespace geos */
