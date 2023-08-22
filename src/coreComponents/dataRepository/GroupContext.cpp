/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file GroupContext.cpp
 */

#include "GroupContext.hpp"

namespace geos
{
namespace dataRepository
{


GroupContext::GroupContext( Group & group, string const & objectName ):
  DataContext( objectName ),
  m_group( group )
{}
GroupContext::GroupContext( Group & group ):
  GroupContext( group, group.getName() )
{}

string GroupContext::toString() const
{
  std::vector< ToStringInfo > parentsInfo;
  for( Group const * group = &m_group; group->hasParent(); group = &group->getParent() )
  {
    parentsInfo.push_back( group->getDataContext().getToStringInfo() );
  }

  std::ostringstream path;
  auto lastFileInfo = std::find_if( parentsInfo.begin(), parentsInfo.end(),
                                    []( ToStringInfo const & i ) { return i.hasInputFileInfo(); } );
  for( auto info = parentsInfo.rbegin(); info != parentsInfo.rend(); ++info )
  {
    path << ( info.base() != lastFileInfo ?
              GEOS_FMT( "/{}", info->m_targetName ) :
              GEOS_FMT( "/{}({},l.{})", info->m_targetName, info->m_filePath, info->m_line ) );
  }
  return path.str();
}

DataContext::ToStringInfo GroupContext::getToStringInfo() const
{ return ToStringInfo( m_targetName ); }


} /* namespace dataRepository */
} /* namespace geos */
