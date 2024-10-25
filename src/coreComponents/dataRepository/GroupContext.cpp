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
 * @file GroupContext.cpp
 */

#include "GroupContext.hpp"

namespace geos
{
namespace dataRepository
{


GroupContext::GroupContext( Group & group, string_view objectName ):
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
    path << ( std::prev( info.base() ) == lastFileInfo ? // Is `info` pointing to the last file info?
              GEOS_FMT( "/{}({},l.{})",
                        info->m_targetName, splitPath( info->m_filePath ).second, info->m_line ) :
              GEOS_FMT( "/{}", info->m_targetName ));
  }
  return path.str();
}

DataContext::ToStringInfo GroupContext::getToStringInfo() const
{ return ToStringInfo( m_targetName ); }


} /* namespace dataRepository */
} /* namespace geos */
