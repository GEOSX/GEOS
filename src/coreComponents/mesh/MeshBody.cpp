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
 * @file MeshBody.cpp
 */

#include "MeshBody.hpp"
#include "MeshLevel.hpp"

namespace geos
{

using namespace dataRepository;

MeshBody::MeshBody( string const & name,
                    Group * const parent ):
  Group( name, parent ),
  m_meshLevels( registerGroup( groupStructKeys::meshLevelsString() ) ),
  m_globalLengthScale( 0 ),
  m_hasParticles( false )
{}

MeshLevel & MeshBody::createMeshLevel( localIndex const newLevel )
{
  return m_meshLevels.registerGroup< MeshLevel >( intToMeshLevelString( newLevel ) );
}

MeshLevel & MeshBody::createMeshLevel( string const & name )
{
  return m_meshLevels.registerGroup< MeshLevel >( name );
}

MeshLevel & MeshBody::createMeshLevel( string const & sourceLevelName,
                                       string const & newLevelName,
                                       int const order )
{
  MeshLevel const & sourceMeshLevel = this->getMeshLevel( sourceLevelName );
  return m_meshLevels.registerGroup( newLevelName,
                                     std::make_unique< MeshLevel >( newLevelName,
                                                                    this,
                                                                    sourceMeshLevel,
                                                                    order ) );
}

MeshLevel & MeshBody::createShallowMeshLevel( string const & sourceLevelName,
                                              string const & newLevelName )
{
  MeshLevel & sourceMeshLevel = this->getMeshLevel( sourceLevelName );

  MeshLevel & rval = m_meshLevels.registerGroup( newLevelName,
                                                 std::make_unique< MeshLevel >( newLevelName,
                                                                                this,
                                                                                sourceMeshLevel ) );
  rval.setRestartFlags( RestartFlags::NO_WRITE );

  return rval;
}


void MeshBody::setGlobalLengthScale( real64 scale )
{
  m_globalLengthScale = scale;
}

string MeshBody::intToMeshLevelString( localIndex const meshLevel )
{
  return GEOS_FMT( "Level{}", meshLevel );
}

void MeshBody::setHasParticles( bool hasParticles )
{
  m_hasParticles = hasParticles;
}


} /* namespace geos */
