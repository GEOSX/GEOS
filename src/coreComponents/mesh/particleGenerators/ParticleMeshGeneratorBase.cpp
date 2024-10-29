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

#include "ParticleMeshGeneratorBase.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "mesh/particleGenerators/ParticleBlockManager.hpp"
#include "common/GeosxMacros.hpp"

namespace geos
{
using namespace dataRepository;

ParticleMeshGeneratorBase::ParticleMeshGeneratorBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

Group * ParticleMeshGeneratorBase::createChild( string const & GEOS_UNUSED_PARAM( childKey ),
                                                string const & GEOS_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void ParticleMeshGeneratorBase::expandObjectCatalogs()
{
  // During schema generation, register one of each type derived from WellGeneratorBase here
  for( auto & catalogIter: WellGeneratorBase::getCatalog())
  {
    createChild( catalogIter.first, catalogIter.first );
  }
}

ParticleMeshGeneratorBase::CatalogInterface::CatalogType & ParticleMeshGeneratorBase::getCatalog()
{
  static ParticleMeshGeneratorBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void ParticleMeshGeneratorBase::generateMesh( Group & parent, SpatialPartition & partition )
{
  MeshBody & meshBody = dynamic_cast< MeshBody & >( parent );
  if( meshBody.hasParticles() )
  {
    ParticleBlockManager & particleBlockManager = parent.registerGroup< ParticleBlockManager >( keys::particleManager );

    MeshLevel & meshLevel0 = meshBody.getBaseDiscretization();
    ParticleManager & particleManager = meshLevel0.getParticleManager();

    fillParticleBlockManager( particleBlockManager, particleManager, partition );
  }
  else
  {
    GEOS_ERROR( "Internal error. MeshBody " << meshBody.getName() << " should have particles." );
  }
}

}
