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
 * @file ParticleRegion.cpp
 */

#include "ParticleRegion.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/SparsityPattern.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{
using namespace dataRepository;

ParticleRegion::ParticleRegion( string const & name, Group * const parent ):
  ParticleRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::sourceParticleBlockNamesString(), &m_particleBlockNames ).
    setInputFlag( InputFlags::OPTIONAL );
}

ParticleRegion::~ParticleRegion()
{}

void ParticleRegion::generateMesh( Group & particleBlocks )
{
  Group & particleSubRegions = this->getGroup( viewKeyStruct::particleSubRegions() );

  for( string const & particleBlockName : this->m_particleBlockNames )
  {
    ParticleSubRegion & subRegion = particleSubRegions.registerGroup< ParticleSubRegion >( particleBlockName );
    ParticleBlockABC & source = particleBlocks.getGroup< ParticleBlockABC >( subRegion.getName() );

    // Copy data from particle blocks into subregions. This also sets particle type.
    subRegion.copyFromParticleBlock( source );

    // Set the rank of particles on each subregion
    int const mpiRank = MpiWrapper::commRank( MPI_COMM_GEOS );
    subRegion.setParticleRank( mpiRank );

    // Set the number of vertices of the particles on each subregion
    subRegion.setNumberOfVerticesPerParticle( subRegion.getParticleType() );
  }
}

// TODO This should be changed to call a ParticleSubRegion::getParticleCoordinates (and/or getParticleCorners) on each subregion such that
// we can access those functions directly if needed
array2d< real64 > ParticleRegion::getParticleCorners() const
{
  int const size = 8 * ( this->size() ); // number of particle corners in this region for a parallelpiped CPDI particle type (TODO add
                                         // support for general types)
  array2d< real64 > coords( size, 3 );
  int index = 0;
  int const signs[8][3] = { { 1, 1, 1},   // These have to be in the order VTK wants
    { 1, 1, -1},
    { 1, -1, 1},
    { 1, -1, -1},
    {-1, 1, 1},
    {-1, 1, -1},
    {-1, -1, 1},
    {-1, -1, -1} };
  this->forParticleSubRegions( [&]( auto & subRegion )
  {
    arrayView2d< real64 const > const particleCenter = subRegion.getParticleCenter();
    if( !subRegion.hasRVectors() ) // if no r-vectors, return coordinates of a cube centered at the particle with side length = volume^(1/3)
    {
      arrayView1d< real64 const > const particleVolume = subRegion.getParticleVolume();
      forAll< serialPolicy >( subRegion.size(), [=, &coords, &index] GEOS_HOST ( localIndex const p )
        {
          real64 a = 0.5 * pow( particleVolume[p], 1.0/3.0 ); // cube half-side-length
          for( int corner=0; corner<8; corner++ )
          {
            for( int i=0; i<3; i++ )
            {
              coords[index][i] = particleCenter[p][i] + signs[corner][i] * a;
            }
            index++;
          }
        } );
    }
    else
    {
      arrayView3d< real64 const > const particleRVectors = subRegion.getParticleRVectors();
      forAll< serialPolicy >( subRegion.size(), [=, &coords, &index] GEOS_HOST ( localIndex const p )
        {
          for( int corner=0; corner<8; corner++ )
          {
            for( int i=0; i<3; i++ )
            {
              coords[index][i] = particleCenter[p][i] + signs[corner][0] * particleRVectors[p][0][i] + signs[corner][1] * particleRVectors[p][1][i] + signs[corner][2] * particleRVectors[p][2][i];
            }
            index++;
          }
        } );
    }
  } );
  return coords;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleRegion, string const &, Group * const )

} /* namespace geos */
