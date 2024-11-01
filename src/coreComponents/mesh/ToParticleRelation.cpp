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
 * @file ToParticleRelation.cpp
 */

#include "ToParticleRelation.hpp"

#include "ParticleManager.hpp"

namespace geos
{


void erase( OrderedVariableToManyParticleRelation & relation,
            localIndex const firstIndex,
            localIndex const er,
            localIndex const esr,
            localIndex const ei )
{
  for( localIndex a=relation.m_toParticleRegion.sizeOfArray( firstIndex )-1; a>=0; --a )
  {
    if( er==relation.m_toParticleRegion[firstIndex][a] &&
        esr==relation.m_toParticleSubRegion[firstIndex][a] &&
        ei==relation.m_toParticleIndex[firstIndex][a] )
    {
      relation.m_numParticles[firstIndex]--;
      relation.m_toParticleRegion.eraseFromArray( firstIndex, a );
      relation.m_toParticleSubRegion.eraseFromArray( firstIndex, a );
      relation.m_toParticleIndex.eraseFromArray( firstIndex, a );
    }
  }
}

void insert( OrderedVariableToManyParticleRelation & relation,
             localIndex const firstIndex,
             localIndex const er,
             localIndex const esr,
             localIndex const ei )
{
  bool alreadyPresent = false;
  for( localIndex a=0; a<relation.m_toParticleRegion.sizeOfArray( firstIndex ); ++a )
  {
    if( er==relation.m_toParticleRegion[firstIndex][a] &&
        esr==relation.m_toParticleSubRegion[firstIndex][a] &&
        ei==relation.m_toParticleIndex[firstIndex][a] )
    {
      alreadyPresent = true;
    }
  }
  if( !alreadyPresent )
  {
    relation.m_numParticles[firstIndex]++;
    relation.m_toParticleRegion.emplaceBack( firstIndex, er );
    relation.m_toParticleSubRegion.emplaceBack( firstIndex, esr );
    relation.m_toParticleIndex.emplaceBack( firstIndex, ei );
  }
}

void fastInsert( OrderedVariableToManyParticleRelation & relation,
                 localIndex const firstIndex,
                 localIndex const er,
                 localIndex const esr,
                 localIndex const ei )
{
  relation.m_numParticles[firstIndex]++;
  relation.m_toParticleRegion.emplaceBack( firstIndex, er );
  relation.m_toParticleSubRegion.emplaceBack( firstIndex, esr );
  relation.m_toParticleIndex.emplaceBack( firstIndex, ei );
}

void insertMany( OrderedVariableToManyParticleRelation & relation,
                 localIndex const firstIndex,
                 std::vector< localIndex > const & erArray,
                 std::vector< localIndex > const & esrArray,
                 std::vector< localIndex > const & eiArray )
{
  relation.m_numParticles[firstIndex] += erArray.size();
  relation.m_toParticleRegion.appendToArray( firstIndex, erArray.begin(), erArray.end() );
  relation.m_toParticleSubRegion.appendToArray( firstIndex, esrArray.begin(), esrArray.end() );
  relation.m_toParticleIndex.appendToArray( firstIndex, eiArray.begin(), eiArray.end() );
}

void reserveNeighbors( OrderedVariableToManyParticleRelation & relation,
                       int const numToReserve )
{
  relation.m_toParticleRegion.reserveValues( numToReserve );
  relation.m_toParticleSubRegion.reserveValues( numToReserve );
  relation.m_toParticleIndex.reserveValues( numToReserve );
}

} /* namespace geos */
