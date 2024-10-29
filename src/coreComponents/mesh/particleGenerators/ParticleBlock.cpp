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
 * @file ParticleBlock.cpp
 *
 */

#include "ParticleBlock.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{
using namespace dataRepository;

ParticleBlock::ParticleBlock( string const & name, Group * const parent ):
  ParticleBlockABC( name, parent ),
  m_externalPropertyNames()
{}

ParticleBlock::~ParticleBlock()
{}

void ParticleBlock::setParticleType( ParticleType const particleType )
{
  m_particleType = particleType;

  switch( m_particleType )
  {
    case ParticleType::SinglePoint:
    {
      // Single Point
      m_hasRVectors = true; // For visualization only
      break;
    }
    case ParticleType::CPDI:
    {
      // CPDI
      m_hasRVectors = true;
      break;
    }
    case ParticleType::CPTI:
    {
      // CPTI
      m_hasRVectors = true;
      break;
    }
    case ParticleType::CPDI2:
    {
      // CPDI2
      m_hasRVectors = false;
      break;
    }
    default:
    {
      GEOS_ERROR( "Invalid particle type: " << m_particleType );
    }
  }
}

void ParticleBlock::resize( dataRepository::indexType const numParticles )
{
  Group::resize( numParticles );

  // Those members are not registered as wrappers because I do not want them
  // to be exposed though the `Group` public interface.
  m_localToGlobalMap.resize( numParticles );
}


}
