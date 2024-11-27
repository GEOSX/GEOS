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


#include "ParticleSubRegion.hpp"

#include "common/TypeDispatch.hpp"
#include "mesh/MeshLevel.hpp"

namespace geos
{
using namespace dataRepository;
using namespace constitutive;

ParticleSubRegion::ParticleSubRegion( string const & name, Group * const parent ):
  ParticleSubRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::constitutiveGroupingString(), &m_constitutiveGrouping ).
    setSizedFromParent( 0 );
}

ParticleSubRegion::~ParticleSubRegion()
{
  // Left blank
}

void ParticleSubRegion::setParticleRank( int const rank )
{
  arrayView1d< int > const particleRank = m_particleRank;
  forAll< serialPolicy >( particleRank.size(), [=] GEOS_HOST_DEVICE ( localIndex const p )
  {
    particleRank[p] = rank;
  } );
}

void ParticleSubRegion::copyFromParticleBlock( ParticleBlockABC & particleBlock )
{
  // Defines the (unique) particle type of this cell particle region,
  // and its associated number of nodes, edges, faces.
  m_particleRank.resize( particleBlock.size() );
  m_particleType = particleBlock.getParticleType();
  m_particleID = particleBlock.getParticleID();
  m_particleGroup = particleBlock.getParticleGroup();
  m_particleSurfaceFlag = particleBlock.getParticleSurfaceFlag();
  m_particleDamage = particleBlock.getParticleDamage();
  m_particlePorosity = particleBlock.getParticlePorosity();
  m_particleTemperature = particleBlock.getParticleTemperature();
  m_particleStrengthScale = particleBlock.getParticleStrengthScale();
  m_particleCenter = particleBlock.getParticleCenter();
  m_particleVelocity = particleBlock.getParticleVelocity();
  m_particleMaterialDirection = particleBlock.getParticleMaterialDirection();
  m_particleVolume = particleBlock.getParticleVolume();
  m_particleRVectors = particleBlock.getParticleRVectors();
  m_hasRVectors = particleBlock.hasRVectors();
  m_particleSurfaceNormal = particleBlock.getParticleSurfaceNormal();
  m_particleSurfacePosition = particleBlock.getParticleSurfacePosition();
  m_particleSurfaceTraction = particleBlock.getParticleSurfaceTraction();

  // We call the `resize` member function of the cell to (nodes, edges, faces) relations,
  // before calling the `ParticleSubRegion::resize` in order to keep the first dimension.
  // Be careful when refactoring.
  this->resize( particleBlock.numParticles() );
  this->m_localToGlobalMap = particleBlock.localToGlobalMap();
  this->constructGlobalToLocalMap();
  particleBlock.forExternalProperties( [&]( WrapperBase const & wrapper )
  {
    types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
    {
      using ArrayType = camp::first< decltype( tupleOfTypes ) >;
      auto const src = Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();
      this->registerWrapper( wrapper.getName(), std::make_unique< ArrayType >( &src ) );
    }, wrapper );
  } );
}

// void ParticleSubRegion::copyFromParticleSubRegion( ParticleSubRegion & particleSubRegion )
// {
//   // Defines the (unique) particle type of this cell particle region,
//   // and its associated number of nodes, edges, faces.
//   m_particleRank = particleSubRegion.getParticleRank();
//   m_particleType = particleSubRegion.getParticleType();
//   m_particleID = particleSubRegion.getParticleID();
//   m_particleGroup = particleSubRegion.getParticleGroup();
//   m_particleSurfaceFlag = particleSubRegion.getParticleSurfaceFlag();
//   m_particleDamage = particleSubRegion.getParticleDamage();
//   m_particleStrengthScale = particleSubRegion.getParticleStrengthScale();
//   m_particleCenter = particleSubRegion.getParticleCenter();
//   m_particleVelocity = particleSubRegion.getParticleVelocity();
//   m_particleMaterialDirection = particleSubRegion.getParticleMaterialDirection();
//   m_particleVolume = particleSubRegion.getParticleVolume();
//   m_particleRVectors = particleSubRegion.getParticleRVectors();
//   m_hasRVectors = particleSubRegion.hasRVectors();
//   m_particleSurfaceNormal = particleSubRegion.getParticleSurfaceNormal();
//   m_particleSurfacePosition = particleSubRegion.getParticleSurfacePosition();
//   m_particleSurfaceTraction = particleSubRegion.getParticleSurfaceTraction();

//   // We call the `resize` member function of the cell to (nodes, edges, faces) relations,
//   // before calling the `ParticleSubRegion::resize` in order to keep the first dimension.
//   // Be careful when refactoring.
//   this->resize( particleSubRegion.size() );
//   this->m_localToGlobalMap = particleSubRegion.localToGlobalMap();
//   this->constructGlobalToLocalMap();

//   // CC: I believe the particleblock assosciated with this subregion is already registered from copyFromParticleBlock on startup so 
//   // there should not be a need to get the particle block of the source sub region during the copy (double check this!!!)

//   // Not sure what the following actually does, but ParticleSubRegion does not have a particleBlock or forExternalProperties
//   // particleBlock.forExternalProperties( [&]( WrapperBase & wrapper ) 
//   // {
//   //   types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), true, [&]( auto array )
//   //   {
//   //     using ArrayType = decltype( array );
//   //     Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
//   //     this->registerWrapper( wrapper.getName(), &wrapperT.reference() );
//   //   } );
//   // } );
// }

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleSubRegion, string const &, Group * const )

} /* namespace geos */
