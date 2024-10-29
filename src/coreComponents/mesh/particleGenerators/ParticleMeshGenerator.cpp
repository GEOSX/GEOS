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
 * @file ParticleMeshGenerator.cpp
 */

#include "ParticleMeshGenerator.hpp"

#include "ParticleBlockManager.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

#include <cmath>

namespace geos
{
using namespace dataRepository;

ParticleMeshGenerator::ParticleMeshGenerator( string const & name, Group * const parent ):
  ParticleMeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::particleFilePathString(), &m_particleFilePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the particle file" );

  // registerWrapper( viewKeyStruct::headerFilePathString(), &m_headerFilePath ).
  //   setInputFlag( InputFlags::REQUIRED ).
  //   setRestartFlags( RestartFlags::NO_WRITE ).
  //   setDescription( "path to the header file" );

  registerWrapper( viewKeyStruct::particleBlockNamesString(), &m_blockNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of each particle block" );

  registerWrapper( viewKeyStruct::particleTypesString(), &m_particleTypes ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Particle types of each particle block" );
}

Group * ParticleMeshGenerator::createChild( string const & GEOS_UNUSED_PARAM( childKey ),
                                            string const & GEOS_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void ParticleMeshGenerator::postInputInitialization()
{
  ParticleMeshGeneratorBase::postInputInitialization();

  GEOS_ERROR_IF(m_blockNames.size() == 0, "No particle blocks were specified! Must specify at least one particle block.");
  GEOS_ERROR_IF(m_blockNames.size() != m_particleTypes.size(), "The particle block and type lists must have the same size.");

  for( int i = 0; i < m_particleTypes.size(); i++)
  {
    for( int j = 0; j < EnumSize<ParticleType>; j++ )
    {
      if( m_particleTypes[i] == EnumStrings< ParticleType >::toString( static_cast< ParticleType >( j ) ) )
      {
        break;
      }

      if( j == EnumSize<ParticleType>-1)
      {
        GEOS_ERROR( "No particle type of " << m_particleTypes[i] << " in particleTypes input. Available options are SinglePoint (0), SinglePointBSpline (1), CPDI (2), CPTI (3), CPDI2 (4)."); // TODO Make available options dynamic
      }
    }
  }
}

void ParticleMeshGenerator::fillParticleBlockManager( ParticleBlockManager & particleBlockManager, ParticleManager & particleManager, SpatialPartition const & partition )
{
  GEOS_MARK_FUNCTION;

  int numBlocks = m_blockNames.size();

  // This should probably be handled elsewhere:
  std::vector<int> blockMaterialIndices( numBlocks );
  for(int b = 0; b < numBlocks; b++ )
  {
    ParticleBlock & particleBlock = particleBlockManager.registerParticleBlock( m_blockNames[b] );
    particleBlock.setParticleType( EnumStrings< ParticleType >::fromString( m_particleTypes[b] ) );
  }

  // Collect material names for each particle block
  // Assumes particleRegions loop follows material name ordering of input file
  int regionIndex = 0;    
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    auto & subRegions = particleRegion.getSubRegions();
    for( int r=0; r < subRegions.size(); ++r)
    {
      ParticleSubRegion & subRegion = dynamicCast< ParticleSubRegion & >( *subRegions[r] );
      GEOS_LOG_RANK( particleRegion.getName() << " | " << subRegion.getName() );
    }

    string_array particleBlockNames = particleRegion.getParticleBlockNames();
    for( auto & particleBlockName : particleBlockNames)
    {
      for( int b = 0; b < numBlocks; b++ )
      {
        if( particleBlockName == m_blockNames[b])
        {
          blockMaterialIndices[b] = regionIndex;                                                     
          break;
        }
      }
    }
    regionIndex++;
  } );

  std::string line; // initialize line variable
  std::ifstream particleFile( m_particleFilePath );
  std::vector< std::vector< std::vector< double > > > particleData( numBlocks );

  // Read column headers for particle data
  int numColumnHeaders = 0;
  std::map< int, int > columnHeaderMap;
  {
    std::getline( particleFile, line );
    std::string token;
  
    std::istringstream lineStream( line );
    while( std::getline( lineStream, token, '\t') )
    {
      columnHeaderMap.insert( std::pair< int, int >( static_cast< int >( EnumStrings< ParticleColumnHeaders >::fromString( token ) ), numColumnHeaders ) );
      numColumnHeaders++;
    }
  }

  // Read in particle data
  int lineNumber = 1; // Since column header takes one line
  while( particleFile.peek() != EOF )
  {
    std::getline( particleFile, line );
    std::istringstream lineStream( line );

    std::vector< double > lineData; // TODO: Not great because we cast all input as doubles, but it all gets re-cast later so maybe it's fine.

    // Read line from particle file and parse columns
    double value;
    int numLineColumns = 0;
    while( lineStream >> value )
    {
      lineData.push_back( value );
      numLineColumns++;
    }

    GEOS_ERROR_IF( numLineColumns != numColumnHeaders, "Particle file line " << lineNumber << " has a different number of terms than the column headers! Was " << numLineColumns << " but should be " << numColumnHeaders );
    lineNumber++;

    // If particle is inside partition add to particleData otherwise ignore and continue parsing file
    bool inPartition = partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionX ) ] ], 0 ) && 
                       partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionY ) ] ], 1 ) && 
                       partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionZ ) ] ], 2 );
    if( !inPartition )
    {
      continue;
    }

    // Reformat particle data and apply defaults to particle fields not specified in particle file
    std::vector< double > lineDataInside;
    // Presently get the number of enums using a sentinel (COUNT), but ideally this should be done more robustly
    for(int c = 0; c < EnumSize<ParticleColumnHeaders>; c++)
    {
      if( columnHeaderMap.find( c ) != columnHeaderMap.end() )
      {
        lineDataInside.push_back( lineData[ columnHeaderMap[ c ] ] );
        continue;
      }

      // Apply default value
      double defaultValue;
      switch( static_cast< ParticleColumnHeaders >( c ) )
      {
        case ParticleColumnHeaders::StrengthScale:
        case ParticleColumnHeaders::MaterialDirectionX:
        case ParticleColumnHeaders::SurfaceNormalX:
          defaultValue = 1.0;
          break;
        case ParticleColumnHeaders::Temperature:
          defaultValue = 300.0;
          break;
        case ParticleColumnHeaders::ParticleType:
          defaultValue = 2.0;
          break;
        case ParticleColumnHeaders::MaterialType:
          defaultValue = 0.0;
          break;
        case ParticleColumnHeaders::ContactGroup:
        case ParticleColumnHeaders::Damage:
        case ParticleColumnHeaders::Porosity:
        case ParticleColumnHeaders::VelocityX:
        case ParticleColumnHeaders::VelocityY:
        case ParticleColumnHeaders::VelocityZ:
        case ParticleColumnHeaders::MaterialDirectionY:
        case ParticleColumnHeaders::MaterialDirectionZ:
        case ParticleColumnHeaders::SurfaceNormalY:
        case ParticleColumnHeaders::SurfaceNormalZ:
        case ParticleColumnHeaders::SurfacePositionX:
        case ParticleColumnHeaders::SurfacePositionY:
        case ParticleColumnHeaders::SurfacePositionZ:
        case ParticleColumnHeaders::SurfaceTractionX:
        case ParticleColumnHeaders::SurfaceTractionY:
        case ParticleColumnHeaders::SurfaceTractionZ:
          defaultValue = 0.0;
          break;
        default:
          GEOS_ERROR( EnumStrings< ParticleColumnHeaders >::toString( static_cast< ParticleColumnHeaders >( c ) ) << " must be specified in particle file!" );
          break;
      }

      lineDataInside.push_back( defaultValue );
    }

    int materialTypeIndex = lineDataInside[static_cast< int >( ParticleColumnHeaders::MaterialType )];
    int particleTypeIndex = lineDataInside[static_cast< int >( ParticleColumnHeaders::ParticleType )];

    // Match particle to block using material and particle types
    int blockIndex = -1;
    for( int b = 0; b < numBlocks; b++)
    {
      if( materialTypeIndex == blockMaterialIndices[b] && 
          particleTypeIndex == static_cast< int >( EnumStrings< ParticleType >::fromString( m_particleTypes[b] ) ) )
      {
        blockIndex = b;
      }
    }

    GEOS_ERROR_IF(blockIndex < 0, "Particle at line  " << lineNumber <<  " with type and material indices of " << particleTypeIndex << " and " << materialTypeIndex << " respectively, does not match any particle block!");

    particleData[blockIndex].push_back( lineDataInside );
  }

  // Loop over blocks and assign particle attributes
  for( int b = 0; b < numBlocks; b++ )
  {
    std::string blockName = m_blockNames[b];
    ParticleBlock & particleBlock = particleBlockManager.getParticleBlock( blockName );

    int npInBlock = particleData[b].size();
    particleBlock.resize( npInBlock );
    array1d< globalIndex > particleID( npInBlock );
    array2d< real64 > particleCenter( npInBlock, 3 );
    array2d< real64 > particleVelocity( npInBlock, 3 );
    array2d< real64 > particleMaterialDirection( npInBlock, 3 );
    array1d< int > particleGroup( npInBlock );
    array1d< int > particleSurfaceFlag( npInBlock );
    array1d< real64 > particleDamage( npInBlock );
    array1d< real64 > particlePorosity( npInBlock );
    array1d< real64 > particleTemperature( npInBlock );
    array1d< real64 > particleVolume( npInBlock );
    array1d< real64 > particleStrengthScale( npInBlock );
    array3d< real64 > particleRVectors( npInBlock, 3, 3 ); // TODO: Flatten the r-vector array into a 1x9 for each particle
    array2d< real64 > particleSurfaceNormal( npInBlock, 3);
    array2d< real64 > particleSurfacePosition( npInBlock, 3 );
    array2d< real64 > particleSurfaceTraction( npInBlock, 3 );

    // Populate particle fields with data
    for( int i = 0; i < npInBlock; i++ )
    {
      // Global ID
      particleID[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::ID )];

      // Position
      particleCenter[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::PositionX )];
      particleCenter[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::PositionY )];
      particleCenter[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::PositionZ )];

      // Velocity
      particleVelocity[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::VelocityX )];
      particleVelocity[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::VelocityY )];
      particleVelocity[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::VelocityZ )];
  
      // Material (set above) is [10]

      // Group
      particleGroup[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::ContactGroup )];

      // surfaceFlag
      particleSurfaceFlag[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceFlag )];

      // Damage
      particleDamage[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::Damage )];

      // Porosity
      particlePorosity[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::Porosity )];

      // Temperature
      particleTemperature[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::Temperature )];

      // strengthScale
      particleStrengthScale[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::StrengthScale )];

      // Volume and R-Vectors
      switch( particleBlock.getParticleType() )
      {
        case ParticleType::SinglePoint:
        {
          particleVolume[i] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXX )];

          double a = std::pow( particleVolume[i], 1.0/3.0 );
          particleRVectors[i][0][0] = a;
          particleRVectors[i][0][1] = 0.0;
          particleRVectors[i][0][2] = 0.0;
          particleRVectors[i][1][0] = 0.0;
          particleRVectors[i][1][1] = a;
          particleRVectors[i][1][2] = 0.0;
          particleRVectors[i][2][0] = 0.0;
          particleRVectors[i][2][1] = 0.0;
          particleRVectors[i][2][2] = a;
          break;
        }
        case ParticleType::SinglePointBSpline:
        {
          GEOS_ERROR("SinglePointBSpline particle type is not implemented!");
        }
        case ParticleType::CPDI:
        {
           double x1, y1, z1, x2, y2, z2, x3, y3, z3;
          x1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXX )];
          y1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXY )];
          z1 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorXZ )];
          x2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYX )];
          y2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYY )];
          z2 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorYZ )];
          x3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZX )];
          y3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZY )];
          z3 = particleData[b][i][static_cast< int >( ParticleColumnHeaders::RVectorZZ )];
          particleRVectors[i][0][0] = x1;
          particleRVectors[i][0][1] = y1;
          particleRVectors[i][0][2] = z1;
          particleRVectors[i][1][0] = x2;
          particleRVectors[i][1][1] = y2;
          particleRVectors[i][1][2] = z2;
          particleRVectors[i][2][0] = x3;
          particleRVectors[i][2][1] = y3;
          particleRVectors[i][2][2] = z3;
          particleVolume[i] = 8.0*std::fabs( -(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3 );
          break;
        }
        case ParticleType::CPTI:
        {
          GEOS_ERROR("CPTI particle type is not implemented!");
        }
        case ParticleType::CPDI2:
        {
          GEOS_ERROR("CPDI2 particle type is not implemented!");
        }
        default:
        {
          GEOS_ERROR( "Invalid particle type: " << particleBlock.getParticleType() );
        }
      }

      // Material Direction
      particleMaterialDirection[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionX )];
      particleMaterialDirection[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionY )];
      particleMaterialDirection[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::MaterialDirectionZ )];

      particleSurfaceNormal[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalX )];
      particleSurfaceNormal[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalY )];
      particleSurfaceNormal[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceNormalZ )];

      particleSurfacePosition[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionX )];
      particleSurfacePosition[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionY )];
      particleSurfacePosition[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfacePositionZ )];

      particleSurfaceTraction[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionX )];
      particleSurfaceTraction[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionY )];
      particleSurfaceTraction[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::SurfaceTractionZ )];
    }

    particleBlock.setParticleID( particleID );
    particleBlock.setParticleCenter( particleCenter );
    particleBlock.setParticleVelocity( particleVelocity );
    particleBlock.setParticleMaterialDirection( particleMaterialDirection );
    particleBlock.setParticleGroup( particleGroup );
    particleBlock.setParticleSurfaceFlag( particleSurfaceFlag );
    particleBlock.setParticleDamage( particleDamage );
    particleBlock.setParticlePorosity( particlePorosity );
    particleBlock.setParticleTemperature( particleTemperature );
    particleBlock.setParticleStrengthScale( particleStrengthScale );
    particleBlock.setParticleVolume( particleVolume );
    particleBlock.setParticleRVectors( particleRVectors );
    particleBlock.setParticleSurfaceNormal( particleSurfaceNormal );
    particleBlock.setParticleSurfacePosition( particleSurfacePosition );
    particleBlock.setParticleSurfaceTraction( particleSurfaceTraction );
  } // loop over particle blocks

  // Resize particle regions
  int numParticles = 0;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    string_array particleBlockNames = particleRegion.getParticleBlockNames();

    int size = 0; 
    for( auto & particleBlockName : particleBlockNames )
    {
      for( int b = 0; b < numBlocks; b++ )
      {
        if( particleBlockName == m_blockNames[b] )
        {
          size += particleData[b].size();
          break;
        }
      }
    }
    numParticles += size;
    particleRegion.resize( size );
  } );

  particleManager.resize( numParticles ); // All this does is change m_size for the particleManager, gives a convenient way to get the total
                                          // number of particles
  GEOS_LOG_RANK( "Total number of particles on this rank: " << particleManager.size() );
}

REGISTER_CATALOG_ENTRY( ParticleMeshGeneratorBase, ParticleMeshGenerator, string const &, Group * const )

} /* namespace geos */
