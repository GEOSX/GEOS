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

#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/PartitionBase.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "ParticleBlockManager.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

#include <cmath>

namespace geos
{
using namespace dataRepository;

ParticleMeshGenerator::ParticleMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_dim( 3 ),
  m_min(),
  m_max()
{
  registerWrapper( viewKeyStruct::particleFilePathString(), &m_particleFilePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the particle file" );

  registerWrapper( viewKeyStruct::headerFilePathString(), &m_headerFilePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the header file" );

  registerWrapper( viewKeyStruct::particleBlockNamesString(), &m_blockNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of each particle block" );

  registerWrapper( viewKeyStruct::particleTypesString(), &m_particleType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Particle types of each particle block" );
}

Group * ParticleMeshGenerator::createChild( string const & GEOS_UNUSED_PARAM( childKey ),
                                            string const & GEOS_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

/**
 * @param partition
 * @param domain
 */
void ParticleMeshGenerator::generateMesh( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  MeshBody & meshBody = domain.getMeshBody( this->getName() );
  //GEOS_LOG_RANK( "Generating mesh for " + meshBody.getName() );
  MeshLevel & meshLevel0 = meshBody.getBaseDiscretization();
  ParticleManager & particleManager = meshLevel0.getParticleManager();

  ParticleBlockManager & particleBlockManager = meshBody.registerGroup< ParticleBlockManager >( keys::particleManager );

  // SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );
  SpatialPartition & partition = dynamic_cast< SpatialPartition & >( domain.getPartition() );

  // This should probably handled elsewhere:
  int aa = 0;

  for( auto & particleBlockName : m_blockNames )
  {
    ParticleBlock & particleBlock = particleBlockManager.registerParticleBlock( particleBlockName );
    particleBlock.setParticleType( EnumStrings< ParticleType >::fromString( m_particleType[aa++] ) );
  }

  //GEOS_LOG_RANK_0( "MPM particle file path: " << m_particleFilePath );
  //GEOS_LOG_RANK_0( "MPM header file path: " << m_headerFilePath );

  int numMaterials, numParticleTypes;
  map< std::string, std::vector< std::vector< double > > > particleData;
  map< std::string, int > particleTypeMap;
  std::vector< std::string > particleTypes; // This is needed because the input file format is such that data associated with each particle
                                            // type is in the same order as the preceding type listing.
                                            // Looping over the particleTypeMap with an iterator does not respect this ordering since a map
                                            // automatically sorts by its keys.
  map< std::string, int > materialMap;

  // Get and process header and particle files
  std::ifstream headerFile( m_headerFilePath );
  std::ifstream particleFile( m_particleFilePath );
  std::string line; // initialize line variable

  // Read in number of materials and particle types
  std::getline( headerFile, line ); // get a line
  std::istringstream iss1( line ); // turn the line into a stream
  iss1 >> numMaterials >> numParticleTypes;
  particleTypes.resize( numParticleTypes );

  //GEOS_LOG_RANK_0( "Number of particle materials: " << numMaterials );
  //GEOS_LOG_RANK_0( "Number of particle types: " << numParticleTypes );

  // Read in the material key
  for( int i=0; i<numMaterials; i++ )
  {
    std::getline( headerFile, line );
    std::istringstream iss2( line );
    std::string key; // Material name
    int value; // Material ID
    iss2 >> key >> value;
    materialMap[key] = value;
    //GEOS_LOG_RANK_0( "Material name/ID: " + key + "/" << value );
  }

  // Read in the particle type key
  for( int i=0; i<numParticleTypes; i++ )
  {
    std::getline( headerFile, line );
    std::istringstream iss2( line );
    std::string particleType; // Particle type
    int np; // Number of particles of that type
    iss2 >> particleType >> np;
    particleTypeMap[particleType] = np;
    particleTypes[i] = particleType;
  }

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
  int lineNumber = 1; // Since colum header takes one line
  for( size_t i=0; i < particleTypes.size(); i++ )
  {
    for( int j=0; j < particleTypeMap[particleTypes[i]]; j++ )
    {
      std::getline( particleFile, line );
      std::istringstream lineStream( line );

      std::vector< double > lineData; // TODO: Not great because we cast all input as doubles, but it all gets re-cast later so maybe it's fine.

      // Read line from particle file and parse columns
      double value;
      int numColumns = 0;
      while( lineStream >> value )
      {
        lineData.push_back( value );
        numColumns++;
      }

      GEOS_ERROR_IF( numColumns != numColumnHeaders, "Particle file line " << lineNumber << " has a different number of terms than the column headers!" );

      lineNumber++;

      // If particle is inside partition add to particleData otherwise ignore and continue parsing file
      bool inPartition = partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionX ) ] ], 0 ) && 
                         partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionY ) ] ], 1 ) && 
                         partition.isCoordInPartition( lineData[ columnHeaderMap[ static_cast< int >( ParticleColumnHeaders::PositionZ ) ] ], 2 );
      if( !inPartition )
      {
        continue;
      }

      // Reformat particle data and apply defaults to fields not specified
      std::vector< double > lineDataInside;
      // CC: TODO: Can you get the number of options from the enum directly? maybe JA
      for(int c = 0; c < static_cast<int>(ParticleColumnHeaders::COUNT); c++)
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
          case ParticleColumnHeaders::MaterialType:
          case ParticleColumnHeaders::ContactGroup:
          case ParticleColumnHeaders::VelocityX:
          case ParticleColumnHeaders::VelocityY:
          case ParticleColumnHeaders::VelocityZ:
          case ParticleColumnHeaders::MaterialDirectionY:
          case ParticleColumnHeaders::MaterialDirectionZ:
          case ParticleColumnHeaders::SurfaceNormalY:
          case ParticleColumnHeaders::SurfaceNormalZ:
          case ParticleColumnHeaders::Damage:
            defaultValue = 0.0;
            break;
          case ParticleColumnHeaders::AccelerationX:
            defaultValue = 0.0;
            break;
          case ParticleColumnHeaders::AccelerationY:
            defaultValue = 0.0;
            break;
          case ParticleColumnHeaders::AccelerationZ:
            defaultValue = 0.0;
            break;
          default:
            GEOS_ERROR( EnumStrings< ParticleColumnHeaders >::toString( static_cast< ParticleColumnHeaders >( c ) ) << " must be specified in particle file!" );
            break;
        }

        lineDataInside.push_back( defaultValue );
      }

      particleData[particleTypes[i]].push_back( lineDataInside );

    }
  }

  // Construct the map from particle blocks to particle regions (more specifically, the regions' associated materials)
  map< std::string, int > blockMaterialMap;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    string_array particleBlockNames = particleRegion.getParticleBlockNames();
    std::string material = particleRegion.getMaterialList()[0]; // We will assume that the material list for a region contains only one
                                                                // material since MPM will only be doing single phase mechanics for now
    for( auto i=0; i < particleBlockNames.size(); i++ )
    {
      blockMaterialMap[particleBlockNames[i]] = materialMap[material];
    }
  } );

  // Distribute particle information to particle blocks
  map< std::string, std::vector< int > > indexMap; // This will keep track of the indices of particleData associated with each particle
                                                   // block. It's populated in the loop that checks for which particles belong to a block.
  map< std::string, int > sizeMap;  // This keeps track of the size of each particle block so we can resize the ParticleRegions later
  for( auto & particleBlockName : m_blockNames )
  {
    ParticleBlock & particleBlock = particleBlockManager.getParticleBlock( particleBlockName );
    std::string particleType = EnumStrings< ParticleType >::toString( particleBlock.getParticleType() );

    int numThisType = particleData[particleType].size(); // I hope this returns zero when particleType isn't in the map...
    int materialID;
    int npInBlock = 0; // Number of particles in this particle block
    for( localIndex i=0; i<numThisType; i++ ) // Find out which particles belong to the current particle block
    {
      materialID = particleData[particleType][i][7]; // The particle file is configured such that the 11th column has the material ID
      if( materialID == blockMaterialMap[particleBlock.getName()] )
      {
        npInBlock++;
        indexMap[particleBlockName].push_back( i );
      }
    }
    sizeMap[particleBlock.getName()] = npInBlock;
    particleBlock.resize( npInBlock );

    array1d< globalIndex > particleID( npInBlock );
    array2d< real64 > particleCenter( npInBlock, 3 );
    array2d< real64 > particleVelocity( npInBlock, 3 );
    array2d< real64 > particleAcceleration( npInBlock, 3 );
    array2d< real64 > particleMaterialDirection( npInBlock, 3 );
    array1d< int > particleGroup( npInBlock );
    array1d< int > particleSurfaceFlag( npInBlock );
    array1d< real64 > particleDamage( npInBlock );
    array1d< real64 > particleVolume( npInBlock );
    array1d< real64 > particleStrengthScale( npInBlock );
    array3d< real64 > particleRVectors( npInBlock, 3, 3 ); // TODO: Flatten the r-vector array into a 1x9 for each particle
    array2d< real64 > particleSurfaceNormal( npInBlock, 3); // TODO:: read from file eventually

    // Assign particle data to the appropriate block.
    std::vector< int > & indices = indexMap[particleBlockName];
    int index = 0;
    for( int i : indices )
    {
      // Global ID
      particleID[index] = particleData[particleType][i][0];

      // Position
      particleCenter[index][0] = particleData[particleType][i][1];
      particleCenter[index][1] = particleData[particleType][i][2];
      particleCenter[index][2] = particleData[particleType][i][3];

      // Velocity
      particleVelocity[index][0] = particleData[particleType][i][4];
      particleVelocity[index][1] = particleData[particleType][i][5];
      particleVelocity[index][2] = particleData[particleType][i][6];

      // Acceleration
      particleAcceleration[i][0] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::AccelerationX )];
      particleAcceleration[i][1] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::AccelerationY )];
      particleAcceleration[i][2] = particleData[b][i][static_cast< int >( ParticleColumnHeaders::AccelerationZ )];
  
      // Material (set above) is [7]

      // Group
      particleGroup[index] = particleData[particleType][i][8];

      // surfaceFlag
      particleSurfaceFlag[index] = particleData[particleType][i][9];

      // Damage
      particleDamage[index] = particleData[particleType][i][10];

      // strengthScale
      particleStrengthScale[index] = particleData[particleType][i][11];

      // Volume and R-Vectors
      if( particleType == "SinglePoint" )
      {
        particleVolume[index] = particleData[particleType][i][12];
        double a = std::pow( particleVolume[index], 1.0/3.0 );
        particleRVectors[index][0][0] = a;
        particleRVectors[index][0][1] = 0.0;
        particleRVectors[index][0][2] = 0.0;
        particleRVectors[index][1][0] = 0.0;
        particleRVectors[index][1][1] = a;
        particleRVectors[index][1][2] = 0.0;
        particleRVectors[index][2][0] = 0.0;
        particleRVectors[index][2][1] = 0.0;
        particleRVectors[index][2][2] = a;
      }
      else if( particleType == "CPDI" )
      {
        double x1, y1, z1, x2, y2, z2, x3, y3, z3;
        x1 = particleData[particleType][i][12];
        y1 = particleData[particleType][i][13];
        z1 = particleData[particleType][i][14];
        x2 = particleData[particleType][i][15];
        y2 = particleData[particleType][i][16];
        z2 = particleData[particleType][i][17];
        x3 = particleData[particleType][i][18];
        y3 = particleData[particleType][i][19];
        z3 = particleData[particleType][i][20];
        particleRVectors[index][0][0] = x1;
        particleRVectors[index][0][1] = y1;
        particleRVectors[index][0][2] = z1;
        particleRVectors[index][1][0] = x2;
        particleRVectors[index][1][1] = y2;
        particleRVectors[index][1][2] = z2;
        particleRVectors[index][2][0] = x3;
        particleRVectors[index][2][1] = y3;
        particleRVectors[index][2][2] = z3;
        particleVolume[index] = 8.0*std::fabs( -(x3*y2*z1) + x2*y3*z1 + x3*y1*z2 - x1*y3*z2 - x2*y1*z3 + x1*y2*z3 );
      }
      else
      {
        GEOS_ERROR( "Invalid particle type specification! Cannot determine particle volume, aborting." );
      }

      // Material Direction
      particleMaterialDirection[index][0] = particleData[particleType][i][21];
      particleMaterialDirection[index][1] = particleData[particleType][i][22];
      particleMaterialDirection[index][2] = particleData[particleType][i][23];

      // Surface Normal ( currently only works assuming that particles are formatted for CPDI data)
      // CC: TODO eventually include header in particle file for column data
      particleSurfaceNormal[index][0] = particleData[particleType][i][24];
      particleSurfaceNormal[index][1] = particleData[particleType][i][25];
      particleSurfaceNormal[index][2] = particleData[particleType][i][26];

      // Increment index
      index++;
    }
    particleBlock.setParticleID( particleID );
    particleBlock.setParticleCenter( particleCenter );
    particleBlock.setParticleVelocity( particleVelocity );
    particleBlock.setParticleAcceleration( particleAcceleration );
    particleBlock.setParticleInitialMaterialDirection( particleMaterialDirection );
    particleBlock.setParticleMaterialDirection( particleMaterialDirection );
    particleBlock.setParticleGroup( particleGroup );
    particleBlock.setParticleSurfaceFlag( particleSurfaceFlag );
    particleBlock.setParticleDamage( particleDamage );
    particleBlock.setParticleStrengthScale( particleStrengthScale );
    particleBlock.setParticleVolume( particleVolume );
    particleBlock.setParticleRVectors( particleRVectors );
    particleBlock.setParticleInitialSurfaceNormal( particleSurfaceNormal );
    particleBlock.setParticleSurfaceNormal( particleSurfaceNormal );
  } // loop over particle blocks

  // Resize particle regions
  int numParticles = 0;
  particleManager.forParticleRegions< ParticleRegion >( [&]( auto & particleRegion )
  {
    string_array particleBlockNames = particleRegion.getParticleBlockNames();
    std::string material = particleRegion.getMaterialList()[0]; // We will assume that the material list for a region contains only one
                                                                // material since MPM will only be doing single phase mechanics for now
    int size = 0;
    for( auto i=0; i<particleBlockNames.size(); i++ )
    {
      size += sizeMap[particleBlockNames[i]];
    }
    numParticles += size;
    particleRegion.resize( size );
  } );

  particleManager.resize( numParticles ); // All this does is change m_size for the particleManager, gives a convenient way to get the total
                                          // number of particles
  //GEOS_LOG_RANK( "Total number of particles on this rank: " << particleManager.size() );
}

void ParticleMeshGenerator::postProcessInput()
{
  //GEOS_LOG_RANK_0( "Someone called ParticleMeshGenerator::postProcessInput!" );
}

void ParticleMeshGenerator::importFieldOnArray( Block block,
                                                string const & blockName,
                                                string const & meshFieldName,
                                                bool isMaterialField,
                                                dataRepository::WrapperBase & wrapper ) const
{
  GEOS_UNUSED_VAR( block );
  GEOS_UNUSED_VAR( blockName );
  GEOS_UNUSED_VAR( meshFieldName );
  GEOS_UNUSED_VAR( isMaterialField );
  GEOS_UNUSED_VAR( wrapper );
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, ParticleMeshGenerator, string const &, Group * const )

} /* namespace geos */
