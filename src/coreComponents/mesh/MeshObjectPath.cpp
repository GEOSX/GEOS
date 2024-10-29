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
 * @file MeshObjectPath.cpp
 */

#include "MeshObjectPath.hpp"
#include "MeshBody.hpp"

#include <fnmatch.h>

namespace geos
{

namespace
{
MeshObjectPath::ObjectTypes extractObjectType( string const & path )
{
  MeshObjectPath::ObjectTypes objectType = MeshObjectPath::ObjectTypes::invalid;
  if( path.find( MeshLevel::groupStructKeys::nodeManagerString() ) != std::string::npos )
  {
    objectType = MeshObjectPath::ObjectTypes::nodes;
  }
  else if( path.find( MeshLevel::groupStructKeys::edgeManagerString() ) != std::string::npos )
  {
    objectType = MeshObjectPath::ObjectTypes::edges;
  }
  else if( path.find( MeshLevel::groupStructKeys::faceManagerString() ) != std::string::npos )
  {
    objectType = MeshObjectPath::ObjectTypes::faces;
  }
  else if( path.find( MeshLevel::groupStructKeys::elemManagerString() ) != std::string::npos )
  {
    objectType = MeshObjectPath::ObjectTypes::elems;
  }
  return objectType;
}

}

MeshObjectPath::MeshObjectPath( string const path,
                                dataRepository::Group const & meshBodies ):
  m_objectType( extractObjectType( path ) ),
  m_pathPermutations()
{
  processPath( path, meshBodies );
}


std::vector< string >
MeshObjectPath::fillPathTokens( string const & path,
                                dataRepository::Group const & meshBodies ) const
{
  std::vector< string > pathTokens = stringutilities::tokenize( path, "/" );

  // find where the object specification is in the path
  auto findObjectIndex = [&]() -> int
  {
    for( size_t a=0; a<pathTokens.size(); ++a )
    {
      if( pathTokens[a] == EnumStrings< ObjectTypes >::toString( m_objectType ) )
      {
        return a;
      }
    }
    return -1;
  };

  int objectIndex = findObjectIndex();

  GEOS_THROW_IF( objectIndex==-1,
                 GEOS_FMT( "Path {} does not contain a valid object type. "
                           "It must contain one of the following: {}, {}, {}, {}",
                           path,
                           MeshLevel::groupStructKeys::nodeManagerString(),
                           MeshLevel::groupStructKeys::edgeManagerString(),
                           MeshLevel::groupStructKeys::faceManagerString(),
                           MeshLevel::groupStructKeys::elemManagerString() ),
                 InputError );

  // No MeshBody or MeshLevels were specified. add all of them
  if( objectIndex==0 )
  {
    pathTokens.insert( pathTokens.begin(), "{*}" );
    pathTokens.insert( pathTokens.begin(), "{*}" );
  }
  // MeshBody OR MeshLevel specified. Check which one, and add all of the other.
  else if( objectIndex==1 )
  {
    string const unidentifiedName = pathTokens[0];
    // if the MeshBody is specified, add all MeshLevels
    if( meshBodies.hasGroup( unidentifiedName ) )
    {
      pathTokens.insert( pathTokens.begin(), unidentifiedName );
      pathTokens.insert( pathTokens.begin()+1, "{*}" );
    }
    // It wasn't the MeshBody that was specified, it was the MeshLevel?? Check and add.
    else
    {
      pathTokens.insert( pathTokens.begin(), "{*}" );

      // searching if the mesh level exists
      bool levelNameFound = false;
      meshBodies.forSubGroups< MeshBody >( [&]( MeshBody const & meshBody )
      {
        meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
        {
          levelNameFound |= ( unidentifiedName==meshLevel.getName() );
        } );
      } );

      if( !levelNameFound )
      {
        string existingMeshBodiesAndLevels;
        meshBodies.forSubGroups< MeshBody >( [&]( MeshBody const & meshBody )
        {
          std::vector< string > meshLevelsNames;
          existingMeshBodiesAndLevels += "  MeshBody "+meshBody.getName() + ": { ";
          meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
          {
            meshLevelsNames.push_back( meshLevel.getName() );
          } );
          existingMeshBodiesAndLevels += stringutilities::join( meshLevelsNames, ", " ) + " }\n";
        } );

        GEOS_THROW( GEOS_FMT( "Path {0} specifies an invalid MeshBody or MeshLevel. ",
                              "existing MeshBodies: \n{1}\n",
                              path,
                              existingMeshBodiesAndLevels ),
                    InputError );
      }
      pathTokens.insert( pathTokens.begin()+1, unidentifiedName );
    }
  }

  // checks to make sure we have properly inserted thus far
  objectIndex = findObjectIndex();
  size_t targetTokenLength = pathTokens.size();

  GEOS_THROW_IF_NE_MSG( objectIndex, 2,
                        "Filling of MeshBody and/or MeshLevel in path has failed. Object Index should be 2",
                        InputError );

  GEOS_THROW_IF( targetTokenLength < 2,
                 "Filling of MeshBody and/or MeshLevel in path has failed. targetTokenLength should be greater than 2",
                 InputError );

  // now we need to fill in any missing region/subregion specifications.

  if( m_objectType == ObjectTypes::elems )
  {
    // there are no regions specified
    if( targetTokenLength == 3 )
    {
      pathTokens.push_back( "{*}" );
      pathTokens.push_back( "{*}" );
    }
    // there are no subregions specified
    else if( targetTokenLength == 4 )
    {
      pathTokens.push_back( "{*}" );
    }
  }


  for( size_t a=0; a<pathTokens.size(); ++a )
  {
    pathTokens[a].erase( std::remove( pathTokens[a].begin(), pathTokens[a].end(), '{' ), pathTokens[a].end());
    pathTokens[a].erase( std::remove( pathTokens[a].begin(), pathTokens[a].end(), '}' ), pathTokens[a].end());
  }

  return pathTokens;
}


template< typename SUBNODE >
static SUBNODE & insertPathNode( std::map< string, SUBNODE > & node, string const & name )
{
  return node[ name ];
}

static string & insertPathNode( std::vector< string > & node, string & name )
{
  node.push_back( name );
  return name;
}


template< typename TYPE, typename NODETYPE, typename CALLBACK >
void processTokenRecursive( dataRepository::Group const & parentGroup,
                            string const & pathToken,
                            NODETYPE & node,
                            CALLBACK && cbfunc )
{
  std::vector< string > namesInRepository;
  parentGroup.forSubGroups< TYPE >( [&]( TYPE const & group )
  {
    namesInRepository.emplace_back( group.getName() );
  } );

  GEOS_THROW_IF( namesInRepository.empty(),
                 GEOS_FMT( "{0} has no children.", parentGroup.getDataContext().toString()),
                 InputError );

  for( string const & inputEntry : stringutilities::tokenize( pathToken, " " ) )
  {
    bool foundMatch = false;
    for( string const & candidateName : namesInRepository )
    {
      string name = candidateName;
      int const fnmatchResult = fnmatch( inputEntry.c_str(), candidateName.c_str(), 0 );
      if( fnmatchResult != FNM_NOMATCH )
      {
        auto & subNode = insertPathNode( node, name );
        foundMatch=true;
        // recursive call
        cbfunc( parentGroup.getGroup< TYPE >( name ),
                subNode );

      }
    }
    GEOS_THROW_IF( !foundMatch,
                   GEOS_FMT( "{0} has no child named {1}.\n"
                             "{0} has the following children: {{ {2} }}",
                             parentGroup.getDataContext().toString(),
                             inputEntry,
                             stringutilities::join( namesInRepository, ", " ) ),
                   InputError );
  }
}


void MeshObjectPath::processPathTokens( std::vector< string > const & pathTokens,
                                        dataRepository::Group const & meshBodies )
{

  processTokenRecursive< MeshBody >( meshBodies,
                                     pathTokens[0],
                                     m_pathPermutations,
                                     [this, &pathTokens] ( MeshBody const & meshBody,
                                                           std::map< string, std::map< string, std::vector< string > > > & meshBodyNode )
  {
    dataRepository::Group const & meshLevels = meshBody.getMeshLevels();
    processTokenRecursive< MeshLevel >( meshLevels,
                                        pathTokens[1],
                                        meshBodyNode,
                                        [this, &pathTokens]( MeshLevel const & meshLevel,
                                                             std::map< string, std::vector< string > > & meshLevelNode )
    {
      if( m_objectType == ObjectTypes::elems )
      {
        dataRepository::Group const & elemRegionGroup = meshLevel.getElemManager().getGroup( ElementRegionManager::groupKeyStruct::elementRegionsGroup() );
        processTokenRecursive< ElementRegionBase >( elemRegionGroup,
                                                    pathTokens[3],
                                                    meshLevelNode,
                                                    [&]( ElementRegionBase const & elemRegion,
                                                         std::vector< string > & elemRegionNode )
        {
          dataRepository::Group const & elemSubRegionGroup = elemRegion.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() );
          processTokenRecursive< ElementSubRegionBase >( elemSubRegionGroup,
                                                         pathTokens[4],
                                                         elemRegionNode,
                                                         [&]( ElementSubRegionBase const &,
                                                              string & )
          {} );
        } );
      }
    } );
  } );
}


void MeshObjectPath::processPath( string const objectPath,
                                  dataRepository::Group const & meshBodies )
{
  std::vector< string > pathTokens = fillPathTokens( objectPath, meshBodies );
  processPathTokens( pathTokens, meshBodies );
}

void MeshObjectPath::printPermutations() const
{
  for( auto const & meshBodyPair : m_pathPermutations )
  {
    std::cout<<meshBodyPair.first<<": "<<std::endl;
    for( auto const & meshLevelPair : meshBodyPair.second )
    {
      std::cout<<"  "<<meshLevelPair.first<<": "<<std::endl;
      for( auto const & elemRegionPair : meshLevelPair.second )
      {
        std::cout<<"    "<<elemRegionPair.first<<": "<<std::endl;
        std::cout<<"      ";
        for( auto const & elemSubRegionName : elemRegionPair.second )
        {
          std::cout<<elemSubRegionName<<", ";
        }
        std::cout<<std::endl;
      }
    }
  }
  std::cout<<std::endl<<std::endl;
}

bool MeshObjectPath::containsMeshLevel( MeshLevel const & meshLevel ) const
{
  bool isMeshLevelInObjectPath = false;
  string const bodyName = meshLevel.getParent().getParent().getName();
  string const levelName = meshLevel.getName();

  auto bodyIter = m_pathPermutations.find( bodyName );
  if( bodyIter != m_pathPermutations.end() )
  {
    auto const levelIter = bodyIter->second.find( levelName );
    isMeshLevelInObjectPath = ( levelIter != bodyIter->second.end() );
  }
  return isMeshLevelInObjectPath;
}

} /* namespace geos */
