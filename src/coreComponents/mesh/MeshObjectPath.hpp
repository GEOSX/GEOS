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
 * @file MeshObjectPath.hpp
 */

#ifndef GEOS_MESH_MESHOBJECTPATH_HPP_
#define GEOS_MESH_MESHOBJECTPATH_HPP_


#include "common/format/EnumStrings.hpp"
#include "MeshLevel.hpp"

namespace geos
{
class MeshBody;
class MeshLevel;

/**
 * @brief Class to hold the path to a collection of mesh objects
 *
 */
class MeshObjectPath
{
public:

  /**
   * @brief The container type that holds the path information
   * The first key is the name of a MeshBody
   * The second key is the name of a MeshLevel
   * The third key is the name of an ElementRegion
   * The third value is a vector of subregion names
   */
  using permutationMapType = std::map< string, std::map< string, std::map< string, std::vector< string > > > >;

  /**
   * @brief Contains enums for the types of objects
   */
  enum class ObjectTypes : int
  {
    nodes,  ///< a NodeManager
    edges,  ///< an EdgeManager
    faces,  ///< a FaceManager
    elems,  ///< a ElementManager
    invalid ///< an invalide object
  };

  /**
   * @brief Construct a new Mesh Object Path object
   *
   * @param path The path string
   * @param meshBodies  The Group that contains all MeshBody objects
   * @throw InputError when the input path is wrong.
   */
  MeshObjectPath( string const path,
                  dataRepository::Group const & meshBodies );

  /**
   * @brief Processes the path string into the permutation container
   *
   * @param path The path string
   * @param meshBodies  The Group that contains all MeshBody objects
   */
  void processPath( string const path,
                    dataRepository::Group const & meshBodies );

  /**
   * @brief Get the Object Type object
   * @return ObjectTypes const&
   */
  ObjectTypes const & getObjectType() const
  {
    return m_objectType;
  }

  /**
   * @brief Get the m_pathPermutations object
   *
   * @return permutationMapType const&
   */
  permutationMapType const & pathPermutations() const
  {
    return m_pathPermutations;
  }

  /**
   * @brief Helper function to decide whether a given meshLevel is in the objectPath
   * @param[in] meshLevel the mesh level that we want to search for in the objectPath
   * @return true if the meshLevel is in the objectPath, false otherwise
   * @details An example use case is in the validation of boundary conditions
   */
  bool containsMeshLevel( MeshLevel const & meshLevel ) const;

  /**
   * @brief LLoop over objects in the path and execute a callback function.
   *
   * @tparam OBJECT_TYPE The type of object to loop over
   * @tparam FUNC The type of function that is executed on the OBJECT_TYPE. Takes a
   *  single OBJECT_TYPE as an argument.
   *  func( dynamic_cast< OBJECT_TYPE & >(object) );
   * @param meshBodies Group that contains the MeshBody objects.
   * @param func The function that is executed on the OBJECT_TYPE
   */
  template< typename OBJECT_TYPE = dataRepository::Group,
            typename FUNC >
  void forObjectsInPath( dataRepository::Group & meshBodies,
                         FUNC && func ) const;

  /**
   * @brief Loop over objects in the path and execute a callback function.
   *
   * @tparam OBJECT_TYPE The type of object to loop over
   * @tparam FUNC The type of function that is executed on the OBJECT_TYPE. Takes a
   *  single OBJECT_TYPE as an argument.
   *  func( dynamic_cast< OBJECT_TYPE & >(object) );
   * @param meshBodies Group that contains the MeshBody objects.
   * @param func The function that is executed on the OBJECT_TYPE
   */
  template< typename OBJECT_TYPE = dataRepository::Group,
            typename FUNC >
  void forObjectsInPath( dataRepository::Group const & meshBodies,
                         FUNC && func ) const;

  /**
   * @brief Loop over objects in the path and execute a callback function.
   *
   * @tparam OBJECT_TYPE The type of object to loop over
   * @tparam FUNC The type of function that is executed on the OBJECT_TYPE. Takes a
   *  single OBJECT_TYPE as an argument.
   *  func( dynamic_cast< OBJECT_TYPE & >(object) );
   * @param level The MeshLevel that contains OBJECT_TYPE to be executed on.
   * @param func The function that is executed on the OBJECT_TYPE
   */
  template< typename OBJECT_TYPE = dataRepository::Group,
            typename FUNC >
  void forObjectsInPath( MeshLevel & level, FUNC && func ) const;

#if defined(MESH_OBJECT_PATH_PRIVATE_FUNCTION_UNIT_TESTING)
  template< typename OBJECT_TYPE >
  bool testCheckObjectTypeConsistency()
  {
    return checkObjectTypeConsistency< OBJECT_TYPE >();
  }

  std::vector< string > testFillPathTokens( string const & path,
                                            dataRepository::Group const & meshBodies )
  {
    return fillPathTokens( path, meshBodies );
  }

#endif

private:

  /**
   * @brief Loop over objects in the path and execute a callback function.
   * @tparam OBJECT_TYPE The type of object to loop over
   * @tparam FUNC The type of function that is executed on the OBJECT_TYPE. Takes a
   *  single OBJECT_TYPE as an argument.
   *  func( dynamic_cast< OBJECT_TYPE & >(object) );
   * @param levelPair an entry for a given level extracted from m_pathPermutations.
   * @param level The MeshLevel that contains OBJECT_TYPE to be executed on.
   * @param func The function that is executed on the OBJECT_TYPE
   */
  template< typename OBJECT_TYPE,
            typename FUNC >
  void forObjectsInPath( std::pair< string const, std::map< string, std::vector< string > > > const & levelPair,
                         MeshLevel & meshLevel,
                         FUNC && func ) const;

  template< typename OBJECT_TYPE,
            typename FUNC >
  void forObjectsInPath( std::pair< string const, std::map< string, std::vector< string > > > const & levelPair,
                         MeshLevel const & meshLevel,
                         FUNC && func ) const;

  /**
   * @brief A logical check for whether or not the m_objecType is consistent
   *  with a specific OBJECT_TYPE
   * @tparam OBJECT_TYPE The type to check m_objectType against.
   * @return true If OBJECT_TYPE is the same type of a base of the type implied by m_objectType.
   * @return false If OBJECT_TYPE is NOT the same type of a base of the type implied by m_objectType.
   */
  template< typename OBJECT_TYPE >
  bool checkObjectTypeConsistency() const;

  /**
   * @brief prints the contents of m_pathPermutations for debugging
   */
  void printPermutations() const;

  /**
   * @brief Create a tokenized version of the path
   * @param path The input path
   * @param meshBodies The Group that contains the MeshBody objects on the domain
   * @return std::vector< string >  A tokenized representation of the path.
   */
  std::vector< string > fillPathTokens( string const & path,
                                        dataRepository::Group const & meshBodies ) const;

  /**
   * @brief Convert the tokenized path into a collection of permutations and fill
   *  m_pathPermutations.
   * @param pathTokens The tokenized path
   * @param meshBodies The Group that contains the MeshBody objects on the domain
   */
  void processPathTokens( std::vector< string > const & pathTokens,
                          dataRepository::Group const & meshBodies );



  /// The type ObjectType for this path
  ObjectTypes const m_objectType;

  /// The path container
  permutationMapType m_pathPermutations;

};

} /* namespace geos */

#include "MeshBody.hpp"

namespace geos
{

/**
 * @brief Create conversion functions for ObjectTypes
 */
ENUM_STRINGS( MeshObjectPath::ObjectTypes,
              MeshLevel::groupStructKeys::nodeManagerString(),
              MeshLevel::groupStructKeys::edgeManagerString(),
              MeshLevel::groupStructKeys::faceManagerString(),
              MeshLevel::groupStructKeys::elemManagerString(),
              "invalid" );

template< typename OBJECT_TYPE >
bool MeshObjectPath::checkObjectTypeConsistency() const
{
  bool consistent = false;
  if( m_objectType == ObjectTypes::nodes )
  {
    consistent = std::is_same< NodeManager, OBJECT_TYPE >::value;
  }
  else if( m_objectType == ObjectTypes::edges )
  {
    consistent = std::is_same< EdgeManager, OBJECT_TYPE >::value;
  }
  else if( m_objectType == ObjectTypes::faces )
  {
    consistent = std::is_same< FaceManager, OBJECT_TYPE >::value;
  }
  else if( m_objectType == ObjectTypes::elems )
  {
    consistent = std::is_base_of< ElementRegionBase, OBJECT_TYPE >::value ||
                 std::is_base_of< ElementSubRegionBase, OBJECT_TYPE >::value;
  }
  return consistent;
}

template< typename OBJECT_TYPE,
          typename FUNC >
void MeshObjectPath::forObjectsInPath( std::pair< string const, std::map< string, std::vector< string > > > const & levelPair,
                                       MeshLevel & meshLevel,
                                       FUNC && func ) const
{
  forObjectsInPath< OBJECT_TYPE >( levelPair, const_cast< MeshLevel const & >( meshLevel ), [&]( OBJECT_TYPE const & object )
  {
    func( const_cast< OBJECT_TYPE & >(object) );
  } );
}

template< typename OBJECT_TYPE,
          typename FUNC >
void MeshObjectPath::forObjectsInPath( std::pair< string const, std::map< string, std::vector< string > > > const & levelPair,
                                       MeshLevel const & meshLevel,
                                       FUNC && func ) const
{
  if( m_objectType == ObjectTypes::nodes )
  {
    func( dynamic_cast< OBJECT_TYPE const & >(meshLevel.getNodeManager() ) );
  }
  else if( m_objectType == ObjectTypes::edges )
  {
    func( dynamic_cast< OBJECT_TYPE const & >(meshLevel.getEdgeManager()) );
  }
  else if( m_objectType == ObjectTypes::faces )
  {
    func( dynamic_cast< OBJECT_TYPE const & >(meshLevel.getFaceManager()) );
  }
  else if( m_objectType == ObjectTypes::elems )
  {
    ElementRegionManager const & elemRegionMan = meshLevel.getElemManager();
    for( auto & elemRegionPair : levelPair.second )
    {
      ElementRegionBase const & elemRegion = elemRegionMan.getRegion( elemRegionPair.first );
      if( std::is_base_of< ElementRegionBase, OBJECT_TYPE >::value )
      {
        func( dynamic_cast< OBJECT_TYPE const & >(elemRegion) );
      }
      else
      {
        for( auto & elemSubRegionName : elemRegionPair.second )
        {
          ElementSubRegionBase const & subRegion = elemRegion.getSubRegion( elemSubRegionName );
          if( std::is_base_of< ElementSubRegionBase, OBJECT_TYPE >::value ||
              std::is_same< dataRepository::Group, OBJECT_TYPE >::value )
          {
            func( dynamic_cast< OBJECT_TYPE const & >(subRegion) );
          }
          else
          {
            GEOS_ERROR( "You shouldn't be here" );
          }
        }
      }
    }
  }
}

template< typename OBJECT_TYPE,
          typename FUNC >
void MeshObjectPath::forObjectsInPath( dataRepository::Group & meshBodies,
                                       FUNC && func ) const
{
  forObjectsInPath< OBJECT_TYPE >( const_cast< dataRepository::Group const & >(meshBodies),
                                   [&]( auto const & object )
  {
    func( const_cast< OBJECT_TYPE & >(object) );
  } );
}

template< typename OBJECT_TYPE, typename FUNC >
void MeshObjectPath::forObjectsInPath( dataRepository::Group const & meshBodies,
                                       FUNC && func ) const
{
  checkObjectTypeConsistency< OBJECT_TYPE >();
  for( auto const & meshBodyPair : m_pathPermutations )
  {
    MeshBody const & meshBody = meshBodies.getGroup< MeshBody >( meshBodyPair.first );
    for( auto const & meshLevelPair : meshBodyPair.second )
    {
      MeshLevel const & meshLevel = meshBody.getMeshLevel( meshLevelPair.first );
      forObjectsInPath< OBJECT_TYPE, FUNC >( meshLevelPair, meshLevel, std::forward< FUNC >( func ));
    }
  }
}

template< typename OBJECT_TYPE,
          typename FUNC >
void MeshObjectPath::forObjectsInPath( MeshLevel & meshLevel,
                                       FUNC && func ) const
{
  string const bodyName = meshLevel.getParent().getParent().getName();
  string const levelName = meshLevel.getName();

  auto bodyIter = m_pathPermutations.find( bodyName );
  if( bodyIter != m_pathPermutations.end() )
  {
    auto const levelIter = bodyIter->second.find( levelName );
    if( levelIter != bodyIter->second.end() )
    {
      // string const objectTypeName = stringutilities::toLower( OBJECT_TYPE::catalogName());
      // if( objectTypeName == stringutilities::toLower( EnumStrings< ObjectTypes >::toString( m_objectType ) ) ||
      //     ( objectTypeName.find( "elem" )!=string::npos && m_objectType==ObjectTypes::elems ) ||
      //     ( objectTypeName == "group") )
      if( checkObjectTypeConsistency< OBJECT_TYPE >() ||
          std::is_same< OBJECT_TYPE, dataRepository::Group >::value )
      {
        forObjectsInPath< OBJECT_TYPE, FUNC >( *levelIter, meshLevel, std::forward< FUNC >( func ) );
      }
    }
  }
}



} /* namespace geos */

#endif /* GEOS_MESH_MESHOBJECTPATH_HPP_ */
