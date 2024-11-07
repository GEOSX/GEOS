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
 * @file Group.hpp
 */

#ifndef GEOS_DATAREPOSITORY_GROUP_HPP_
#define GEOS_DATAREPOSITORY_GROUP_HPP_

#include "InputFlags.hpp"
#include "ObjectCatalog.hpp"
#include "MappedVector.hpp"
#include "RestartFlags.hpp"
#include "Wrapper.hpp"
#include "xmlWrapper.hpp"
#include "LogLevelsInfo.hpp"
#include "LogLevelsRegistry.hpp"


#include <iostream>

#ifndef NOCHARTOSTRING_KEYLOOKUP
/// macro definition to enable/disable char * lookups
#define NOCHARTOSTRING_KEYLOOKUP 0
#endif

/**
 * namespace to encapsulate GEOSX
 */
namespace geos
{

/**
 * Encapsulates all dataRepository classes and functionality.
 */
namespace dataRepository
{

//START_SPHINX_INCLUDE_00
/// The default key type for entries in the hierarchy.
using keyType = string;

/// The default index type for entries the hierarchy.
using indexType = localIndex;
//END_SPHINX_INCLUDE_00

/**
 * @class Group
 *
 * The Group class serves as a "node" in a hierarchy of the dataRepository. The data structure is built as a
 * hierarchy of Group objects, or objects derived from group objects.
 */
class Group
{
public:
  //START_SPHINX_INCLUDE_01
  /// The template specialization of MappedVector to use for the collection of sub-Group objects.
  using subGroupMap = MappedVector< Group, Group *, keyType, indexType >;

  /// The template specialization of MappedVector to use for the collection wrappers objects.
  using wrapperMap = MappedVector< WrapperBase, WrapperBase *, keyType, indexType >;
  //END_SPHINX_INCLUDE_01

  /**
   * @name Constructors/destructor
   */
  ///@{

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  explicit Group( string const & name,
                  Group * const parent );

  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param rootNode The root node of the data repository.
   * @note This Group will not have a parent group.
   */
  explicit Group( string const & name,
                  conduit::Node & rootNode );

  /**
   * @brief Move constructor
   * @param[in] source source Group
   */
  Group( Group && source ) = default;

  /**
   * @brief Destructor, deletes all Groups and Wrappers owned by this Group
   */
  virtual ~Group();

  /**
   * @brief Deleted default constructor.
   */
  Group() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  Group( Group const & ) = delete;


  /**
   * @brief Deleted copy assignment operator.
   * @return
   */
  Group & operator=( Group const & ) = delete;

  /**
   * @brief Deleted move assignment operator.
   * @return
   */
  Group & operator=( Group && ) = delete;

  ///@}


  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Type alias for catalog interface used by this class. See CatalogInterface.
   */
  using CatalogInterface = dataRepository::CatalogInterface< Group, string const &, Group * const >;

  /**
   * @brief Get the singleton catalog for this Group.
   * @return reference to the catalog object
   */
  static CatalogInterface::CatalogType & getCatalog();

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  /**
   * @brief Prints the data hierarchy recursively.
   * @param[in] indent The level of indentation to add to this level of output.
   */
  void printDataHierarchy( integer indent = 0 ) const;

  /**
   * @brief @return a table formatted string containing all input options.
   */
  string dumpInputOptions() const;

  /**
   * @brief @return a comma separated string containing all sub groups name.
   */
  string dumpSubGroupsNames() const;

  /**
   * @brief @return a comma separated string containing all wrappers name.
   */
  string dumpWrappersNames() const;

  ///@}

  //START_SPHINX_INCLUDE_REGISTER_GROUP
  /**
   * @name Sub-group registration interface
   */
  ///@{

  /**
   * @brief Register a new Group as a sub-group of current Group.
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param[in] name      The name of the group to use as a string key.
   * @param[in] newObject A unique_ptr to the object that is being registered.
   * @return              A pointer to the newly registered Group.
   *
   * Registers a Group or class derived from Group as a subgroup of this Group and takes ownership.
   */
  template< typename T = Group >
  T & registerGroup( string const & name, std::unique_ptr< T > newObject )
  {
    newObject->m_parent = this;
    return dynamicCast< T & >( *m_subGroups.insert( name, newObject.release(), true ) );
  }

  /**
   * @brief @copybrief registerGroup(string const &,std::unique_ptr<T>)
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param[in] name          The name of the group to use as a string key.
   * @param[in] newObject     A unique_ptr to the object that is being registered.
   * @return                  A pointer to the newly registered Group.
   *
   * Registers a Group or class derived from Group as a subgroup of this Group but does not take ownership.
   */
  template< typename T = Group >
  T & registerGroup( string const & name, T * newObject )
  { return dynamicCast< T & >( *m_subGroups.insert( name, newObject, false ) ); }


  /**
   * @brief @copybrief registerGroup(string const &,std::unique_ptr<T>)
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param[in] name The name of the group to use as a string key.
   * @return         A pointer to the newly registered Group.
   *
   * Creates and registers a Group or class derived from Group as a subgroup of this Group.
   */
  template< typename T = Group >
  T & registerGroup( string const & name )
  { return registerGroup< T >( name, std::make_unique< T >( name, this ) ); }

  /**
   * @brief @copybrief registerGroup(string const &,std::unique_ptr<T>)
   *
   * @tparam T The type of the Group to add/register. This should be a type that derives from Group.
   * @param keyIndex A KeyIndexT object that will be used to specify the name of
   *   the new group. The index of the KeyIndex will also be set.
   * @return A pointer to the newly registered Group, or @c nullptr if no group was registered.
   *
   * Creates and registers a Group or class derived from Group as a subgroup of this Group.
   */
  template< typename T = Group >
  T & registerGroup( subGroupMap::KeyIndex const & keyIndex )
  {
    T & rval = registerGroup< T >( keyIndex.key(), std::make_unique< T >( keyIndex.key(), this ) );
    keyIndex.setIndex( m_subGroups.getIndex( keyIndex.key() ) );
    return rval;
  }

  /**
   * @brief Removes a child group from this group.
   * @param name the name of the child group to remove from this group.
   */
  void deregisterGroup( string const & name );

  /**
   * @brief Creates a new sub-Group using the ObjectCatalog functionality.
   * @param[in] childKey The name of the new object type's key in the
   *                     ObjectCatalog.
   * @param[in] childName The name of the new object in the collection of
   *                      sub-Groups.
   * @return A pointer to the new Group created by this function.
   */
  virtual Group * createChild( string const & childKey, string const & childName );

  ///@}

  //END_SPHINX_INCLUDE_REGISTER_GROUP

  //START_SPHINX_INCLUDE_GET_GROUP
  /**
   * @name Sub-group retrieval methods.
   *
   * This collection of functions are used to get a sub-Group from the current group. Various methods
   * for performing the lookup are provided (localIndex, string, KeyIndex), and each have their
   * advantages and costs. The lowest cost lookup is the "localIndex" lookup. The KeyIndex lookup
   * will add a cost for checking to make sure the index stored in KeyIndex is valid (a string
   * compare, and a hash if it is incorrect). The string lookup is the full cost hash lookup every
   * time that it is called.
   *
   * The template parameter specifies the "type" that the caller expects to lookup, and thus attempts
   * to cast the pointer that is stored in m_subGroups to a pointer of the desired type. If this
   * cast fails, then a @p nullptr is returned. If no template parameter is specified then a default
   * type of Group is assumed.
   */
  ///@{

  /**
   * @brief Return a pointer to a sub-group of the current Group.
   * @tparam T The type of subgroup.
   * @tparam KEY The type of the lookup.
   * @param key The key used to perform the lookup.
   * @return A pointer to @p T that refers to the sub-group, if the Group does not exist or it
   *   has an incompatible type a @c nullptr is returned.
   */
  template< typename T = Group, typename KEY = void >
  T * getGroupPointer( KEY const & key )
  { return dynamicCast< T * >( m_subGroups[ key ] ); }

  /**
   * @copydoc getGroupPointer(KEY const &)
   */
  template< typename T = Group, typename KEY = void >
  T const * getGroupPointer( KEY const & key ) const
  { return dynamicCast< T const * >( m_subGroups[ key ] ); }

  /**
   * @brief Return a reference to a sub-group of the current Group.
   * @tparam T The type of subgroup.
   * @tparam KEY The type of the lookup.
   * @param key The key used to perform the lookup.
   * @return A reference to @p T that refers to the sub-group.
   * @throw std::domain_error If the Group does not exist is thrown.
   */
  template< typename T = Group, typename KEY = void >
  T & getGroup( KEY const & key )
  {
    Group * const child = m_subGroups[ key ];
    GEOS_THROW_IF( child == nullptr,
                   "Group " << getDataContext() << " has no child named " << key << std::endl
                            << dumpSubGroupsNames(),
                   std::domain_error );
    T * const castedChild = dynamicCast< T * >( child );
    GEOS_THROW_IF( castedChild == nullptr,
                   GEOS_FMT( "{} was expected to be a '{}'.",
                             child->getDataContext(), LvArray::system::demangleType< T >() ),
                   BadTypeError );
    return *castedChild;
  }

  /**
   * @copydoc getGroup( KEY const & )
   */
  template< typename T = Group, typename KEY = void >
  T const & getGroup( KEY const & key ) const
  {
    Group const * const child = m_subGroups[ key ];
    GEOS_THROW_IF( child == nullptr,
                   "Group " << getDataContext() << " has no child named " << key << std::endl
                            << dumpSubGroupsNames(),
                   std::domain_error );
    T const * const castedChild = dynamicCast< T const * >( child );
    GEOS_THROW_IF( castedChild == nullptr,
                   GEOS_FMT( "{} was expected to be a '{}'.",
                             child->getDataContext(), LvArray::system::demangleType< T >() ),
                   BadTypeError );
    return *castedChild;
  }

  /**
   * @brief Retrieve a group from the hierarchy using a path.
   * @tparam T type of subgroup
   * @param[in] path a unix-style string (absolute, relative paths valid)
   *                 to lookup the Group to return. Absolute paths search
   *                 from the tree root, while relative - from current group.
   * @return A reference to @p T that refers to the sub-group.
   * @throw std::domain_error If the Group doesn't exist.
   */
  template< typename T = Group >
  T & getGroupByPath( string const & path )
  { return dynamicCast< T & >( const_cast< Group & >( getBaseGroupByPath( path ) ) ); }

  /**
   * @copydoc getGroupByPath(string const &)
   */
  template< typename T = Group >
  T const & getGroupByPath( string const & path ) const
  { return dynamicCast< T const & >( getBaseGroupByPath( path ) ); }

  //END_SPHINX_INCLUDE_GET_GROUP

  /**
   * @brief Get the subgroups object
   * @return a reference to the sub-group map.
   */
  subGroupMap & getSubGroups()
  { return m_subGroups; }

  /**
   * @brief Get the subgroups object
   * @return a reference to const that points to the sub-group map.
   */
  subGroupMap const & getSubGroups() const
  { return m_subGroups; }

  /**
   * @brief return the number of sub groups in this Group
   * @return number of sub groups in this Group
   */
  localIndex numSubGroups() const { return m_subGroups.size(); }

  /**
   * @return An array containing all sub groups keys
   */
  std::vector< string > getSubGroupsNames() const;

  /**
   * @brief Check whether a sub-group exists.
   * @param name the name of sub-group to search for
   * @return @p true if sub-group exists, @p false otherwise
   */
  template< typename T = Group >
  bool hasGroup( string const & name ) const
  { return dynamicCast< T const * >( m_subGroups[ name ] ) != nullptr; }

  /**
   * @brief Check whether a sub-group exists by type.
   * @tparam T The type of sub-group to search for
   * @return @p true if sub-group of type T exists, @p false otherwise
   */
  template< typename T >
  bool hasSubGroupOfType( ) const
  {
    bool hasSubGroup = false;
    // since forSubGroups only applies the lambda to groups matching the type,
    //   any calls to the lambda indicates that we have a subgroup of the correct type.
    forSubGroups< T >( [&]( T const & ){ hasSubGroup = true; } );
    return hasSubGroup;
  }

  ///@}

  /**
   * @name Functor application helpers
   *
   * This function is useful when trying to apply a functor that passes a pointer to an container,
   * but it is desired that the functor is only executed if the container can be casted to a certain
   * type. The variadic list consisting of CASTTYPE/S will be used recursively to check if the
   * container is able to be casted to the one of these types. The first type in the CASTTYPE/S list
   * will be used to execute the functor, and the function will return true.
   */
  ///@{

  /** @cond DO_NOT_DOCUMENT */
  template< typename CASTTYPE, typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE & container, LAMBDA && lambda )
  {
    using T = std::conditional_t< std::is_const< CONTAINERTYPE >::value, CASTTYPE const, CASTTYPE >;
    T * const castedContainer = dynamic_cast< T * >( &container );
    if( castedContainer != nullptr )
    {
      lambda( *castedContainer );
      return true;
    }

    return false;
  }
  /** @endcond */

  /**
   * @brief Apply a given functor to a container if the container can be
   *        cast to one of the specified types.
   * @tparam CASTTYPE      the first type that will be used in the attempted casting of container
   * @tparam CASTTYPES     a variadic list of types that will be used in the attempted casting of container
   * @tparam CONTAINERTYPE the type of container
   * @tparam LAMBDA        the type of lambda function to call in the function
   * @param[in] container  a pointer to the container which will be passed to the lambda function
   * @param[in] lambda     the lambda function to call in the function
   * @return               a boolean to indicate whether the lambda was successfully applied to the container.
   */
  template< typename T0, typename T1, typename ... CASTTYPES, typename CONTAINERTYPE, typename LAMBDA >
  static bool applyLambdaToContainer( CONTAINERTYPE & container, LAMBDA && lambda )
  {
    using T = std::conditional_t< std::is_const< CONTAINERTYPE >::value, T0 const, T0 >;
    T * const castedContainer = dynamic_cast< T * >( &container );

    if( castedContainer != nullptr )
    {
      lambda( *castedContainer );
      return true;
    }

    return applyLambdaToContainer< T1, CASTTYPES... >( container, std::forward< LAMBDA >( lambda ) );
  }
  ///@}


  //START_SPHINX_INCLUDE_LOOP_INTERFACE
  /**
   * @name Functor-based subgroup iteration
   *
   * These functions loop over sub-groups and executes a functor that uses the sub-group as an
   * argument. The functor is only executed if the group can be cast to a certain type specified
   * by the @p ROUPTYPE/S pack. The variadic list consisting of @p GROUPTYPE/S will be used recursively
   * to check if the group is able to be cast to the one of these types. The first type in the
   * @p GROUPTYPE/S list will be used to execute the functor, and the next sub-group will be processed.
   */
  ///@{

  /**
   * @brief Apply the given functor to subgroups that can be casted to one of specified types.
   * @tparam GROUPTYPE  the first type that will be used in the attempted casting of group.
   * @tparam GROUPTYPES a variadic list of types that will be used in the attempted casting of group.
   * @tparam LAMBDA     the type of functor to call
   * @param[in] lambda  the functor to call on subgroups
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forSubGroups( LAMBDA && lambda )
  {
    for( auto & subGroupIter : m_subGroups )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( *subGroupIter.second, [&]( auto & castedSubGroup )
      {
        lambda( castedSubGroup );
      } );
    }
  }

  /**
   * @copydoc forSubGroups(LAMBDA &&)
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forSubGroups( LAMBDA && lambda ) const
  {
    for( auto const & subGroupIter : m_subGroups )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( *subGroupIter.second, [&]( auto const & castedSubGroup )
      {
        lambda( castedSubGroup );
      } );
    }
  }


  /**
   * @brief Apply the given functor to subgroups that can be casted to one of specified types.
   * @tparam GROUPTYPE  the first type that will be used in the attempted casting of group.
   * @tparam GROUPTYPES a variadic list of types that will be used in the attempted casting of group.
   * @tparam LAMBDA     the type of functor to call
   * @param[in] lambda  the functor to call on subgroups
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forSubGroupsIndex( LAMBDA && lambda )
  {
    localIndex counter = 0;
    for( auto & subGroupIter : m_subGroups )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( *subGroupIter.second,
                                                          [&]( auto & castedSubGroup )
      {
        lambda( counter, castedSubGroup );
      } );
      ++counter;
    }
  }

  /**
   * @copydoc forSubGroupsIndex(LAMBDA &&)
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LAMBDA >
  void forSubGroupsIndex( LAMBDA && lambda ) const
  {
    localIndex counter = 0;
    for( auto const & subGroupIter : m_subGroups )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( *subGroupIter.second,
                                                          [&]( auto const & castedSubGroup )
      {
        lambda( counter, castedSubGroup );
      } );
      ++counter;
    }
  }

  /**
   * @copybrief forSubGroups(LAMBDA &&)
   * @tparam GROUPTYPE        the first type that will be used in the attempted casting of group.
   * @tparam GROUPTYPES       a variadic list of types that will be used in the attempted casting of group.
   * @tparam LOOKUP_CONTAINER type of container of subgroup lookup keys (names or indices), must support range-based for
   * loop
   * @tparam LAMBDA           type of functor callable with an index in lookup container and a reference to casted
   * subgroup
   * @param[in] subGroupKeys  container with subgroup lookup keys (e.g. names or indices) to apply the functor to
   * @param[in] lambda        the functor to call
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forSubGroups( LOOKUP_CONTAINER const & subGroupKeys, LAMBDA && lambda )
  {
    localIndex counter = 0;

    for( auto const & subgroup : subGroupKeys )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( getGroup( subgroup ), [&]( auto & castedSubGroup )
      {
        lambda( counter, castedSubGroup );
      } );
      ++counter;
    }
  }

  /**
   * @copybrief forSubGroups(LAMBDA &&)
   * @tparam GROUPTYPE        the first type that will be used in the attempted casting of group.
   * @tparam GROUPTYPES       a variadic list of types that will be used in the attempted casting of group.
   * @tparam LOOKUP_CONTAINER type of container of subgroup lookup keys (names or indices), must support range-based for
   * loop
   * @tparam LAMBDA           type of functor callable with an index in lookup container and a reference to casted
   * subgroup
   * @param[in] subGroupKeys  container with subgroup lookup keys (e.g. names or indices) to apply the functor to
   * @param[in] lambda        the functor to call
   */
  template< typename GROUPTYPE = Group, typename ... GROUPTYPES, typename LOOKUP_CONTAINER, typename LAMBDA >
  void forSubGroups( LOOKUP_CONTAINER const & subGroupKeys, LAMBDA && lambda ) const
  {
    localIndex counter = 0;
    for( auto const & subgroup : subGroupKeys )
    {
      applyLambdaToContainer< GROUPTYPE, GROUPTYPES... >( getGroup( subgroup ), [&]( auto const & castedSubGroup )
      {
        lambda( counter, castedSubGroup );
      } );
      ++counter;
    }
  }
  ///@}

  /**
   * @name Functor-based wrapper iteration
   *
   * These functions loop over the wrappers contained in this group, and executes a functor that
   * uses the Wrapper as an argument. The functor is only executed if the Wrapper can be casted to
   * a certain type specified by the @p TYPE/S pack. The variadic list consisting of
   * @p TYPE/S will be used recursively to check if the Wrapper is able to be casted to the
   * one of these types. The first type in the @p WRAPPERTYPE/S list will be used to execute the
   * functor, and the next Wrapper will be processed.
   */
  ///@{

  /**
   * @brief Apply the given functor to wrappers.
   * @tparam LAMBDA the type of functor to call
   * @param[in] lambda  the functor to call
   */
  template< typename LAMBDA >
  void forWrappers( LAMBDA && lambda )
  {
    for( auto & wrapperIter : m_wrappers )
    {
      lambda( *wrapperIter.second );
    }
  }

  /**
   * @copydoc forWrappers(LAMBDA &&)
   */
  template< typename LAMBDA >
  void forWrappers( LAMBDA && lambda ) const
  {
    for( auto const & wrapperIter : m_wrappers )
    {
      lambda( *wrapperIter.second );
    }
  }

  /**
   * @brief Apply the given functor to wrappers that can be cast to one of specified types.
   * @tparam TYPE   the first type that will be used in the attempted casting of Wrapper
   * @tparam TYPES  a variadic list of types that will be used in the attempted casting of Wrapper
   * @tparam LAMBDA the type of functor to call
   * @param[in] lambda  the functor to call
   */
  template< typename TYPE, typename ... TYPES, typename LAMBDA >
  void forWrappers( LAMBDA && lambda )
  {
    for( auto & wrapperIter : m_wrappers )
    {
      applyLambdaToContainer< Wrapper< TYPE >, Wrapper< TYPES >... >( *wrapperIter.second,
                                                                      std::forward< LAMBDA >( lambda ));
    }
  }

  /**
   * @brief Apply the given functor to wrappers that can be cast to one of specified types.
   * @tparam TYPE   the first type that will be used in the attempted casting of Wrapper
   * @tparam TYPES  a variadic list of types that will be used in the attempted casting of Wrapper
   * @tparam LAMBDA the type of functor to call
   * @param[in] lambda  the functor to call
   */
  template< typename TYPE, typename ... TYPES, typename LAMBDA >
  void forWrappers( LAMBDA && lambda ) const
  {
    for( auto const & wrapperIter : m_wrappers )
    {
      applyLambdaToContainer< Wrapper< TYPE >, Wrapper< TYPES >... >( *wrapperIter.second,
                                                                      std::forward< LAMBDA >( lambda ));
    }
  }

  ///@}
  //END_SPHINX_INCLUDE_LOOP_INTERFACE

  /**
   * @name Initialization and data registration
   */
  ///@{
  /**
   * @brief initialization post generation of the mesh.
   */
  virtual void initialize_postMeshGeneration();


  /**
   * @brief Run initialization functions on this and all subgroups.
   * @param[in] root A group that is passed in to the initialization functions
   *   in order to facilitate the initialization.
   *
   * This function will first call initializePreSubGroups() on this Group, then
   * loop over all subgroups and call Initialize() on them, then
   * call InitializePostSubGroups() on this Group.
   *
   * @note The order in which the sub-Groups are iterated over is defined by
   * InitializationOrder().
   */
  void initialize();

  /**
   * @brief Sets the initialization order for sub-Groups.
   * @param[out] order An array of strings that define the iteration order.
   *
   * This function will fill the @p order array that is used to specify the
   * order in which the Initialize() function loops over sub-Groups. If a
   * custom order is required by a derived type, this function should be
   * overridden with a implementation that specifies the desired order.
   */
  virtual void initializationOrder( string_array & order );


  /**
   * @brief Initialization routine to be called after calling
   *        ApplyInitialConditions().
   *
   * This function provides a capability for post-initial condition problem
   * initialization. First the InitializePostInitialConditions_PreSubGroups()
   * function is called on this Group. Then there is a loop over all subgroups
   * and InitializePostInitialConditions() is called on them. Finally, the
   * InitializePostInitialConditions_PostSubGroups() function is called this
   * Group.
   *
   * @note The order in which the sub-Groups are iterated over is defined by
   * InitializationOrder().
   */
  void initializePostInitialConditions();

  /**
   * @brief Initialization routine to be called after calling reading a restart file.
   *
   * This functions recurses and calls postRestartInitialization() on nested sub-groups
   * at any depth, providing a capability to add custom post-restart initialization.
   */
  void postRestartInitializationRecursive();

  static string processInputName( xmlWrapper::xmlNode const & targetNode,
                                  xmlWrapper::xmlNodePos const & targetNodePos,
                                  string_view parentNodeName,
                                  xmlWrapper::xmlNodePos const & parentNodePos,
                                  std::set< string > & siblingNames );

  /**
   * @brief Recursively read values using ProcessInputFile() from the input
   * file and put them into the wrapped values for this group.
   * Also add the includes content to the xmlDocument when `Include` nodes are encountered.
   * @param[in] xmlDocument the XML document that contains the targetNode.
   * @param[in] targetNode the XML node that to extract input values from.
   */
  void processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                  xmlWrapper::xmlNode & targetNode );
  /**
   * @brief Same as processInputFileRecursive(xmlWrapper::xmlDocument &, xmlWrapper::xmlNode &)
   * but allow to reuse an existing xmlNodePos.
   * @param[in] xmlDocument the XML document that contains the targetNode.
   * @param[in] targetNode the XML node that to extract input values from.
   * @param[in] nodePos the target node position, typically obtained with xmlDocument::getNodePosition().
   */
  void processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                  xmlWrapper::xmlNode & targetNode,
                                  xmlWrapper::xmlNodePos const & targetNodePos );

  /**
   * @brief Recursively call postInputInitialization() to apply post processing after
   * reading input values.
   */
  void postInputInitializationRecursive();

  ///@}

  //START_SPHINX_INCLUDE_REGISTER_WRAPPER
  /**
   * @name Wrapper registration interface
   */
  ///@{

  /**
   * @brief Create and register a Wrapper around a new object.
   * @tparam T The type of the object allocated.
   * @tparam TBASE The type of the object that the Wrapper holds.
   * @param[in] name the name of the wrapper to use as a string key
   * @param[out] rkey a pointer to a index type that will be filled with the new
   *   Wrapper index in this Group
   * @return A reference to the newly registered/created Wrapper
   */
  template< typename T, typename TBASE=T >
  Wrapper< TBASE > & registerWrapper( string const & name,
                                      wrapperMap::KeyIndex::index_type * const rkey = nullptr );

  /**
   * @copybrief registerWrapper(string const &,wrapperMap::KeyIndex::index_type * const)
   * @tparam T the type of the wrapped object
   * @tparam TBASE the base type to cast the returned wrapper to
   * @param[in] viewKey The KeyIndex that contains the name of the new Wrapper.
   * @return A reference to the newly registered/created Wrapper
   */
  template< typename T, typename TBASE=T >
  Wrapper< TBASE > & registerWrapper( Group::wrapperMap::KeyIndex const & viewKey );

  /**
   * @brief Register a Wrapper around a given object and take ownership.
   * @tparam T the type of the wrapped object
   * @param[in] name the name of the wrapper to use as a string key
   * @param[in] newObject an owning pointer to the object that is being registered
   * @return A reference to the newly registered/created Wrapper
   * @note Not intended to register a @p WrapperBase instance. Use dedicated member function instead.
   */
  template< typename T >
  Wrapper< T > & registerWrapper( string const & name, std::unique_ptr< T > newObject );

  /**
   * @brief Register a Wrapper around an existing object, does not take ownership of the object.
   * @tparam T the type of the wrapped object
   * @param[in] name the name of the wrapper to use as a string key
   * @param[in] newObject a pointer to the object that is being registered
   * @return A reference to the newly registered/created Wrapper
   * @note Not intended to register a @p WrapperBase instance. Use dedicated member function instead.
   */
  template< typename T >
  Wrapper< T > & registerWrapper( string const & name,
                                  T * newObject );

  /**
   * @brief Register and take ownership of an existing Wrapper.
   * @param wrapper A pointer to the an existing wrapper.
   * @return An un-typed pointer to the newly registered/created wrapper
   */
  WrapperBase & registerWrapper( std::unique_ptr< WrapperBase > wrapper );

  /**
   * @brief Removes a Wrapper from this group.
   * @param name the name of the Wrapper to remove from this group.
   */
  void deregisterWrapper( string const & name );

  ///@}
  //END_SPHINX_INCLUDE_REGISTER_WRAPPER

  /**
   * @brief Append a levelCondition and a log description to the description of the wrapped object given a log info struct.
   * Must be called in constructor.
   * @tparam LOG_LEVEL_INFO The log documentation to add.
   * @return void if the trait is verified.
   */
  template< typename LOG_LEVEL_INFO >
  std::enable_if_t< geos::is_log_level_info< LOG_LEVEL_INFO >, void >
  addLogLevel();

  /**
   * @name Schema generation methods
   */
  ///@{

  /**
   * @brief Build a complete datastructure for schema generation.
   * @param level indent level for printing out the structure
   */
  void generateDataStructureSkeleton( integer const level )
  {
    expandObjectCatalogs();
    string indent( level*2, ' ' );

    for( auto const & subGroupIter : m_subGroups )
    {
      std::cout << indent << subGroupIter.second->getName() << std::endl;
      subGroupIter.second->generateDataStructureSkeleton( level + 1 );
    }
  }

  /**
   * @brief Expand any catalogs in the data structure.
   */
  virtual void expandObjectCatalogs() {}

  /**
   * @brief Inform the schema generator of any deviations between the xml and GEOS data structures.
   * @param schemaRoot        XML node corresponding to the root
   * @param schemaParent      XML node for the parent node
   * @param documentationType type of XML schema generated
   */
  virtual void setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                    xmlWrapper::xmlNode schemaParent,
                                    integer documentationType )
  {
    GEOS_UNUSED_VAR( schemaRoot );
    GEOS_UNUSED_VAR( schemaParent );
    GEOS_UNUSED_VAR( documentationType );
  }

  ///@}

  /**
   * @name Mesh data registration
   */
  ///@{

  /**
   * @brief Calls RegisterDataOnMesh() recursively.
   * @param[in,out] meshBodies the group of MeshBody objects to register data on.
   */
  virtual void registerDataOnMeshRecursive( Group & meshBodies );

  /**
   * @brief Register data on mesh entities.
   * @param[in,out] meshBodies the group of MeshBody objects to register data on.
   *
   * This function is used to register data on mesh entities such as the NodeManager,
   * FaceManager...etc.
   */
  virtual void registerDataOnMesh( Group & meshBodies )
  {
    GEOS_UNUSED_VAR( meshBodies );
  }

  ///@}

  /**
   * @name Packing/upacking methods
   */
  ///@{

  /**
   * @brief Get the size required to pack a list of wrappers.
   * @param[in] wrapperNames an array that contains the names of the wrappers to pack.
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] onDevice    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @param[out] events      a collection of events to poll for completion of async
   *                         packing kernels ( device packing is incomplete until all
   *                         events are finalized )
   * @return                 the size of the buffer required to pack the wrappers.
   */
  virtual localIndex packSize( string_array const & wrapperNames,
                               integer const recursive,
                               bool onDevice,
                               parallelDeviceEvents & events ) const;

  /**
   * @brief Get the size required to pack a list of indices within a list of wrappers.
   * @param[in] wrapperNames an array that contains the names of the wrappers to pack.
   * @param[in] packList     the list of indices to pack
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] onDevice     whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @param[out] events      a collection of events to poll for completion of async
   *                         packing kernels ( device packing is incomplete until all
   *                         events are finalized )
   * @return                 the size of the buffer required to pack the wrapper indices.
   */
  virtual localIndex packSize( string_array const & wrapperNames,
                               arrayView1d< localIndex const > const & packList,
                               integer const recursive,
                               bool onDevice,
                               parallelDeviceEvents & events ) const;

  /**
   * @brief Get the size required to pack a list of indices for all registered wrappers.
   * @param[in] packList     the list of indices to pack
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] onDevice     whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @param[out] events      a collection of events to poll for completion of async
   *                         packing kernels ( device packing is incomplete until all
   *                         events are finalized )
   * @return                 the size of the buffer required to pack the wrapper indices.
   */
  localIndex packSize( arrayView1d< localIndex const > const & packList,
                       integer const recursive,
                       bool onDevice,
                       parallelDeviceEvents & events ) const;

  /**
   * @brief Pack a list of wrappers to a buffer.
   * @param[in,out] buffer   the buffer that will be packed.
   * @param[in] wrapperNames an array that contains the names of the wrappers to pack.
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] onDevice    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @param[out] events      a collection of events to poll for completion of async
   *                         packing kernels ( device packing is incomplete until all
   *                         events are finalized )
   * @return                 the size of data packed to the buffer.
   *
   * This function takes in a reference to a pointer @p buffer, and packs data specified by
   * @p wrapperNames, and @p recursive to that pointer location. The
   * pointer is altered and returned to the new location corresponding the
   * original value of @p buffer plus the size of data packed to the buffer.
   *
   */
  virtual localIndex pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           integer const recursive,
                           bool onDevice,
                           parallelDeviceEvents & events ) const;

  /**
   * @brief Pack a list of indices within a list of wrappers.
   * @param[in,out] buffer   the buffer that will be packed.
   * @param[in] wrapperNames an array that contains the names of the wrappers to pack.
   * @param[in] packList     the list of indices to pack
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] onDevice     whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @param[out] events      a collection of events to poll for completion of async
   *                         packing kernels ( device packing is incomplete until all
   *                         events are finalized )
   * @return                 the size of data packed to the buffer.
   *
   * This function takes in a reference to a pointer @p buffer, and packs data specified by
   * @p wrapperNames, @p packList, and @p recursive to that pointer location. The
   * pointer is altered and returned to the new location corresponding the
   * original value of @p buffer plus the size of data packed to the buffer.
   */
  virtual localIndex pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           arrayView1d< localIndex const > const & packList,
                           integer const recursive,
                           bool onDevice,
                           parallelDeviceEvents & events ) const;

  /**
   * @brief Pack a list of indices for all registered wrappers.
   * @param[in,out] buffer   the buffer that will be packed.
   * @param[in] packList     the list of indices to pack
   * @param[in] recursive    whether or not to perform a recursive pack.
   * @param[in] onDevice     whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @param[out] events      a collection of events to poll for completion of async
   *                         packing kernels ( device packing is incomplete until all
   *                         events are finalized )
   * @return                 the size of data packed to the buffer.
   *
   * This function takes in a reference to a pointer @p buffer, and packs data specified by
   * @p wrapperNames, @p packList, and @p recursive to that pointer location. The
   * pointer is altered and returned to the new location corresponding the
   * original value of @p buffer plus the size of data packed to the buffer.
   */
  localIndex pack( buffer_unit_type * & buffer,
                   arrayView1d< localIndex const > const & packList,
                   integer const recursive,
                   bool onDevice,
                   parallelDeviceEvents & events ) const;

  /**
   * @brief Unpack a buffer.
   * @param[in,out] buffer   the buffer to unpack
   * @param[in,out] packList the list of indices that will be unpacked.
   * @param[in] recursive    whether or not to perform a recursive unpack.
   * @param[in] onDevice    whether to use device-based packing functions
   *                         (buffer must be either pinned or a device pointer)
   * @param[out] events      a collection of events to poll for completion of async
   *                         packing kernels ( device packing is incomplete until all
   *                         events are finalized )
   * @param[in] op           the operation to perform while unpacking
   * @return                 the number of bytes unpacked.
   *
   * This function takes a reference to a pointer to const buffer type, and
   * unpacks data from that buffer into the current Group. If @p packList
   * is non-empty, then a check is made to ensure that the data that is
   * unpacked matches @p packList. If @p packList is empty, the values
   * of the indices that are unpacked are stored and returned in packList.
   */
  virtual localIndex unpack( buffer_unit_type const * & buffer,
                             arrayView1d< localIndex > & packList,
                             integer const recursive,
                             bool onDevice,
                             parallelDeviceEvents & events,
                             MPI_Op op=MPI_REPLACE );

  ///@}

  //***********************************************************************************************

  //START_SPHINX_INCLUDE_GET_WRAPPER
  /**
   * @name Untyped wrapper retrieval methods
   *
   * These functions query the collection of Wrapper objects for the given
   * index/name/KeyIndex and returns a WrapperBase pointer to the object if
   * it exists. If it is not found, nullptr is returned.
   */
  ///@{

  /**
   * @brief Return a reference to a WrapperBase stored in this group.
   * @tparam KEY The lookup type.
   * @param key The value used to lookup the wrapper.
   * @return A reference to the WrapperBase that resulted from the lookup.
   * @throw std::domain_error if the wrapper doesn't exist.
   */
  template< typename KEY >
  WrapperBase const & getWrapperBase( KEY const & key ) const
  {
    WrapperBase const * const wrapper = m_wrappers[ key ];
    GEOS_THROW_IF( wrapper == nullptr,
                   "Group " << getDataContext() << " has no wrapper named " << key << std::endl
                            << dumpWrappersNames(),
                   std::domain_error );

    return *wrapper;
  }

  /**
   * @copydoc getWrapperBase(KEY const &) const
   */
  template< typename KEY >
  WrapperBase & getWrapperBase( KEY const & key )
  {
    WrapperBase * const wrapper = m_wrappers[ key ];
    GEOS_THROW_IF( wrapper == nullptr,
                   "Group " << getDataContext() << " has no wrapper named " << key << std::endl
                            << dumpWrappersNames(),
                   std::domain_error );

    return *wrapper;
  }

  /**
   * @brief
   * @param name
   * @return
   */
  indexType getWrapperIndex( string const & name ) const
  { return m_wrappers.getIndex( name ); }

  /**
   * @brief Get access to the internal wrapper storage.
   * @return a reference to wrapper map
   */
  wrapperMap const & wrappers() const
  { return m_wrappers; }

  /**
   * @copydoc wrappers() const
   */
  wrapperMap & wrappers()
  { return m_wrappers; }

  /**
   * @brief Return the number of wrappers.
   * @return The number of wrappers.
   */
  indexType numWrappers() const
  { return m_wrappers.size(); }

  /**
   * @return An array containing all wrappers keys
   */
  std::vector< string > getWrappersNames() const;

  ///@}

  /**
   * @name Typed wrapper retrieval methods
   *
   * These functions query the collection of Wrapper objects for the given
   * index/key and returns a Wrapper<T> pointer to the object if
   * it exists. The template parameter @p T is used to perform a cast
   * on the WrapperBase pointer that is returned by the lookup, into
   * a Wrapper<T> pointer. If the wrapper is not found, or the
   * WrapperBase pointer cannot be cast to a Wrapper<T> pointer, then nullptr
   * is returned.
   */
  ///@{

  /**
   * @brief Check if a wrapper exists
   * @tparam LOOKUP_TYPE the type of key used to perform the lookup.
   * @param[in] lookup a lookup value used to search the collection of wrappers
   * @return @p true if wrapper exists (regardless of type), @p false otherwise
   */
  template< typename LOOKUP_TYPE >
  bool hasWrapper( LOOKUP_TYPE const & lookup ) const
  { return m_wrappers[ lookup ] != nullptr; }

  /**
   * @brief Retrieve a Wrapper stored in this group.
   * @tparam T the object type contained in the Wrapper
   * @tparam LOOKUP_TYPE the type of key used to perform the lookup
   * @param[in] index    a lookup value used to search the collection of wrappers
   * @return A reference to the Wrapper<T> that resulted from the lookup.
   * @throw std::domain_error if the Wrapper doesn't exist.
   */
  template< typename T, typename LOOKUP_TYPE >
  Wrapper< T > const & getWrapper( LOOKUP_TYPE const & index ) const
  {
    WrapperBase const & wrapper = getWrapperBase( index );
    return dynamicCast< Wrapper< T > const & >( wrapper );
  }

  /**
   * @copydoc getWrapper(LOOKUP_TYPE const &) const
   */
  template< typename T, typename LOOKUP_TYPE >
  Wrapper< T > & getWrapper( LOOKUP_TYPE const & index )
  {
    WrapperBase & wrapper = getWrapperBase( index );
    return dynamicCast< Wrapper< T > & >( wrapper );
  }

  /**
   * @brief Retrieve a Wrapper stored in this group.
   * @tparam T the object type contained in the Wrapper
   * @tparam LOOKUP_TYPE the type of key used to perform the lookup
   * @param[in] index a lookup value used to search the collection of wrappers
   * @return A pointer to the Wrapper<T> that resulted from the lookup, if the Wrapper
   *   doesn't exist or has a different type a @c nullptr is returned.
   */
  template< typename T, typename LOOKUP_TYPE >
  Wrapper< T > const * getWrapperPointer( LOOKUP_TYPE const & index ) const
  { return dynamicCast< Wrapper< T > const * >( m_wrappers[ index ] ); }

  /**
   * @copydoc getWrapperPointer(LOOKUP_TYPE const &) const
   */
  template< typename T, typename LOOKUP_TYPE >
  Wrapper< T > * getWrapperPointer( LOOKUP_TYPE const & index )
  { return dynamicCast< Wrapper< T > * >( m_wrappers[ index ] ); }

  ///@}

  /**
   * @name Wrapper data access methods.
   *
   * These functions can be used to get referece/pointer access to the data
   * stored by wrappers in this group. They are essentially just shortcuts for
   * @p Group::getWrapper() and @p Wrapper<T>::getReference().
   * An additional template parameter can be provided to cast the return pointer
   * or reference to a base class pointer or reference (e.g. Array to ArrayView).
   */
  ///@{

  /**
   * @brief Look up a wrapper and get reference to wrapped object.
   * @tparam T return value type
   * @tparam WRAPPEDTYPE wrapped value type (by default, same as return)
   * @tparam LOOKUP_TYPE type of value used for wrapper lookup
   * @param lookup       value for wrapper lookup
   * @return reference to @p T
   * @throw A std::domain_error if the Wrapper does not exist.
   */
  template< typename T, typename LOOKUP_TYPE >
  GEOS_DECLTYPE_AUTO_RETURN
  getReference( LOOKUP_TYPE const & lookup ) const
  { return getWrapper< T >( lookup ).reference(); }

  /**
   * @copydoc getReference(LOOKUP_TYPE const &) const
   */
  template< typename T, typename LOOKUP_TYPE >
  T & getReference( LOOKUP_TYPE const & lookup )
  { return getWrapper< T >( lookup ).reference(); }

  //END_SPHINX_INCLUDE_GET_WRAPPER

  ///@}

  /**
   * @name Size/capacity management
   */
  /// @{

  /**
   * @brief Resize the group and all contained wrappers that resize with parent.
   * @param newSize the new size of the group
   */
  virtual void resize( localIndex const newSize );

  /**
   * @brief Set the new capacity and reserve it in all wrappers that resize with parent.
   * @param newsize the new capacity of the group
   */
  virtual void reserve( indexType const newsize );

  /**
   * @brief Get the "capacity" of the group, which determines the capacity of resizable wrappers.
   * @return capacity of this group
   */
  inline localIndex capacity() const
  { return m_capacity; }

  /**
   * @brief Get the "size" of the group, which determines the number of elements in resizable wrappers.
   * @return size of this group
   */
  inline localIndex size() const
  { return m_size; }

  /// @}

  /**
   * @name Basic group properties
   */
  ///@{

  /**
   * @brief Get group name.
   * @return group name
   */
  inline string const & getName() const
  { return m_name; }

  /**
   * @brief Return the path of this Group in the data repository.
   * Starts with '/' followed by the hierarchy of the children of the "Problem" in which the Group is.
   * @return The path of this group in the data repository.
   */
  string getPath() const;

  /**
   * @return DataContext object that that stores contextual information on this group that can be
   * used in output messages.
   */
  DataContext const & getDataContext() const
  { return *m_dataContext; }

  /**
   * @return DataContext object that that stores contextual information on a wrapper contained by
   * this group that can be used in output messages.
   * @tparam KEY The lookup type.
   * @param key The value used to lookup the wrapper.
   * @throw std::domain_error if the wrapper doesn't exist.
   */
  template< typename KEY >
  DataContext const & getWrapperDataContext( KEY key ) const
  { return getWrapperBase< KEY >( key ).getDataContext(); }

  /**
   * @brief Access the group's parent.
   * @return reference to parent Group
   * @throw std::domain_error if the Group doesn't have a parent.
   */
  Group & getParent()
  {
    GEOS_THROW_IF( m_parent == nullptr, "Group at " << getDataContext() << " does not have a parent.", std::domain_error );
    return *m_parent;
  }

  /**
   * @copydoc getParent()
   */
  Group const & getParent() const
  {
    GEOS_THROW_IF( m_parent == nullptr, "Group at " << getDataContext() << " does not have a parent.", std::domain_error );
    return *m_parent;
  }

  /**
   * @return true if this group has a parent.
   */
  bool hasParent() const
  { return m_parent != nullptr; }

  /**
   * @brief Get the group's index within its parent group
   * @return integral index of current group within its parent
   */
  localIndex getIndexInParent() const
  { return m_parent->getSubGroups().getIndex( m_name ); }

  /**
   * @brief Get the index of a sub-Group within this group.
   * @param key The key of the sub-Group
   * @return The index of the sub-Group.
   */
  localIndex getSubGroupIndex( keyType const & key ) const;

  /**
   * @brief Check whether this Group is resized when its parent is resized.
   * @return @p true if Group is resized with parent group, @p false otherwise
   */
  int sizedFromParent() const
  { return m_sizedFromParent; }

  /**
   * @brief Set whether this wrapper is resized when its parent is resized.
   * @param val an int that is converted into a bool
   * @return a pointer to this Group
   */
  Group & setSizedFromParent( int val )
  {
    m_sizedFromParent = val;
    return *this;
  }
  /**
   * @brief Get flags that control restart output of this group.
   * @return the current value of restart flags
   */
  RestartFlags getRestartFlags() const { return m_restart_flags; }

  /**
   * @brief Set flags that control restart output of this group.
   * @param flags the new value of restart flags
   */
  void setRestartFlags( RestartFlags flags ) { m_restart_flags = flags; }

  /**
   * @brief Get input flags for schema generation.
   * @return the current value of input flags
   */
  InputFlags getInputFlags() const { return m_input_flags; }

  /**
   * @brief Set input flags for schema generation.
   * @param flags the new value of input flags
   */
  void setInputFlags( InputFlags flags ) { m_input_flags = flags; }

  /**
   * @brief Structure to hold scoped key names
   */
  struct viewKeyStruct
  {
    /// @return String for the logLevel wrapper
    static constexpr char const * logLevelString() { return "logLevel"; }
  };

  ///@}

  /**
   * @brief Register a callback function on the group
   * @param func the function to register
   * @param funcType the type of the function to register
   * @return true if successful, false else
   */
  virtual bool registerCallback( void * func, const std::type_info & funcType )
  {
    GEOS_UNUSED_VAR( func );
    GEOS_UNUSED_VAR( funcType );
    return false;
  }

  /**
   * @name Restart output methods
   */
  ///@{

  /**
   * @brief Return the Conduit node object associated with this group.
   * @return The Conduit node object associated with this group.
   */
  conduit::Node & getConduitNode()
  { return m_conduitNode; }

  /// @copydoc getConduitNode()
  conduit::Node const & getConduitNode() const
  { return m_conduitNode; }

  /**
   * @brief Register the group and its wrappers with Conduit.
   */
  void prepareToWrite();

  /**
   * @brief Write the group and its wrappers into Conduit.
   */
  void finishWriting();

  /**
   * @brief Read the group and its wrappers from Conduit.
   */
  void loadFromConduit();

  /**
   * @deprecated will be remove and replace by addLogLevel
   */
  void enableLogLevelInput();

  /**
   * @brief Set verbosity level
   * @param logLevel new verbosity level value
   */
  void setLogLevel( integer const logLevel ) { m_logLevel = logLevel; }

  /// @return The verbosity level
  integer getLogLevel() const { return m_logLevel; }
  ///@}

  /**
   * @brief Performs re-initialization of certain variable depending on the solver being used.
   */
  virtual void reinit() {}

  /**
   * @brief Return PyGroup type.
   * @return Return PyGroup type.
   */
#if defined(GEOS_USE_PYGEOSX)
  virtual PyTypeObject * getPythonType() const;
#endif

protected:

  /**
   * @name Overridable initialization hooks
   *
   * These methods can be overridden by derived classes to customize
   * input post processing and object initialization.
   */
  ///@{

  /**
   * This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void postInputInitialization() {}

  /**
   * @brief Called by Initialize() prior to initializing sub-Groups.
   */
  virtual void initializePreSubGroups()
  {}

  /**
   * @brief Called by Initialize() after to initializing sub-Groups.
   */
  virtual void initializePostSubGroups()
  {}

  /**
   * @brief Called by InitializePostInitialConditions() prior to initializing sub-Groups.
   */
  virtual void initializePostInitialConditionsPreSubGroups()
  {}

  /**
   * @brief Called by InitializePostInitialConditions() after to initializing sub-Groups.
   */
  virtual void initializePostInitialConditionsPostSubGroups()
  {}

  /**
   * @brief Performs initialization required after reading from a restart file.
   */
  virtual void postRestartInitialization()
  {}



  ///@}

private:
  /**
   * @brief Read values from the input file and put them into the
   *   wrapped values for this group.
   * @param[in] xmlDocument the XML document that contains the targetNode
   * @param[in] targetNode the XML node that to extract input values from
   * @param[in] nodePos the target node position, typically obtained with xmlDocument::getNodePosition()
   */
  virtual void processInputFile( xmlWrapper::xmlNode const & targetNode,
                                 xmlWrapper::xmlNodePos const & nodePos );

  Group const & getBaseGroupByPath( string const & path ) const;

  /**
   * @brief Concrete implementation of the packing method.
   * @tparam DO_PACKING A template parameter to discriminate between actually packing or only computing the packing size.
   * @param[in,out] buffer The buffer that will receive the packed data.
   * @param[in] wrapperNames The names of the wrapper to be packed. If empty, all the wrappers will be packed.
   * @param[in] packList The element we want packed. If empty, all the elements will be packed.
   * @param[in] recursive Recursive pack or not.
   * @param[in] onDevice Whether to use device-based packing functions
   *                     (buffer must be either pinned or a device pointer)
   * @param[out] events A collection of events to poll for completion of async
   *                    packing kernels ( device packing is incomplete until all
   *                    events are finalized )
   * @return The packed size.
   */
  template< bool DO_PACKING >
  localIndex packImpl( buffer_unit_type * & buffer,
                       array1d< string > const & wrapperNames,
                       arrayView1d< localIndex const > const & packList,
                       integer const recursive,
                       bool onDevice,
                       parallelDeviceEvents & events ) const;

  //START_SPHINX_INCLUDE_02
  /// The parent Group that contains "this" Group in its "sub-Group" collection.
  Group * m_parent = nullptr;

  /// Specification that this group will have the same m_size as m_parent.
  integer m_sizedFromParent;

  /// The container for the collection of all wrappers continued in "this" Group.
  wrapperMap m_wrappers;

  /// The container for the collection of all sub-groups contained in "this" Group.
  subGroupMap m_subGroups;

  /// The size/length of this Group...and all Wrapper<> that are are specified to have the same size as their
  /// owning group.
  indexType m_size;

  /// The capacity for wrappers in this group...and all Wrapper<> that are specified to have the same size as their
  /// owning group.
  indexType m_capacity;

  /// The name/key of this Group in its parent collection of sub-Groups.
  string m_name;

  /// Verbosity flag for group logs
  integer m_logLevel;


  //END_SPHINX_INCLUDE_02

  /// Restart flag for this group... and subsequently all wrappers in this group.
  RestartFlags m_restart_flags;

  /// Input flag for this group.
  InputFlags m_input_flags;

  /// Reference to the conduit::Node that mirrors this group
  conduit::Node & m_conduitNode;

  // Keep track of log levels & descriptions
  std::unique_ptr< LogLevelsRegistry > m_logLevelsRegistry;

  /// A DataContext object used to provide contextual information on this Group,
  /// if it is created from an input XML file, the line or offset in that file.
  std::unique_ptr< DataContext > m_dataContext;

};

/**
 * @brief Type alias for KeyIndexT type used for sub-group lookups.
 */
using GroupKey = Group::subGroupMap::KeyIndex;

/**
 * @brief Type alias for KeyIndexT type used for wrapper lookups.
 */
using ViewKey = Group::wrapperMap::KeyIndex;

// Doxygen bug - sees this as a separate function
/// @cond DO_NOT_DOCUMENT
template< typename T, typename TBASE >
Wrapper< TBASE > & Group::registerWrapper( string const & name,
                                           ViewKey::index_type * const rkey )
{
  std::unique_ptr< TBASE > newObj = std::make_unique< T >();
  m_wrappers.insert( name,
                     new Wrapper< TBASE >( name, *this, std::move( newObj ) ),
                     true );

  if( rkey != nullptr )
  {
    *rkey = m_wrappers.getIndex( name );
  }

  Wrapper< TBASE > & rval = getWrapper< TBASE >( name );
  if( rval.sizedFromParent() == 1 )
  {
    rval.resize( size());
  }
  return rval;
}
/// @endcond

template< typename T, typename TBASE >
Wrapper< TBASE > & Group::registerWrapper( ViewKey const & viewKey )
{
  ViewKey::index_type index;
  Wrapper< TBASE > & rval = registerWrapper< T, TBASE >( viewKey.key(), &index );
  viewKey.setIndex( index );

  return rval;
}


template< typename T >
Wrapper< T > & Group::registerWrapper( string const & name,
                                       std::unique_ptr< T > newObject )
{
  static_assert( !std::is_base_of< WrapperBase, T >::value, "This function should not be used for `WrapperBase`. Use the dedicated `registerWrapper` instead." );
  m_wrappers.insert( name,
                     new Wrapper< T >( name, *this, std::move( newObject ) ),
                     true );

  Wrapper< T > & rval = getWrapper< T >( name );
  if( rval.sizedFromParent() == 1 )
  {
    rval.resize( size());
  }
  return rval;
}

template< typename T >
Wrapper< T > & Group::registerWrapper( string const & name,
                                       T * newObject )
{
  static_assert( !std::is_base_of< WrapperBase, T >::value, "This function should not be used for `WrapperBase`. Use the dedicated `registerWrapper` instead." );
  m_wrappers.insert( name,
                     new Wrapper< T >( name, *this, newObject ),
                     true );

  Wrapper< T > & rval = getWrapper< T >( name );
  if( rval.sizedFromParent() == 1 )
  {
    rval.resize( size());
  }
  return rval;
}

template< typename LOG_LEVEL_INFO >
std::enable_if_t< geos::is_log_level_info< LOG_LEVEL_INFO >, void >
Group::addLogLevel()
{
  GEOS_ERROR_IF( m_logLevelsRegistry == nullptr, "You cannot call addLogLevel after schema generation" );

  Wrapper< integer > * wrapper = getWrapperPointer< integer >( viewKeyStruct::logLevelString() );
  if( wrapper == nullptr )
  {
    wrapper = &registerWrapper( viewKeyStruct::logLevelString(), &m_logLevel );
    wrapper->setApplyDefaultValue( 0 );
    wrapper->setInputFlag( InputFlags::OPTIONAL );
  }
  m_logLevelsRegistry->addEntry( LOG_LEVEL_INFO::getMinLogLevel(),
                                 LOG_LEVEL_INFO::getDescription() );
  wrapper->setDescription( m_logLevelsRegistry->buildLogLevelDescription());
}

} /* end namespace dataRepository */
} /* end namespace geos */

#endif /* GEOS_DATAREPOSITORY_GROUP_HPP_ */
