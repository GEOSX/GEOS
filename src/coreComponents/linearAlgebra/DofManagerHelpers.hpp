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
 * @file DofManagerHelpers.hpp
 */

#ifndef GEOS_LINEARALGEBRA_DOFMANAGERHELPERS_HPP
#define GEOS_LINEARALGEBRA_DOFMANAGERHELPERS_HPP

#include "mesh/utilities/MeshMapUtilities.hpp"

namespace geos
{

// unnamed namespace to avoid needless external linkage
namespace
{

/**
 * @brief A struct to abstract away some details of mesh access.
 * @tparam LOC type of mesh location
 */
template< FieldLocation LOC >
struct MeshHelper;

template<>
struct MeshHelper< FieldLocation::Node >
{
  using ManagerType = NodeManager;
  using LocalIndexType = localIndex;

  static LocalIndexType constexpr invalid_local_index{ -1 };

  static constexpr char const * managerGroupName() { return MeshLevel::groupStructKeys::nodeManagerString(); }
  static constexpr char const * mapViewKey() { return ElementSubRegionBase::viewKeyStruct::nodeListString(); }
  static constexpr char const * syncObjName = "node";

  template< typename MANAGER >
  using MapType = typename MANAGER::NodeMapType;
};

template<>
struct MeshHelper< FieldLocation::Edge >
{
  using ManagerType = EdgeManager;
  using LocalIndexType = localIndex;

  static LocalIndexType constexpr invalid_local_index{ -1 };

  static constexpr char const * managerGroupName() { return MeshLevel::groupStructKeys::edgeManagerString(); }
  static constexpr char const * mapViewKey() { return ElementSubRegionBase::viewKeyStruct::edgeListString(); }
  static constexpr char const * syncObjName = "edge";

  template< typename MANAGER >
  using MapType = typename MANAGER::EdgeMapType;
};

template<>
struct MeshHelper< FieldLocation::Face >
{
  using ManagerType = FaceManager;
  using LocalIndexType = localIndex;

  static LocalIndexType constexpr invalid_local_index{ -1 };

  static constexpr char const * managerGroupName() { return MeshLevel::groupStructKeys::faceManagerString(); }
  static constexpr char const * mapViewKey() { return ElementSubRegionBase::viewKeyStruct::faceListString(); }
  static constexpr char const * syncObjName = "face";

  template< typename MANAGER >
  using MapType = typename MANAGER::FaceMapType;
};

template<>
struct MeshHelper< FieldLocation::Elem >
{
  using ManagerType = ElementSubRegionBase;
  using LocalIndexType = std::array< localIndex, 3 >;

  static LocalIndexType constexpr invalid_local_index{ -1, -1, -1 };

  static constexpr auto managerGroupName() { return MeshLevel::groupStructKeys::elemManagerString(); }
  static constexpr auto syncObjName = "elems";

  template< typename MANAGER >
  using MapType = typename MANAGER::ElemMapType;
};

/**
 * @brief Compile-time limits on incidence between mesh objects.
 * @tparam LOC primary mesh object (location) type
 * @tparam LOC_ADJ adjacent mesh pbject type
 *
 * Typical usages:
 * - preallocate space in sparsity patterns
 * - choose size for stack arrays at compile time
 *
 * These limits may need to increase as we encounter more complex meshes.
 */
template< FieldLocation LOC, FieldLocation LOC_ADJ >
struct MeshIncidence {};

/// Self-adjacency is always 1.
template< FieldLocation LOC >
struct MeshIncidence< LOC, LOC >
{
  static localIndex constexpr max = 1;
};

/// Shortcut macro for declaring adjacency.
#define SET_MAX_MESH_INCIDENCE( LOC, LOC_ADJ, MAX ) \
  template<> \
  struct MeshIncidence< FieldLocation::LOC, FieldLocation::LOC_ADJ > \
  { \
    static localIndex constexpr max = MAX; \
  }

SET_MAX_MESH_INCIDENCE( Node, Edge, 10 );
SET_MAX_MESH_INCIDENCE( Node, Face, 10 );
SET_MAX_MESH_INCIDENCE( Node, Elem, 12 );

SET_MAX_MESH_INCIDENCE( Edge, Node, 2 );
SET_MAX_MESH_INCIDENCE( Edge, Face, 8 );
SET_MAX_MESH_INCIDENCE( Edge, Elem, 8 );

SET_MAX_MESH_INCIDENCE( Face, Node, 8 );
SET_MAX_MESH_INCIDENCE( Face, Edge, 4 );
SET_MAX_MESH_INCIDENCE( Face, Elem, 2 );

SET_MAX_MESH_INCIDENCE( Elem, Node, 8 );
SET_MAX_MESH_INCIDENCE( Elem, Edge, 12 );
SET_MAX_MESH_INCIDENCE( Elem, Face, 6 );

#undef SET_MAX_MESH_INCIDENCE

/**
 * @brief Switchyard for Location (Node/Edge/Face/Elem)
 * @tparam LAMBDA generic functor type
 * @param loc location type
 * @param lambda functor to be called
 */
template< typename LAMBDA >
bool LocationSwitch( FieldLocation const loc,
                     LAMBDA lambda )
{
  switch( loc )
  {
    case  FieldLocation::Node:
    {
      lambda( std::integral_constant< FieldLocation, FieldLocation::Node >() );
      return true;
    }
    case  FieldLocation::Edge:
    {
      lambda( std::integral_constant< FieldLocation, FieldLocation::Edge >() );
      return true;
    }
    case  FieldLocation::Face:
    {
      lambda( std::integral_constant< FieldLocation, FieldLocation::Face >() );
      return true;
    }
    case  FieldLocation::Elem:
    {
      lambda( std::integral_constant< FieldLocation, FieldLocation::Elem >() );
      return true;
    }
    default:
      return false;
  }
}

template< typename LAMBDA >
bool LocationSwitch( FieldLocation const loc1,
                     FieldLocation const loc2,
                     LAMBDA lambda )
{
  bool ret2;
  bool const ret1 =
    LocationSwitch( loc1, [&]( auto const loc_type1 )
  {
    ret2 =
      LocationSwitch( loc2, [&]( auto const loc_type2 )
    {
      lambda( loc_type1, loc_type2 );
    } );
  } );
  return ret1 && ret2;
}

template< FieldLocation LOC >
typename MeshHelper< LOC >::ManagerType const & getObjectManager( MeshLevel const & mesh )
{
  using ObjectManager = typename MeshHelper< LOC >::ManagerType;
  return mesh.getGroup< ObjectManager >( MeshHelper< LOC >::managerGroupName() );
}

template< FieldLocation LOC >
typename MeshHelper< LOC >::ManagerType & getObjectManager( MeshLevel & mesh )
{
  using ObjectManager = typename MeshHelper< LOC >::ManagerType;
  return mesh.getGroup< ObjectManager >( MeshHelper< LOC >::managerGroupName() );
}

ObjectManagerBase const & getObjectManager( FieldLocation const loc,
                                            MeshLevel const & mesh )
{
  ObjectManagerBase const * manager = nullptr;
  LocationSwitch( loc, [&]( auto const LOC )
  {
    manager = &getObjectManager< decltype(LOC)::value >( mesh );
  } );
  GEOS_ASSERT( manager != nullptr );
  return *manager;
}

ObjectManagerBase & getObjectManager( FieldLocation const loc,
                                      MeshLevel & mesh )
{
  return const_cast< ObjectManagerBase & >( getObjectManager( loc, const_cast< MeshLevel const & >( mesh ) ) );
}

/**
 * @brief Helper alias for accessing mesh maps.
 * @tparam LOC typeof object (location)
 * @tparam MANAGER type of source object manager
 */
template< FieldLocation LOC, typename MANAGER >
using MapType = typename MeshHelper< LOC >::template MapType< MANAGER >;

/**
 * @brief A helper struct implementing visitation of mesh objects and connected (adjacent) objects
 * @tparam LOC type of primary mesh object (location) to visit
 * @tparam CONN_LOC type of adjacent mesh object to visit from primary; can be the same as primary,
 *                  in which case each is visited only once and the same index is passed for both.
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 */
template< FieldLocation LOC, FieldLocation CONN_LOC, bool VISIT_GHOSTS >
struct MeshVisitor;

/**
 * @brief A specialization of MeshLoopHelper used when adjacent objects are NOT needed
 *        and the primary objects are NOT elements.
 * @tparam LOC type of primary mesh object (location) to visit
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 *
 * The algorithm will visit objects in order of increasing local index.
 */
template< FieldLocation LOC, bool VISIT_GHOSTS >
struct MeshVisitor< LOC, LOC, VISIT_GHOSTS >
{
  template< typename POLICY, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER, typename LAMBDA >
  static void visit( MeshLevel const & meshLevel,
                     REGIONS_CONTAINER const & regions,
                     LAMBDA && lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper< LOC >::ManagerType;

    // get access to location ghost rank (we don't want to visit ghosted locations
    ObjectManagerLoc const & objectManager = getObjectManager< LOC >( meshLevel );
    arrayView1d< integer const > const & ghostRank = objectManager.ghostRank();

    // create an array to track previously visited locations (to avoid multiple visits)
    array1d< integer > visited( objectManager.size() );

    meshLevel.getElemManager().
      forElementSubRegions< SUBREGIONTYPES... >( regions, [&]( localIndex const, auto const & subRegion )
    {
      // derive some more useful, subregion-dependent type aliases
      using ElemToLocMapType = MapType< LOC, TYPEOFREF( subRegion ) >;

      // get access to element-to-location map
      auto const & elemToLocMap =
        subRegion.template getReference< ElemToLocMapType >( MeshHelper< LOC >::mapViewKey() ).toViewConst();

      // loop over all elements (including ghosts, which may be necessary to access some locally owned locations)
      forAll< parallelHostPolicy >( subRegion.size(), [=, visited = visited.toView()]( localIndex const ei )
      {
        // loop over all locations incident on an element and increment visitation counter
        auto const locList = elemToLocMap[ei];
        for( localIndex a = 0; a < locList.size(); ++a )
        {
          localIndex const locIdx = locList[a];
          if( VISIT_GHOSTS || ghostRank[locIdx] < 0 )
          {
            RAJA::atomicInc< parallelHostAtomic >( &visited[locIdx] );
          }
        }
      } );
    } );

    // turn marker array into a list of indices to visit (this is inherently serial, but fast)
    array1d< localIndex > locations;
    locations.reserve( objectManager.size() );
    for( localIndex i = 0; i < visited.size(); ++i )
    {
      if( visited[i] > 0 )
      {
        locations.emplace_back( i );
      }
    }

    // optimize for the common case (e.g. all nodes of the mesh) to avoid extra memory access
    if( locations.size() == objectManager.size() )
    {
      forAll< POLICY >( objectManager.size(), [=]( localIndex const locIdx )
      {
        lambda( locIdx, locIdx, 0 );
      } );
    }
    else
    {
      forAll< POLICY >( locations.size(), [=, locations = locations.toViewConst()]( localIndex const i )
      {
        localIndex const locIdx = locations[i];
        lambda( locIdx, locIdx, 0 );
      } );
    }
  }
};

/**
 * @brief A specialization of MeshLoopHelper used when NEITHER primary NOR adjacent objects are elements
 * @tparam LOC type of primary mesh object (location) to visit
 * @tparam CONN_LOC type of adjacent mesh object to visit from primary
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 */
template< FieldLocation LOC, FieldLocation CONN_LOC, bool VISIT_GHOSTS >
struct MeshVisitor
{
  template< typename POLICY, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER, typename LAMBDA >
  static void visit( MeshLevel const & meshLevel,
                     REGIONS_CONTAINER const & regions,
                     LAMBDA && lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper< LOC >::ManagerType;
    using LocToConnMapType = MapType< CONN_LOC, ObjectManagerLoc >;

    ObjectManagerLoc const & objectManager = getObjectManager< LOC >( meshLevel );

    // get access to location-to-connected map
    auto const & locToConnMap =
      objectManager.template getReference< LocToConnMapType >( MeshHelper< CONN_LOC >::mapViewKey() ).toViewConst();

    // call the specialized version first, then add an extra loop over connected objects
    MeshVisitor< LOC, LOC, VISIT_GHOSTS >::
    template visit< POLICY, SUBREGIONTYPES... >( meshLevel, regions,
                                                 [=]( localIndex const locIdx,
                                                      localIndex const,
                                                      localIndex const )
    {
      // loop over all connected locations
      auto const connList = locToConnMap[locIdx];
      for( localIndex b = 0; b < connList.size(); ++b )
      {
        lambda( locIdx, connList[b], b );
      }
    } );
  }
};

/**
 * @brief A specialization of MeshLoopHelper used when primary objects are NOT elements, but adjacent are
 * @tparam LOC type of primary mesh object (location) to visit
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 *
 * @note We have to have this specialization because to-element maps are different from all other
 *       map types in that they are stored in data repository as three separate arrays.
 */
template< FieldLocation LOC, bool VISIT_GHOSTS >
struct MeshVisitor< LOC, FieldLocation::Elem, VISIT_GHOSTS >
{
  template< typename POLICY, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER, typename LAMBDA >
  static void visit( MeshLevel const & meshLevel,
                     REGIONS_CONTAINER const & regions,
                     LAMBDA && lambda )
  {
    // derive some useful type aliases
    using ObjectManagerLoc = typename MeshHelper< LOC >::ManagerType;
    using ToElemMapType = typename MapType< FieldLocation::Elem, ObjectManagerLoc >::base_type;

    // get mesh object manager for LOC to access maps
    ObjectManagerLoc const & objectManager = getObjectManager< LOC >( meshLevel );

    // access to location-to-element map
    using keys = typename ObjectManagerLoc::viewKeyStruct;
    auto const & elemRegionList =
      objectManager.template getReference< ToElemMapType >( keys::elementRegionListString() ).toViewConst();
    auto const & elemSubRegionList =
      objectManager.template getReference< ToElemMapType >( keys::elementSubRegionListString() ).toViewConst();
    auto const & elemIndexList =
      objectManager.template getReference< ToElemMapType >( keys::elementListString() ).toViewConst();

    // call the specialized version first, then add an extra loop over connected elements
    MeshVisitor< LOC, LOC, VISIT_GHOSTS >::
    template visit< POLICY, SUBREGIONTYPES... >( meshLevel, regions,
                                                 [=]( localIndex const locIdx,
                                                      localIndex const,
                                                      localIndex const )
    {
      // loop over all connected locations
      auto const elemRegions = elemRegionList[ locIdx ];
      auto const elemSubRegions = elemSubRegionList[ locIdx ];
      auto const elemIndices = elemIndexList[ locIdx ];
      for( localIndex b = 0; b < elemIndices.size(); ++b )
      {
        MeshHelper< FieldLocation::Elem >::LocalIndexType const elemIdx =
        { elemRegions[b], elemSubRegions[b], elemIndices[b] };

        if( elemIdx[0] >= 0 && elemIdx[1] >= 0 && elemIdx[2] >= 0 )
        {
          lambda( locIdx, elemIdx, b );
        }
      }
    } );
  }
};

/**
 * @brief A specialization of MeshLoopHelper used when primary objects are elements and adjacent are NOT
 * @tparam CONN_LOC type of adjacent mesh objects (locations)
 * @tparam VISIT_GHOSTS whether to visit ghosted primary locations
 *
 * @note We have to have this specialization because various subregion types might have different
 *       map types, so the correct map type must be derived within the subregion loop
 */
template< FieldLocation CONN_LOC, bool VISIT_GHOSTS >
struct MeshVisitor< FieldLocation::Elem, CONN_LOC, VISIT_GHOSTS >
{
  template< typename POLICY, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER, typename LAMBDA >
  static void visit( MeshLevel const & meshLevel,
                     REGIONS_CONTAINER const & regions,
                     LAMBDA && lambda )
  {
    meshLevel.getElemManager().
      forElementSubRegionsComplete< SUBREGIONTYPES... >( regions, [&]( localIndex const,
                                                                       localIndex const er,
                                                                       localIndex const esr,
                                                                       ElementRegionBase const &,
                                                                       auto const & subRegion )
    {
      // derive subregion-dependent map type
      using ElemToConnMapType = MapType< CONN_LOC, TYPEOFREF( subRegion ) >;

      // get access to element-to-location map
      auto const & elemToConnMap =
        subRegion.template getReference< ElemToConnMapType >( MeshHelper< CONN_LOC >::mapViewKey() ).toViewConst();

      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

      forAll< POLICY >( subRegion.size(), [=]( localIndex const ei )
      {
        if( VISIT_GHOSTS || ghostRank[ei] < 0 )
        {
          auto const elemIdx = MeshHelper< FieldLocation::Elem >::LocalIndexType{ er, esr, ei };
          auto const connList = elemToConnMap[ei];

          for( localIndex a = 0; a < connList.size(); ++a )
          {
            lambda( elemIdx, connList[a], a );
          }
        }
      } );
    } );
  }
};

/**
 * @brief A specialization of MeshLoopHelper used when only loop over elements is required
 * @tparam VISIT_GHOSTS whether to visit ghosted elements
 */
template< bool VISIT_GHOSTS >
struct MeshVisitor< FieldLocation::Elem, FieldLocation::Elem, VISIT_GHOSTS >
{
  template< typename POLICY, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER, typename LAMBDA >
  static void visit( MeshLevel const & meshLevel,
                     REGIONS_CONTAINER const & regions,
                     LAMBDA && lambda )
  {
    meshLevel.getElemManager().
      forElementSubRegionsComplete< SUBREGIONTYPES... >( regions, [&]( localIndex const,
                                                                       localIndex const er,
                                                                       localIndex const esr,
                                                                       ElementRegionBase const &,
                                                                       ElementSubRegionBase const & subRegion )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

      forAll< POLICY >( subRegion.size(), [=]( localIndex const ei )
      {
        if( VISIT_GHOSTS || ghostRank[ei] < 0 )
        {
          MeshHelper< FieldLocation::Elem >::LocalIndexType elemIdx{ er, esr, ei };
          lambda( elemIdx, elemIdx, 0 );
        }
      } );
    } );
  }
};

/**
 * @brief Visit mesh locations within active regions and their connected locations (e.g. nodes of an element)
 *
 * @tparam LOC type of location (Node/Edge/Face/Elem)
 * @tparam CONN_LOC type of connected location (Node/Edge/Face/Elem)
 * @tparam SUBREGIONTYPES variadic pack of subregion types to visit;
 *         if omitted, uses ElementRegionManager's default behavior
 * @tparam REGIONS_CONTAINER type of region list container
 * @tparam LAMBDA type of user-provided functor or lambda to call
 * @param mesh the mesh to loop over
 * @param regions list of region names to visit (assumed unique)
 * @param lambda functor or lambda to call
 *
 * Calls user-provided lambda with 3 parameters:
 * - mesh local index of primary location (one index for node/edge/face, reg/subreg/index triplet for elems)
 * - mesh local index of connected location
 * - index of connected location inside primary location's map (e.g. which node of an element)
 * If CONN_LOC == LOC, lambda is called only once per location, as lambda(loc, loc, 0)
 *
 * @note This function will exclude ghosted locations, but will not perform a similar check for
 *       connected locations - these may include ghosts, it is up to the user to filter as needed.
 *       Similarly, while primary loop is limited to @p regions, adjacent locations may not belong.
 */
template< FieldLocation LOC, FieldLocation CONN_LOC, bool VISIT_GHOSTS,
          typename POLICY, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER, typename LAMBDA >
void forMeshLocation( MeshLevel const & mesh,
                      REGIONS_CONTAINER const & regions,
                      LAMBDA && lambda )
{
  MeshVisitor< LOC, CONN_LOC, VISIT_GHOSTS >::
  template visit< POLICY, SUBREGIONTYPES... >( mesh,
                                               regions,
                                               std::forward< LAMBDA >( lambda ) );
}

/**
 * @brief A shortcut for forMeshLocation with CONN_LOC == LOC
 */
template< FieldLocation LOC, bool VISIT_GHOSTS,
          typename POLICY, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER, typename LAMBDA >
void forMeshLocation( MeshLevel const & mesh,
                      REGIONS_CONTAINER const & regions,
                      LAMBDA && lambda )
{
  forMeshLocation< LOC, LOC, VISIT_GHOSTS, POLICY, SUBREGIONTYPES... >( mesh,
                                                                        regions,
                                                                        [=]( auto const locIdx, auto, auto )
  {
    lambda( locIdx );
  } );
}

/**
 * @brief Count objects of given type within a set of regions
 * @tparam LOC type of mesh objects (locations)
 * @tparam VISIT_GHOSTS whether ghosted objects should be counted
 * @tparam SUBREGIONTYPES types of subregions to include in counting
 * @tparam REGIONS_CONTAINER type of region list container
 * @param mesh the mesh level to use
 * @param regions a list of region names (assumed to be unique)
 * @return number of mesh objects adjacent to given regions
 */
template< FieldLocation LOC, bool VISIT_GHOSTS, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER >
localIndex countMeshObjects( MeshLevel const & mesh,
                             REGIONS_CONTAINER const & regions )
{
  using countPolicy = parallelHostPolicy;
  RAJA::ReduceSum< ReducePolicy< countPolicy >, localIndex > count( 0 );
  forMeshLocation< LOC, VISIT_GHOSTS, countPolicy, SUBREGIONTYPES... >( mesh, regions, [=]( auto )
  {
    count += 1;
  } );
  return count.get();
}

/**
 * @brief Count objects of given type within a set of regions
 * @tparam VISIT_GHOSTS whether ghosted objects should be counted
 * @tparam SUBREGIONTYPES types of subregions to include in counting
 * @tparam REGIONS_CONTAINER type of region list container
 * @param mesh the mesh level to use
 * @param regions a list of region names (assumed to be unique)
 * @param loc type of mesh objects (locations)
 * @return number of mesh objects adjacent to given regions
 */
template< bool VISIT_GHOSTS, typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER >
localIndex countMeshObjects( FieldLocation const location,
                             MeshLevel const & mesh,
                             REGIONS_CONTAINER const & regions )
{
  localIndex count = 0;
  bool const success = LocationSwitch( location, [&]( auto const loc )
  {
    FieldLocation constexpr LOC = decltype(loc)::value;
    count = countMeshObjects< LOC, VISIT_GHOSTS, SUBREGIONTYPES... >( mesh, regions );
  } );
  GEOS_ERROR_IF( !success, "Invalid location type: " << static_cast< int >( location ) );
  return count;
}

/**
 * @brief This template abstracts some differences between index arrays that live
 *        on elements vs other mesh objects: (de)registration, value access, etc.
 * @tparam T type of index value
 * @tparam LOC target location
 */
template< typename T, FieldLocation LOC >
struct ArrayHelper
{
  using ArrayType = array1d< std::remove_const_t< T > >;
  using ViewType = arrayView1d< T >;
  using Accessor = ViewType;

  template< typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER >
  static void
  create( MeshLevel & mesh,
          string const & key,
          string const & description,
          REGIONS_CONTAINER const & GEOS_UNUSED_PARAM( regions ) )
  {
    ObjectManagerBase & baseManager = getObjectManager< LOC >( mesh );
    baseManager.registerWrapper< ArrayType >( key ).
      setApplyDefaultValue( -1 ).
      setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
      setRestartFlags( dataRepository::RestartFlags::NO_WRITE ).
      setDescription( description );
  }

  static Accessor get( add_const_if_t< MeshLevel, std::is_const< T >::value > & mesh, string const & key )
  {
    return getObjectManager< LOC >( mesh ).template getReference< ArrayType >( key );
  }

  static inline T value( Accessor const & indexArray, typename MeshHelper< LOC >::LocalIndexType const i )
  {
    return indexArray[i];
  }

  static inline T & reference( Accessor const & indexArray, typename MeshHelper< LOC >::LocalIndexType const i )
  {
    return indexArray[i];
  }

  template< typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER >
  static void
  remove( MeshLevel & mesh,
          string const & key,
          REGIONS_CONTAINER const & GEOS_UNUSED_PARAM( regions ) )
  {
    getObjectManager< LOC >( mesh ).deregisterWrapper( key );
  }
};

/**
 * @brief Specialized version of ArrayHelper for element-based arrays.
 * @tparam T type of index value
 */
template< typename T >
struct ArrayHelper< T, FieldLocation::Elem >
{
  using ArrayType = array1d< std::remove_const_t< T > >;
  using ViewType = arrayView1d< T >;
  using Accessor = ElementRegionManager::ElementViewAccessor< ViewType >;

  template< typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER >
  static void
  create( MeshLevel & mesh,
          string const & key,
          string const & description,
          REGIONS_CONTAINER const & regions )
  {
    mesh.getElemManager().template forElementSubRegions< SUBREGIONTYPES... >( regions,
                                                                              [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< ArrayType >( key ).
        setApplyDefaultValue( -1 ).
        setPlotLevel( dataRepository::PlotLevel::LEVEL_1 ).
        setRestartFlags( dataRepository::RestartFlags::NO_WRITE ).
        setDescription( description );
    } );
  }

  static Accessor get( add_const_if_t< MeshLevel, std::is_const< T >::value > & mesh, string const & key )
  {
    return mesh.getElemManager().template constructViewAccessor< ArrayType, ViewType >( key );
  }

  static inline T value( Accessor const & indexArray,
                         MeshHelper< FieldLocation::Elem >::LocalIndexType const & e )
  {
    if( indexArray[std::get< 0 >( e )].empty() || indexArray[std::get< 0 >( e )][std::get< 1 >( e )].empty() )
    {
      return -1;
    }
    return indexArray[std::get< 0 >( e )][std::get< 1 >( e )][std::get< 2 >( e )];
  }

  static inline T & reference( Accessor const & indexArray,
                               MeshHelper< FieldLocation::Elem >::LocalIndexType const & e )
  {
    return indexArray[std::get< 0 >( e )][std::get< 1 >( e )][std::get< 2 >( e )];
  }

  template< typename ... SUBREGIONTYPES, typename REGIONS_CONTAINER >
  static void
  remove( MeshLevel & mesh,
          string const & key,
          REGIONS_CONTAINER const & regions )
  {
    mesh.getElemManager().template forElementSubRegions< SUBREGIONTYPES... >( regions,
                                                                              [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      subRegion.deregisterWrapper( key );
    } );
  }
};

} // namespace

} // namespace geos

#endif //GEOS_LINEARALGEBRA_DOFMANAGERHELPERS_HPP
