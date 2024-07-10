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
 * @file testDofManagerUtils.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP_
#define GEOS_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP_

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"

#include <gtest/gtest.h>

namespace geos
{

namespace testing
{

/**
 * @brief Set up a problem from an xml input buffer
 * @param problemManager the target problem manager
 * @param xmlInput       the XML input string
 */
void setupProblemFromXML( ProblemManager * const problemManager, char const * const xmlInput )
{
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.loadString( xmlInput );
  if( !xmlResult )
  {
    GEOS_LOG_RANK_0( "XML parsed with errors!" );
    GEOS_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  int mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  dataRepository::Group & commandLine =
    problemManager->getGroup< dataRepository::Group >( problemManager->groupKeys.commandLine );
  commandLine.registerWrapper< integer >( problemManager->viewKeys.xPartitionsOverride.key() ).
    setApplyDefaultValue( mpiSize );

  xmlWrapper::xmlNode xmlProblemNode = xmlDocument.getChild( dataRepository::keys::ProblemManager );
  problemManager->processInputFileRecursive( xmlDocument, xmlProblemNode );

  // Open mesh levels
  DomainPartition & domain = problemManager->getDomainPartition();
  MeshManager & meshManager = problemManager->getGroup< MeshManager >( problemManager->groupKeys.meshManager );
  meshManager.generateMeshLevels( domain );

  ElementRegionManager & elementManager = domain.getMeshBody( 0 ).getBaseDiscretization().getElemManager();
  xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( elementManager.getName().c_str() );
  elementManager.processInputFileRecursive( xmlDocument, topLevelNode );
  elementManager.postInputInitializationRecursive();

  problemManager->problemSetup();
  problemManager->applyInitialConditions();
}

/**
 * @brief Get a region list from an input list (empty region lists to mean all regions in the mesh).
 * @param mesh  pointer to the mesh
 * @param input input list of region names (may be empty)
 * @return the list of region names (same as input unless input is empty)
 *
 * Mainly used to allow empty region lists to mean all regions.
 */
std::vector< DofManager::FieldSupport > getRegions( DomainPartition const & domain, std::vector< DofManager::FieldSupport > const & input )
{
  std::vector< DofManager::FieldSupport > regions( input.begin(), input.end() );
  for( DofManager::FieldSupport & support : regions )
  {
    if( support.regionNames.empty() )
    {
      MeshBody const & meshBody = domain.getMeshBody( support.meshBodyName );
      MeshLevel const & meshLevel = meshBody.getMeshLevel( support.meshLevelName );
      meshLevel.getElemManager().forElementRegions( [&]( ElementRegionBase const & region )
      {
        support.regionNames.insert( region.getName() );
      } );
    }
  }
  return regions;
}

namespace internal
{

template< FieldLocation LOC >
struct testMeshHelper {};

template<>
struct testMeshHelper< FieldLocation::Node >
{
  static constexpr auto managerKey()
  { return MeshLevel::groupStructKeys::nodeManagerString(); }

  static constexpr auto elemMapKey()
  { return ElementSubRegionBase::viewKeyStruct::nodeListString(); }

  template< typename SUBREGION > using ElemToObjMap = typename SUBREGION::NodeMapType;
};

template<>
struct testMeshHelper< FieldLocation::Face >
{
  static constexpr auto managerKey()
  { return MeshLevel::groupStructKeys::faceManagerString(); }

  static constexpr auto elemMapKey()
  { return ElementSubRegionBase::viewKeyStruct::faceListString(); }

  template< typename SUBREGION > using ElemToObjMap = typename SUBREGION::FaceMapType;
};

template< int USD >
localIndex size1( arrayView2d< localIndex const, USD > const & map,
                  localIndex const GEOS_UNUSED_PARAM( i0 ) )
{
  return map.size( 1 );
}

localIndex size1( arrayView1d< arrayView1d< localIndex const > const > const & map, localIndex const i0 )
{
  return map[i0].size();
}

localIndex size1( ArrayOfArraysView< localIndex const > const & map, localIndex const i0 )
{
  return map.sizeOfArray( i0 );
}

template< FieldLocation LOC >
struct forLocalObjectsImpl
{
  template< typename LAMBDA, typename REGIONS_CONTAINER >
  static void f( MeshLevel const & mesh,
                 REGIONS_CONTAINER const & regions,
                 LAMBDA lambda )
  {
    using helper = testMeshHelper< LOC >;
    ObjectManagerBase const & manager = mesh.getGroup< ObjectManagerBase >( helper::managerKey() );

    arrayView1d< integer const > ghostRank = manager.ghostRank();

    array1d< bool > visited( ghostRank.size() );

    mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex const, auto const & subRegion )
    {
      using MapType = typename helper::template ElemToObjMap< std::remove_reference_t< decltype( subRegion ) > >;

      traits::ViewTypeConst< MapType > const
      elemToObjMap = subRegion.template getReference< MapType >( helper::elemMapKey() );

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( localIndex a = 0; a < size1( elemToObjMap, k ); ++a )
        {
          localIndex const idx = elemToObjMap[k][a];
          if( ghostRank[idx] < 0 && !std::exchange( visited[idx], true ) )
          {
            lambda( idx );
          }
        }
      }
    } );
  }
};

template<>
struct forLocalObjectsImpl< FieldLocation::Elem >
{
  template< typename LAMBDA, typename REGIONS_CONTAINER >
  static void f( MeshLevel const & mesh,
                 REGIONS_CONTAINER const & regions,
                 LAMBDA lambda )
  {
    // make a list of regions
    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegionsComplete< ElementSubRegionBase >( regions,
                                                                      [&]( localIndex const,
                                                                           localIndex const er,
                                                                           localIndex const esr,
                                                                           ElementRegionBase const &,
                                                                           ElementSubRegionBase const & subRegion )
    {
      arrayView1d< integer const > ghostRank = subRegion.ghostRank();

      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( ghostRank[ei] < 0 )
        {
          lambda( std::array< localIndex, 3 > { er, esr, ei} );
        }
      }
    } );
  }
};

} // namespace internal

/**
 * @brief Apply a lambda to all locally owned mesh objects in the mesh.
 * @tparam LOC type of mesh location (Node, Element, etc.)
 * @tparam LAMBDA type of lambda
 * @param mesh    pointer to the mesh
 * @param regions list of input region names to loop over
 * @param lambda  the lambda to apply
 */
template< FieldLocation LOC, typename REGIONS_CONTAINER, typename LAMBDA >
void forLocalObjects( MeshLevel const & mesh,
                      REGIONS_CONTAINER const & regions,
                      LAMBDA && lambda )
{
  internal::forLocalObjectsImpl< LOC >::template f( mesh, regions, std::forward< LAMBDA >( lambda ) );
}

/**
 * @brief Count the number of local objects in the mesh.
 * @tparam LOC type of mesh location (Node, Element, etc.)
 * @param mesh    pointer to the mesh
 * @param regions list of input region names to loop over
 * @return the number of locally owned objects (e.g. nodes)
 */
template< FieldLocation LOC >
localIndex countLocalObjects( DomainPartition const & domain, std::vector< DofManager::FieldSupport > const & support )
{
  localIndex numLocal = 0;
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & meshLevel = meshBody.getMeshLevel( regions.meshLevelName );

    forLocalObjects< LOC >( meshLevel, regions.regionNames, [&]( auto const & ) { ++numLocal; } );
  }
  return numLocal;
}

/**
 * @brief Create a TPFA-type sparsity pattern.
 * @param mesh        pointer to the mesh
 * @param dofIndexKey the DofManager key for the dof index array
 * @param regions     list of region names to include (if empty, all regions are used)
 * @param numComp     number of components per cell
 * @param sparsity    the matrix to be populated, must be properly sized.
 */
void makeSparsityTPFA( DomainPartition const & domain,
                       string const & dofIndexKey,
                       std::vector< DofManager::FieldSupport > const & support,
                       globalIndex const rankOffset,
                       localIndex const numComp,
                       CRSMatrix< real64 > & sparsity )
{
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & mesh = meshBody.getMeshLevel( regions.meshLevelName );

    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofIndex =
      elemManager.constructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofIndexKey );

    // Make a set of target region indices to check face fluxes.
    SortedArray< localIndex > regionSet;
    elemManager.forElementRegions( regions.regionNames, [&]( localIndex const, ElementRegionBase const & region )
    {
      regionSet.insert( region.getIndexInParent() );
    } );

    // prepare data for assembly loop
    FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
    localIndex const numElem = faceToElem.size( 1 );

    array1d< globalIndex > localDofIndex( numElem * numComp );

    array2d< real64 > localValues( numElem * numComp, numElem * numComp );
    localValues.setValues< serialPolicy >( 1.0 );

    // Loop over faces and assemble TPFA-style "flux" contributions
    for( localIndex kf=0; kf<faceManager.size(); ++kf )
    {
      localIndex const er0 = faceToElem.m_toElementRegion[kf][0];
      localIndex const er1 = faceToElem.m_toElementRegion[kf][1];

      if( er0 >= 0 && er1 >= 0 && regionSet.contains( er0 ) && regionSet.contains( er1 ) )
      {
        localIndex count = 0;
        for( localIndex ke = 0; ke < numElem; ++ke )
        {
          localIndex const er  = faceToElem.m_toElementRegion[kf][ke];
          localIndex const esr = faceToElem.m_toElementSubRegion[kf][ke];
          localIndex const ei  = faceToElem.m_toElementIndex[kf][ke];

          if( er>-1 && esr>-1 && ei>-1 )
          {
            for( localIndex c = 0; c < numComp; ++c )
            {
              localDofIndex[count++] = elemDofIndex[er][esr][ei] + c;
            }
          }
        }
        std::sort( localDofIndex.begin(), localDofIndex.end() );
        for( localIndex row=0; row<count; ++row )
        {
          localIndex const localDofNumber = localDofIndex[row] - rankOffset;
          if( localDofNumber < 0 || localDofNumber >= sparsity.numRows() )
            continue;
          sparsity.insertNonZeros( localDofNumber,
                                   localDofIndex.data(),
                                   localValues.data(),
                                   count );
        }
      }
    }
  }
}

/**
 * @brief Populate a FEM-type sparsity pattern.
 * @param mesh        pointer to the mesh
 * @param dofIndexKey the DofManager key for the dof index array
 * @param regions     list of region names to include (if empty, all regions are used)
 * @param numComp     number of components per node
 * @param sparsity    the matrix to be populated, must be properly sized.
 */
void makeSparsityFEM( DomainPartition const & domain,
                      string const & dofIndexKey,
                      std::vector< DofManager::FieldSupport > const & support,
                      globalIndex const rankOffset,
                      localIndex const numComp,
                      CRSMatrix< real64 > & sparsity )
{
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & mesh = meshBody.getMeshLevel( regions.meshLevelName );

    ElementRegionManager const & elemManager = mesh.getElemManager();
    NodeManager const & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > nodeDofIndex = nodeManager.getReference< array1d< globalIndex > >( dofIndexKey );

    // perform assembly loop over elements
    elemManager.forElementSubRegions( regions.regionNames, [&]( localIndex const, auto const & subRegion )
    {
      using NodeMapType = typename TYPEOFREF( subRegion ) ::NodeMapType;
      traits::ViewTypeConst< NodeMapType > const
      nodeMap = subRegion.template getReference< NodeMapType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );

      localIndex const numNode = subRegion.numNodesPerElement();
      array1d< globalIndex > localDofIndex( numNode * numComp );
      array2d< real64 > localValues( numNode * numComp, numNode * numComp );
      localValues.setValues< serialPolicy >( 1.0 );

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( localIndex a = 0; a < numNode; ++a )
        {
          for( localIndex c = 0; c < numComp; ++c )
          {
            localDofIndex[a * numComp + c] = nodeDofIndex[nodeMap[k][a]] + c;
          }
        }
        std::sort( localDofIndex.begin(), localDofIndex.end() );
        for( localIndex row=0; row<(numNode * numComp); ++row )
        {
          localIndex const localDofNumber = localDofIndex[row] - rankOffset;
          if( localDofNumber < 0 || localDofNumber >= sparsity.numRows() )
            continue;
          sparsity.insertNonZeros( localDofNumber,
                                   localDofIndex.data(),
                                   localValues.data(),
                                   (numNode * numComp) );
        }
      }
    } );
  }
}

/**
 * @brief Populate a FEM/FVM coupling sparsity.
 * @param mesh            pointer to the mesh
 * @param dofIndexKeyNode the DofManager key for the node-based dof index array
 * @param dofIndexKeyElem the DofManager key for the element-based dof index array
 * @param regions         list of region names to include (if empty, all regions are used)
 * @param numCompNode     number of components per node
 * @param numCompElem     number of components per element
 * @param sparsity        the matrix to be populated, must be properly sized.
 */
void makeSparsityFEM_FVM( DomainPartition const & domain,
                          string const & dofIndexKeyNode,
                          string const & dofIndexKeyElem,
                          std::vector< DofManager::FieldSupport > const & support,
                          globalIndex const rankOffset,
                          localIndex const numCompNode,
                          localIndex const numCompElem,
                          CRSMatrix< real64 > & sparsity )
{
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & mesh = meshBody.getMeshLevel( regions.meshLevelName );

    ElementRegionManager const & elemManager = mesh.getElemManager();
    NodeManager const & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > nodeDofIndex =
      nodeManager.getReference< array1d< globalIndex > >( dofIndexKeyNode );

    // perform assembly loop over elements
    elemManager.forElementSubRegions( regions.regionNames, [&]( localIndex const, auto const & subRegion )
    {
      using NodeMapType = typename TYPEOFREF( subRegion ) ::NodeMapType;
      traits::ViewTypeConst< NodeMapType > const
      nodeMap = subRegion.template getReference< NodeMapType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );

      arrayView1d< globalIndex const > elemDofIndex =
        subRegion.template getReference< array1d< globalIndex > >( dofIndexKeyElem );

      localIndex const numNode = subRegion.numNodesPerElement();

      array1d< globalIndex > localNodeDofIndex( numNode * numCompNode );
      array1d< globalIndex > localElemDofIndex( numCompElem );
      array2d< real64 > localValues1( numNode * numCompNode, numCompElem );
      localValues1.setValues< serialPolicy >( 1.0 );
      array2d< real64 > localValues2( numCompElem, numNode * numCompNode );
      localValues2.setValues< serialPolicy >( 1.0 );

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( localIndex c = 0; c < numCompElem; ++c )
        {
          localElemDofIndex[c] = elemDofIndex[k] + c;
        }
        for( localIndex a = 0; a < numNode; ++a )
        {
          for( localIndex c = 0; c < numCompNode; ++c )
          {
            localNodeDofIndex[a * numCompNode + c] = nodeDofIndex[nodeMap[k][a]] + c;
          }
        }

        std::sort( localNodeDofIndex.begin(), localNodeDofIndex.end() );
        std::sort( localElemDofIndex.begin(), localElemDofIndex.end() );

        for( localIndex row=0; row<(numNode * numCompNode); ++row )
        {
          localIndex const localDofNumber = localNodeDofIndex[row] - rankOffset;
          if( localDofNumber < 0 || localDofNumber >= sparsity.numRows() )
            continue;
          sparsity.insertNonZeros( localDofNumber,
                                   localElemDofIndex.data(),
                                   localValues1.data(),
                                   numCompElem );
        }

        for( localIndex row=0; row<(numCompElem); ++row )
        {
          localIndex const localDofNumber = localElemDofIndex[row] - rankOffset;
          if( localDofNumber < 0 || localDofNumber >= sparsity.numRows() )
            continue;
          sparsity.insertNonZeros( localDofNumber,
                                   localNodeDofIndex.data(),
                                   localValues2.data(),
                                   (numNode * numCompNode) );
        }

      }
    } );
  }
}

/**
 * @brief Create a mass matrix-type sparsity pattern (diagonal)
 * @param mesh        pointer to the mesh
 * @param dofIndexKey the DofManager key for the dof index array
 * @param regions     list of region names to include (if empty, all regions are used)
 * @param numComp     number of components per cell
 * @param sparsity    the matrix to be populated
 */
void makeSparsityMass( DomainPartition const & domain,
                       string const & dofIndexKey,
                       std::vector< DofManager::FieldSupport > const & support,
                       globalIndex const rankOffset,
                       localIndex const numComp,
                       CRSMatrix< real64 > & sparsity )
{

  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & mesh = meshBody.getMeshLevel( regions.meshLevelName );

    ElementRegionManager const & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofIndex =
      elemManager.constructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofIndexKey );

    array1d< globalIndex > localDofIndex( numComp );
    array2d< real64 > localValues( numComp, numComp );
    localValues.setValues< serialPolicy >( 1.0 );

    forLocalObjects< FieldLocation::Elem >( mesh, regions.regionNames, [&]( auto const & idx )
    {
      for( localIndex c = 0; c < numComp; ++c )
      {
        localDofIndex[c] = elemDofIndex[idx[0]][idx[1]][idx[2]] + c;
      }

      std::sort( localDofIndex.begin(), localDofIndex.end() );
      for( localIndex row=0; row<numComp; ++row )
      {
        localIndex const localDofNumber = localDofIndex[row] - rankOffset;
        if( localDofNumber < 0 || localDofNumber >= sparsity.numRows() )
          continue;
        sparsity.insertNonZeros( localDofNumber,
                                 localDofIndex.data(),
                                 localValues.data(),
                                 numComp );
      }

    } );
  }
}

/**
 * @brief Create a flux sparsity pattern (face-based dofs coupled via elements)
 * @param mesh        pointer to the mesh
 * @param dofIndexKey the DofManager key for the dof index array
 * @param regions     list of region names to include (if empty, all regions are used)
 * @param numComp     number of components per cell
 * @param sparsity    the matrix to be populated
 */
void makeSparsityFlux( DomainPartition const & domain,
                       string const & dofIndexKey,
                       std::vector< DofManager::FieldSupport > const & support,
                       globalIndex const rankOffset,
                       localIndex const numComp,
                       CRSMatrix< real64 > & sparsity )
{
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & mesh = meshBody.getMeshLevel( regions.meshLevelName );
    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    arrayView1d< globalIndex const > faceDofIndex =
      faceManager.getReference< array1d< globalIndex > >( dofIndexKey );

    // perform assembly loop over elements
    elemManager.forElementSubRegions( regions.regionNames, [&]( localIndex const, auto const & subRegion )
    {
      using FaceMapType = typename TYPEOFREF( subRegion ) ::FaceMapType;
      traits::ViewTypeConst< FaceMapType > const
      faceMap = subRegion.template getReference< FaceMapType >( ElementSubRegionBase::viewKeyStruct::faceListString() );

      localIndex const numFace = subRegion.numFacesPerElement();
      array1d< globalIndex > localDofIndex( numFace * numComp );
      array2d< real64 > localValues( numFace * numComp, numFace * numComp );
      localValues.setValues< serialPolicy >( 1.0 );

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( localIndex a = 0; a < numFace; ++a )
        {
          for( localIndex c = 0; c < numComp; ++c )
          {
            localDofIndex[a * numComp + c] = faceDofIndex[faceMap[k][a]] + c;
          }
        }

        std::sort( localDofIndex.begin(), localDofIndex.end() );

        for( localIndex row=0; row<(numFace*numComp); ++row )
        {
          localIndex const localDofNumber = localDofIndex[row] - rankOffset;
          if( localDofNumber < 0 || localDofNumber >= sparsity.numRows() )
            continue;
          sparsity.insertNonZeros( localDofNumber,
                                   localDofIndex.data(),
                                   localValues.data(),
                                   (numFace*numComp) );
        }

      }
    } );
  }
}

} // namespace testing
} // namespace geos

#endif //GEOS_LINEARALGEBRA_UNITTESTS_TESTDOFMANAGERUTILS_HPP_
