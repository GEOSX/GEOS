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
 * @file BlueprintOutput.cpp
 */

/// Source includes
#include "BlueprintOutput.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"
#include "dataRepository/ConduitRestart.hpp"

// TPL includes
#include <conduit.hpp>
#include <conduit_blueprint.hpp>
#include <conduit_relay.hpp>

namespace geos
{
namespace internal
{

/**
 * @brief @return The Blueprint shape from the GEOSX element type string.
 * @param elementType the elementType to look up.
 */
string toBlueprintShape( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Tetrahedron: return "tet";
    case ElementType::Hexahedron: return "hex";
    default:
    {
      GEOS_ERROR( "No Blueprint type for element type: " << elementType );
      return {};
    }
  }
}

static std::vector< int > getBlueprintNodeOrdering( ElementType const elementType )
{
  // Same as VTK, but kept separate for flexibility
  switch( elementType )
  {
    case ElementType::Vertex:        return { 0 };
    case ElementType::Line:          return { 0, 1 };
    case ElementType::Triangle:      return { 0, 1, 2 };
    case ElementType::Quadrilateral: return { 0, 1, 2, 3 }; // TODO check
    case ElementType::Polygon:       return { }; // TODO
    case ElementType::Tetrahedron:   return { 1, 0, 2, 3 };
    case ElementType::Pyramid:       return { 0, 3, 2, 1, 4, 0, 0, 0 };
    case ElementType::Wedge:         return { 0, 4, 2, 1, 5, 3, 0, 0 };
    case ElementType::Hexahedron:    return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ElementType::Prism5:        return { }; // TODO
    case ElementType::Prism6:        return { }; // TODO
    case ElementType::Prism7:        return { }; // TODO
    case ElementType::Prism8:        return { }; // TODO
    case ElementType::Prism9:        return { }; // TODO
    case ElementType::Prism10:       return { }; // TODO
    case ElementType::Prism11:       return { }; // TODO
    case ElementType::Polyhedron:    return { }; // TODO
  }
  return {};
}

/**
 * @brief Outputs the element to node map of @p subRegion to @p connectivity in VTK order.
 * @param subRegion The sub-region to output.
 * @param connectivity The Conduit Node to output to.
 */
void reorderElementToNodeMap( CellElementSubRegion const & subRegion, conduit::Node & connectivity )
{
  GEOS_MARK_FUNCTION;

  arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodeMap = subRegion.nodeList();
  localIndex const numElems = elemToNodeMap.size( 0 );
  localIndex const numNodesPerElem = elemToNodeMap.size( 1 );

  std::vector< int > const vtkOrdering = getBlueprintNodeOrdering( subRegion.getElementType() );
  GEOS_ERROR_IF_NE( localIndex( vtkOrdering.size() ), numNodesPerElem );

  constexpr int conduitTypeID = dataRepository::conduitTypeInfo< localIndex >::id;
  conduit::DataType const dtype( conduitTypeID, elemToNodeMap.size() );
  connectivity.set( dtype );

  localIndex * const reorderedConnectivity = connectivity.value();
  forAll< serialPolicy >( numElems, [reorderedConnectivity, numNodesPerElem, elemToNodeMap, &vtkOrdering] ( localIndex const i )
  {
    for( localIndex j = 0; j < numNodesPerElem; ++j )
    {
      reorderedConnectivity[ i * numNodesPerElem + j ] = elemToNodeMap( i, vtkOrdering[ j ] );
    }
  } );
}

} /// namespace internal;

///////////////////////////////////////////////////////////////////////////////////////////////////
BlueprintOutput::BlueprintOutput( string const & name,
                                  dataRepository::Group * const parent ):
  OutputBase( name, parent )
{
  registerWrapper( "plotLevel", &m_plotLevel ).
    setApplyDefaultValue( dataRepository::PlotLevel::LEVEL_1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setRTTypeName( rtTypes::CustomTypes::plotLevel ).
    setDescription( "Determines which fields to write." );

  registerWrapper( "outputFullQuadratureData", &m_outputFullQuadratureData ).
    setApplyDefaultValue( false ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "If true writes out data associated with every quadrature point." );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
bool BlueprintOutput::execute( real64 const time,
                               real64 const,
                               integer const cycle,
                               integer const,
                               real64 const,
                               DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getBaseDiscretization();

  conduit::Node meshRoot;
  conduit::Node & mesh = meshRoot[ "mesh" ];
  conduit::Node & coordset = mesh[ "coordsets/nodes" ];
  conduit::Node & topologies = mesh[ "topologies" ];

  mesh[ "state/time" ] = time;
  mesh[ "state/cycle" ] = cycle;

  addNodalData( meshLevel.getNodeManager(), coordset, topologies, mesh[ "fields" ] );

  dataRepository::Group averagedElementData( "averagedElementData", this );
  addElementData( meshLevel.getElemManager(), coordset, topologies, mesh[ "fields" ], averagedElementData );

  /// The Blueprint will complain if the fields node is present but empty.
  if( mesh[ "fields" ].number_of_children() == 0 )
  {
    mesh.remove( "fields" );
  }

  /// Verify that the mesh conforms to the Blueprint.
  conduit::Node info;
  GEOS_ASSERT_MSG( conduit::blueprint::verify( "mesh", meshRoot, info ), info.to_json() );

  /// Generate the Blueprint index.
  conduit::Node fileRoot;
  conduit::Node & index = fileRoot[ "blueprint_index/mesh" ];
  conduit::blueprint::mesh::generate_index( mesh, "mesh", MpiWrapper::commSize(), index );

  /// Verify that the index conforms to the Blueprint.
  info.reset();
  GEOS_ASSERT_MSG( conduit::blueprint::mesh::index::verify( index, info ), info.to_json() );

  /// Write out the root index file, then write out the mesh.
  string const completePath = GEOS_FMT( "{}/blueprintFiles/cycle_{:07}", OutputBase::getOutputDirectory(), cycle );
  string const filePathForRank = dataRepository::writeRootFile( fileRoot, completePath );
  conduit::relay::io::save( meshRoot, filePathForRank, "hdf5" );

  return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void BlueprintOutput::addNodalData( NodeManager const & nodeManager,
                                    conduit::Node & coordset,
                                    conduit::Node & topologies,
                                    conduit::Node & fields )
{
  GEOS_MARK_FUNCTION;

  /// Populate the coordset group
  coordset[ "type" ] = "explicit";
  dataRepository::wrapperHelpers::populateMCArray( nodeManager.referencePosition(),
                                                   coordset[ "values" ],
                                                   { "x", "y", "z" } );

  /// Create the points topology
  string const coordsetName = coordset.name();
  conduit::Node & nodeTopology = topologies[ coordsetName ];
  nodeTopology[ "coordset" ] = coordsetName;

  /// TODO: Once VisIT supports the implicit "points" topology we can just do the following.
  /// See https://github.com/visit-dav/visit/issues/4593
  // nodeTopology[ "type" ] = "points";

  nodeTopology[ "type" ] = "unstructured";
  nodeTopology[ "elements/shape" ] = "point";
  conduit::Node & connectivity = nodeTopology[ "elements/connectivity" ];

  localIndex const numNodes = nodeManager.size();
  constexpr int conduitTypeID = dataRepository::conduitTypeInfo< localIndex >::id;
  conduit::DataType const dtype( conduitTypeID, numNodes );
  connectivity.set( dtype );

  localIndex * const nodeIDs = connectivity.value();
  forAll< serialPolicy >( numNodes, [nodeIDs] ( localIndex const i )
  {
    nodeIDs[ i ] = i;
  } );

  /// Write out the fields.
  writeOutWrappersAsFields( nodeManager, fields, coordsetName );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void BlueprintOutput::addElementData( ElementRegionManager const & elemRegionManager,
                                      conduit::Node & coordset,
                                      conduit::Node & topologies,
                                      conduit::Node & fields,
                                      dataRepository::Group & averagedElementData )
{
  GEOS_MARK_FUNCTION;

  elemRegionManager.forElementSubRegionsComplete< CellElementSubRegion >(
    [&] ( localIndex, localIndex, ElementRegionBase const & region, CellElementSubRegion const & subRegion )
  {
    string const topologyName = region.getName() + "-" + subRegion.getName();

    /// Create the topology representing the sub-region.
    conduit::Node & topology = topologies[ topologyName ];
    topology[ "coordset" ] = coordset.name();
    topology[ "type" ] = "unstructured";
    topology[ "elements/shape" ] = internal::toBlueprintShape( subRegion.getElementType() );
    internal::reorderElementToNodeMap( subRegion, topology[ "elements/connectivity" ] );

    /// Write out the fields.
    writeOutWrappersAsFields( subRegion, fields, topologyName );

    /// Write out the quadrature averaged constitutive data and the full data if requested.
    Group & averagedSubRegionData = averagedElementData.registerGroup( topologyName );
    subRegion.getConstitutiveModels().forSubGroups( [&]( dataRepository::Group const & constitutiveModel )
    {
      writeOutConstitutiveData( constitutiveModel, fields, topologyName, averagedSubRegionData );

      if( m_outputFullQuadratureData )
      {
        writeOutWrappersAsFields( constitutiveModel, fields, topologyName, constitutiveModel.getName() );
      }
    } );
  } );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void BlueprintOutput::writeOutWrappersAsFields( Group const & group,
                                                conduit::Node & fields,
                                                string const & topology,
                                                string const & prefix )
{
  GEOS_MARK_FUNCTION;

  group.forWrappers( [&] ( dataRepository::WrapperBase const & wrapper )
  {
    if( wrapper.getPlotLevel() <= m_plotLevel && wrapper.sizedFromParent() )
    {
      string const name = prefix.empty() ? wrapper.getName() : prefix + "-" + wrapper.getName();

      // conduit::Node & field = fields[ name ];
      // field[ "association" ] = "element";
      // field[ "volume_dependent" ] = "false";
      // field[ "topology" ] = topology;
      // wrapper.populateMCArray( field[ "values" ] );

      /// TODO: Replace with code above once https://github.com/visit-dav/visit/issues/4637 is fixed and released.
      wrapper.addBlueprintField( fields, name, topology );
    }
  } );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void BlueprintOutput::writeOutConstitutiveData( dataRepository::Group const & constitutiveModel,
                                                conduit::Node & fields,
                                                string const & topology,
                                                dataRepository::Group & averagedSubRegionData )
{
  GEOS_MARK_FUNCTION;

  Group & averagedConstitutiveData = averagedSubRegionData.registerGroup( constitutiveModel.getName() );

  constitutiveModel.forWrappers( [&] ( dataRepository::WrapperBase const & wrapper )
  {
    if( wrapper.getPlotLevel() <= m_plotLevel && wrapper.sizedFromParent() )
    {
      string const fieldName = constitutiveModel.getName() + "-quadrature-averaged-" + wrapper.getName();
      averagedConstitutiveData.registerWrapper( wrapper.averageOverSecondDim( fieldName, averagedConstitutiveData ) ).
        addBlueprintField( fields, fieldName, topology );
    }
  } );
}



REGISTER_CATALOG_ENTRY( OutputBase, BlueprintOutput, string const &, dataRepository::Group * const )

} /* namespace geos */
