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

// Source includes
#include "VTKPolyDataWriterInterface.hpp"

#include "common/logger/Logger.hpp"
#include "common/TypeDispatch.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/DomainPartition.hpp"
#include "fileIO/Outputs/OutputUtilities.hpp"

// TPL includes
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPassThrough.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

// System includes
#include <numeric>
#include <unordered_set>

namespace geos
{

using namespace dataRepository;

namespace vtk
{

VTKPolyDataWriterInterface::VTKPolyDataWriterInterface( string name ):
  m_outputDir( "." ),
  m_outputName( std::move( name ) ),
  m_pvd( m_outputName + ".pvd" ),
  m_writeGhostCells( false ),
  m_plotLevel( PlotLevel::LEVEL_1 ),
  m_requireFieldRegistrationCheck( true ),
  m_previousCycle( -1 ),
  m_outputMode( VTKOutputMode::BINARY ),
  m_outputRegionType( VTKRegionTypes::ALL ),
  m_writeFaceElementsAs3D( false )
{}

static int
toVTKCellType( ElementType const elementType, localIndex const numNodes )
{
  switch( elementType )
  {
    case ElementType::Vertex:        return VTK_VERTEX;
    case ElementType::Line:          return VTK_LINE;
    case ElementType::Triangle:      return VTK_TRIANGLE;
    case ElementType::Quadrilateral: return VTK_QUAD;
    case ElementType::Polygon:       return VTK_POLYGON;
    case ElementType::Tetrahedron:   return VTK_TETRA;
    case ElementType::Pyramid:       return VTK_PYRAMID;
    case ElementType::Wedge:         return VTK_WEDGE;
    case ElementType::Hexahedron:
      switch( numNodes )
      {
        case 8:
          return VTK_HEXAHEDRON;
        case 27:
          return VTK_QUADRATIC_HEXAHEDRON;
        default:
          return VTK_HEXAHEDRON;
      }
    case ElementType::Prism5:        return VTK_PENTAGONAL_PRISM;
    case ElementType::Prism6:        return VTK_HEXAGONAL_PRISM;
    case ElementType::Prism7:        return VTK_POLYHEDRON;
    case ElementType::Prism8:        return VTK_POLYHEDRON;
    case ElementType::Prism9:        return VTK_POLYHEDRON;
    case ElementType::Prism10:       return VTK_POLYHEDRON;
    case ElementType::Prism11:       return VTK_POLYHEDRON;
    case ElementType::Polyhedron:    return VTK_POLYHEDRON;
  }
  return VTK_EMPTY_CELL;
}

static int
toVTKCellType( ParticleType const particleType )
{
  switch( particleType )
  {
    case ParticleType::SinglePoint:   return VTK_HEXAHEDRON;
    case ParticleType::CPDI:          return VTK_HEXAHEDRON;
    case ParticleType::CPDI2:         return VTK_HEXAHEDRON;
    case ParticleType::CPTI:          return VTK_TETRA;
  }
  return VTK_EMPTY_CELL;
}

static std::vector< int >
getVtkToGeosxNodeOrdering( ParticleType const particleType )
{
  switch( particleType )
  {
    case ParticleType::SinglePoint:   return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ParticleType::CPDI:          return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ParticleType::CPDI2:         return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ParticleType::CPTI:          return { 0, 1, 2, 3 };
  }
  return {};
}

/**
 * @brief Provide the local list of nodes or face streams for the corresponding VTK element
 *
 * @param elementType geos element type
 * @return list of nodes or face streams
 *
 * For geos element with existing standard VTK element the corresponding list of nodes is provided.
 * For Prism7+, the geos element is converted to VTK_POLYHEDRON. The vtkUnstructuredGrid
 * stores polyhedron cells as face streams of the following format:
 * [numberOfCellFaces,
 * (numberOfPointsOfFace0, pointId0, pointId1, ... ),
 * (numberOfPointsOfFace1, pointId0, pointId1, ... ),
 * ...]
 * We use the same format except that the number of faces and the number of nodes per faces
 * are provided as negative values. This convention provides a simple way to isolate the local
 * nodes for mapping purpose while keeping a face streams data structure. The negative values are
 * converted to positives when generating the VTK_POLYHEDRON. Check getVtkCells() for more details.
 */
static std::vector< int > getVtkConnectivity( ElementType const elementType, localIndex const numNodes )
{
  switch( elementType )
  {
    case ElementType::Vertex:        return { 0 };
    case ElementType::Line:          return { 0, 1 };
    case ElementType::Triangle:      return { 0, 1, 2 };
    case ElementType::Quadrilateral: return { 0, 1, 2, 3 };
    case ElementType::Polygon:       return { };  // TODO
    case ElementType::Tetrahedron:   return { 0, 1, 2, 3 };
    case ElementType::Pyramid:       return { 0, 1, 3, 2, 4 };
    case ElementType::Wedge:         return { 0, 4, 2, 1, 5, 3 };
    case ElementType::Hexahedron:
      switch( numNodes )
      {
        case 8:
          return { 0, 1, 3, 2, 4, 5, 7, 6 };
          break;
        case 27:
          // Numbering convention changed between VTK writer 1.0 and 2.2, see
          // https://discourse.julialang.org/t/writevtk-node-numbering-for-27-node-lagrange-hexahedron/93698/8
          // GEOS uses 1.0 API
          return { 0, 2, 8, 6, 18, 20, 26, 24, 1, 5, 7, 3, 19, 23, 25, 21, 9, 11, 17, 15, 12, 14, 10, 16, 4, 22, 17 };
          break;
        default:
          return { }; // TODO
      }
    case ElementType::Prism5:        return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    case ElementType::Prism6:        return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    case ElementType::Prism7:        return {-9,
                                             -7, 0, 1, 2, 3, 4, 5, 6,
                                             -7, 7, 8, 9, 10, 11, 12, 13,
                                             -4, 0, 1, 8, 7,
                                             -4, 1, 2, 9, 8,
                                             -4, 2, 3, 10, 9,
                                             -4, 3, 4, 11, 10,
                                             -4, 4, 5, 12, 11,
                                             -4, 5, 6, 13, 12,
                                             -4, 6, 0, 7, 13 };
    case ElementType::Prism8:        return {-10,
                                             -8, 0, 1, 2, 3, 4, 5, 6, 7,
                                             -8, 8, 9, 10, 11, 12, 13, 14, 15,
                                             -4, 0, 1, 9, 8,
                                             -4, 1, 2, 10, 9,
                                             -4, 2, 3, 11, 10,
                                             -4, 3, 4, 12, 11,
                                             -4, 4, 5, 13, 12,
                                             -4, 5, 6, 14, 13,
                                             -4, 6, 7, 15, 14,
                                             -4, 7, 0, 8, 15 };
    case ElementType::Prism9:        return {-11,
                                             -9, 0, 1, 2, 3, 4, 5, 6, 7, 8,
                                             -9, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                                             -4, 0, 1, 10, 9,
                                             -4, 1, 2, 11, 10,
                                             -4, 2, 3, 12, 11,
                                             -4, 3, 4, 13, 12,
                                             -4, 4, 5, 14, 13,
                                             -4, 5, 6, 15, 14,
                                             -4, 6, 7, 16, 15,
                                             -4, 7, 8, 17, 16,
                                             -4, 8, 0, 9, 17 };
    case ElementType::Prism10:       return {-12,
                                             -10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                             -10, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                             -4, 0, 1, 11, 10,
                                             -4, 1, 2, 12, 11,
                                             -4, 2, 3, 13, 12,
                                             -4, 3, 4, 14, 13,
                                             -4, 4, 5, 15, 14,
                                             -4, 5, 6, 16, 15,
                                             -4, 6, 7, 17, 16,
                                             -4, 7, 8, 18, 17,
                                             -4, 8, 9, 19, 18,
                                             -4, 9, 0, 10, 19 };
    case ElementType::Prism11:       return {-13,
                                             -11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                             -11, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                                             -4, 0, 1, 12, 11,
                                             -4, 1, 2, 13, 12,
                                             -4, 2, 3, 14, 13,
                                             -4, 3, 4, 15, 14,
                                             -4, 4, 5, 16, 15,
                                             -4, 5, 6, 17, 16,
                                             -4, 6, 7, 18, 17,
                                             -4, 7, 8, 19, 18,
                                             -4, 8, 9, 20, 19,
                                             -4, 9, 10, 21, 20,
                                             -4, 10, 0, 11, 21 };
    case ElementType::Polyhedron:    return { };  // TODO
  }
  return {};
}

/**
 * @brief Gets the vertices coordinates as a VTK Object for @p nodeManager
 * @param[in] nodeManager the NodeManager associated with the domain being written
 * @return a VTK object storing all nodes of the mesh
 */
static vtkSmartPointer< vtkPoints >
getVtkPoints( NodeManager const & nodeManager,
              arrayView1d< localIndex const > const & nodeIndices )
{
  localIndex const numNodes = LvArray::integerConversion< localIndex >( nodeIndices.size() );
  auto points = vtkSmartPointer< vtkPoints >::New();
  points->SetNumberOfPoints( numNodes );
  auto const coord = nodeManager.referencePosition().toViewConst();
  forAll< parallelHostPolicy >( numNodes, [=, pts = points.GetPointer()]( localIndex const k )
  {
    localIndex const v = nodeIndices[k];
    pts->SetPoint( k, coord[v][0], coord[v][1], coord[v][2] );
  } );
  return points;
}

/**
 * @brief Gets the vertex coordinates as a VTK Object for @p particleManager
 * @param[in] particleManager the ParticleManager associated with the domain being written
 * @return a VTK object storing all particle centers/corners of the mesh
 */
static vtkSmartPointer< vtkPoints >
getVtkPoints( ParticleRegion const & particleRegion ) // TODO: Loop over the subregions owned by this region and operate on them directly
{
  // Particles are plotted as polyhedron with the geometry determined by the particle
  // type.  CPDI particles are parallelepiped (8 corners and 6 faces).
  // TODO: add support for CPTI (tet) and single point (cube or sphere) geometries

  localIndex const numCornersPerParticle = 8; // Each CPDI particle has 8 corners. TODO: add support for other particle types.
  localIndex const numCorners = numCornersPerParticle * particleRegion.getNumberOfParticles();
  auto points = vtkSmartPointer< vtkPoints >::New();
  points->SetNumberOfPoints( numCorners );
  array2d< real64 > const coord = particleRegion.getParticleCorners();
  forAll< parallelHostPolicy >( numCorners, [=, pts = points.GetPointer()]( localIndex const k )
  {
    pts->SetPoint( k, coord[k][0], coord[k][1], coord[k][2] );
  } );
  return points;
}

struct ElementData
{
  std::vector< int > cellTypes;
  vtkSmartPointer< vtkCellArray > cells;
  vtkSmartPointer< vtkPoints > points;
};

/**
 * @brief Gets the cell connectivities and the vertices coordinates as VTK objects for a specific WellElementSubRegion.
 * @param[in] subRegion the WellElementSubRegion to be output
 * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
 * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
 * and a VTKCellArray (with the cell connectivities).
 */
static ElementData
getWell( WellElementSubRegion const & subRegion,
         NodeManager const & nodeManager )
{
  // some notes about WellElementSubRegion:
  // - if the well represented by this subRegion is not on this rank, esr.size() = 0
  // - otherwise, esr.size() is equal to the number of well elements of the well on this rank
  // Each well element has two nodes, shared with the previous and next well elements, respectively
  auto points = vtkSmartPointer< vtkPoints >::New();
  // if esr.size() == 0, we set the number of points and cells to zero
  // if not, we set the number of points to esr.size()+1 and the number of cells to esr.size()
  localIndex const numPoints = subRegion.size() > 0 ? subRegion.size() + 1 : 0;
  points->SetNumberOfPoints( numPoints );
  auto cellsArray = vtkSmartPointer< vtkCellArray >::New();
  cellsArray->SetNumberOfCells( subRegion.size() );
  localIndex const numberOfNodesPerElement = subRegion.numNodesPerElement();
  GEOS_ERROR_IF_NE( numberOfNodesPerElement, 2 );
  std::vector< vtkIdType > connectivity( numberOfNodesPerElement );

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition = nodeManager.referencePosition();

  // note that if esr.size() == 0, we don't have any point or cell to add below and we just return

  for( localIndex edge = 0; edge < subRegion.size(); edge++ )
  {
    localIndex const firstPoint = subRegion.nodeList()[edge][0];
    auto point = referencePosition[firstPoint];
    points->SetPoint( edge, point[0], point[1], point[2] );
    connectivity[0] = edge;
    connectivity[1] = edge + 1; // We can do that because of the pattern in which the wells are stored
    cellsArray->InsertNextCell( numberOfNodesPerElement, connectivity.data() );
  }

  if( subRegion.size() > 0 )
  {
    localIndex const lastPoint = subRegion.nodeList()[subRegion.size() - 1][1];
    auto point = referencePosition[lastPoint];
    points->SetPoint( subRegion.size(), point[0], point[1], point[2] );
  }

  std::vector< int > cellTypes( subRegion.size(), VTK_LINE );
  return { cellTypes, cellsArray, points };
}

/**
 * @brief Gets the cell connectivities and the vertices coordinates as VTK objects for a specific FaceElementSubRegion.
 * @param[in] subRegion the FaceElementSubRegion to be output
 * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
 * @param[in] faceManager the faceManager associated with the DomainPartition being written.
 * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
 * and a VTKCellArray (with the cell connectivities).
 */
static ElementData
getSurface( FaceElementSubRegion const & subRegion,
            NodeManager const & nodeManager,
            FaceManager const & faceManager,
            bool const writeFaceElementsAs3D )
{
  // Get unique node set composing the surface
  auto & elemToFaces = subRegion.faceList();
  auto & elemToNodes = subRegion.nodeList();
  auto & faceToNodes = faceManager.nodeList();

  auto cellArray = vtkSmartPointer< vtkCellArray >::New();
  cellArray->SetNumberOfCells( subRegion.size() );
  std::vector< int > cellTypes;
  cellTypes.reserve( subRegion.size() );

  std::unordered_map< localIndex, localIndex > geos2VTKIndexing;
  geos2VTKIndexing.reserve( subRegion.size() * subRegion.numNodesPerElement() );
  localIndex nodeIndexInVTK = 0;
  // FaceElementSubRegion being heterogeneous, the size of the connectivity vector may vary for each element.
  // In order not to allocate a new vector every time, we combine the usage of `clear` and `push_back`.
  std::vector< vtkIdType > connectivity;

  for( localIndex ei = 0; ei < subRegion.size(); ei++ )
  {
    // we use the nodes of face 0
    auto const & nodes = !writeFaceElementsAs3D ? faceToNodes[elemToFaces( ei, 0 )] : elemToNodes[ei];
    auto const numNodes = nodes.size();

    ElementType const elementType = subRegion.getElementType( ei );
    std::vector< int > vtkOrdering;
    if( elementType == ElementType::Polygon || writeFaceElementsAs3D )
    {
      vtkOrdering.resize( numNodes );
      std::iota( vtkOrdering.begin(), vtkOrdering.end(), 0 );
    }
    else
    {
      vtkOrdering = getVtkConnectivity( elementType, numNodes );
    }

    connectivity.clear();
    for( int const & ordering : vtkOrdering )
    {
      auto const & VTKIndexPos = geos2VTKIndexing.find( nodes[ordering] );
      if( VTKIndexPos == geos2VTKIndexing.end() )
      {
        /// If the node is not found in the geos2VTKIndexing map:
        /// 1. we assign the current value of nodeIndexInVTK to this node in the map (geos2VTKIndexing[nodes[ordering]] =
        /// nodeIndexInVTK++).
        /// 2. we increment nodeIndexInVTK to ensure the next new node gets a unique index.
        /// 3. we add this new VTK node index to the connectivity vector (connectivity.push_back).
        connectivity.push_back( geos2VTKIndexing[nodes[ordering]] = nodeIndexInVTK++ );
      }
      else
      {
        connectivity.push_back( VTKIndexPos->second );
      }
    }

    cellArray->InsertNextCell( vtkOrdering.size(), connectivity.data() );
    cellTypes.emplace_back( toVTKCellType( elementType, numNodes ) );
  }

  auto points = vtkSmartPointer< vtkPoints >::New();
  points->SetNumberOfPoints( geos2VTKIndexing.size() );
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const referencePosition = nodeManager.referencePosition();

  for( auto nodeIndex: geos2VTKIndexing )
  {
    auto point = referencePosition[nodeIndex.first];
    points->SetPoint( nodeIndex.second, point[0], point[1], point[2] );
  }

  return { cellTypes, cellArray, points };
}

/**
 * @brief Gets the cell connectivities and the vertices coordinates as VTK objects for a specific
 * EmbeddedSurafaceSubRegion.
 * @param[in] subRegion the EmbeddedSurfaceSubRegion to be output
 * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
 * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
 * and a VTKCellArray (with the cell connectivities).
 */
static ElementData
getEmbeddedSurface( EmbeddedSurfaceSubRegion const & subRegion,
                    EmbeddedSurfaceNodeManager const & nodeManager )
{
  auto cellsArray = vtkSmartPointer< vtkCellArray >::New();
  auto points = vtkSmartPointer< vtkPoints >::New();

  localIndex const numNodes = nodeManager.size();
  auto const intersectionPoints = nodeManager.referencePosition();

  points->SetNumberOfPoints( numNodes );
  for( localIndex pointIndex = 0; pointIndex < numNodes; ++pointIndex )
  {
    auto const pointCoords = intersectionPoints[pointIndex];
    points->SetPoint( pointIndex, pointCoords[0], pointCoords[1], pointCoords[2] );
  }

  auto const toNodesMap = subRegion.nodeList().toViewConst();
  std::vector< vtkIdType > connectivity( 10 );
  for( localIndex cellIndex = 0; cellIndex < subRegion.size(); ++cellIndex )
  {
    auto const nodes = toNodesMap[cellIndex];
    connectivity.resize( nodes.size() );
    for( localIndex i = 0; i < nodes.size(); ++i )
    {
      connectivity[i] = nodes[i];
    }
    cellsArray->InsertNextCell( nodes.size(), connectivity.data() );
  }

  std::vector< int > cellTypes( subRegion.size(), VTK_POLYGON );
  return { cellTypes, cellsArray, points };
}

struct CellData
{
  std::vector< int > cellTypes;
  vtkSmartPointer< vtkCellArray > cells;
  array1d< localIndex > nodes;
};

/**
 * @brief Gets the cell connectivities as a VTK object for the CellElementRegion @p region
 * @param[in] region the CellElementRegion to be written
 * @param[in] numNodes number of local nodes
 * @return a struct consisting of:
 *         - a list of types for each cell,
 *         - a VTK object containing the connectivity information
 *         - a list of relevant node indices in order in which they must be stored
 */
static CellData
getVtkCells( CellElementRegion const & region,
             localIndex const numNodes )
{
  localIndex const numElems = region.getNumberOfElements< CellElementSubRegion >();
  if( numElems == 0 )
  {
    return { {}, vtkSmartPointer< vtkCellArray >::New(), {} };
  }

  // 1. Mark (in parallel) relevant nodes
  std::vector< localIndex > newNodeIndices( numNodes ); // temporary use as a marker array
  region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    auto const nodeList = subRegion.nodeList().toViewConst();
    forAll< parallelHostPolicy >( subRegion.size(), [&, nodeList]( localIndex const cellIdx )
    {
      auto const nodes = nodeList[cellIdx];
      for( localIndex i = 0; i < nodes.size(); ++i )
      {
        // use atomic write to avoid technical UB
        RAJA::atomicExchange< parallelHostAtomic >( &newNodeIndices[nodes[i]], 1 );
      }
    } );
  } );

  // 2. Assign new node IDs (serial step)
  array1d< localIndex > relevantNodes;
  relevantNodes.reserve( numNodes );
  localIndex newNodeIdx = 0;
  for( localIndex nodeIdx = 0; nodeIdx < numNodes; ++nodeIdx )
  {
    if( newNodeIndices[nodeIdx] > 0 )
    {
      relevantNodes.emplace_back( nodeIdx );
      newNodeIndices[nodeIdx] = newNodeIdx++;
    }
  }

  // 3. Write connectivity using new node IDs
  localIndex const numConns = [&]
  {
    localIndex numConn = 0;
    region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
    {
      numConn += subRegion.size() * getVtkConnectivity( subRegion.getElementType(), subRegion.nodeList().size( 1 ) ).size();
    } );
    return numConn;
  }();

  std::vector< int > cellTypes;
  cellTypes.reserve( numElems );

  auto const offsets = vtkSmartPointer< vtkIdTypeArray >::New();
  offsets->SetNumberOfTuples( numElems + 1 );

  auto const connectivity = vtkSmartPointer< vtkIdTypeArray >::New();
  connectivity->SetNumberOfTuples( numConns );

  // 4. Write connectivity using new node IDs
  localIndex elemOffset = 0;
  localIndex connOffset = 0;
  region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    auto const nodeList = subRegion.nodeList().toViewConst();
    auto subRegionNumNodes = nodeList.size( 1 );
    cellTypes.insert( cellTypes.end(), subRegion.size(), toVTKCellType( subRegion.getElementType(), subRegionNumNodes ) );
    std::vector< int > const vtkOrdering = getVtkConnectivity( subRegion.getElementType(), subRegionNumNodes );
    localIndex const numVtkData = vtkOrdering.size();

    // For all geos element, the corresponding VTK data are copied in "connectivity".
    // Local nodes are mapped to global indices. Any negative value in "vtkOrdering"
    // corresponds to the number of faces or the number of nodes per faces, and they
    // are copied as positive values.
    // Here we privilege code simplicity. This can be more efficient (less tests) if the code is
    // specialized for each type of subregion.
    // This is not a time sensitive part of the code. Can be optimized later if needed.
    forAll< parallelHostPolicy >( subRegion.size(), [&]( localIndex const c )
    {
      localIndex const elemConnOffset = connOffset + c * numVtkData;
      auto const nodes = nodeList[c];
      for( localIndex i = 0; i < numVtkData; ++i )
      {
        if( vtkOrdering[i] < 0 )
        {
          connectivity->SetTypedComponent( elemConnOffset + i, 0, -vtkOrdering[i] );
        }
        else
        {
          connectivity->SetTypedComponent( elemConnOffset + i, 0, newNodeIndices[nodes[vtkOrdering[i]]] );
        }
      }
      offsets->SetTypedComponent( elemOffset + c, 0, elemConnOffset );
    } );

    elemOffset += subRegion.size();
    connOffset += subRegion.size() * numVtkData;
  } );
  offsets->SetTypedComponent( elemOffset, 0, connOffset );

  auto cellsArray = vtkSmartPointer< vtkCellArray >::New();
  cellsArray->SetData( offsets, connectivity );

  return { std::move( cellTypes ), cellsArray, std::move( relevantNodes ) };
}

using ParticleData = std::pair< std::vector< int >, vtkSmartPointer< vtkCellArray > >;
/**
 * @brief Gets the cell connectivities as a VTK object for the ParticleRegion @p region
 * @param[in] region the ParticleRegion to be written
 * @return a standard pair consisting of:
 *         - a list of types for each cell,
 *         - a VTK object containing the connectivity information
 */
static ParticleData
getVtkCells( ParticleRegion const & region )
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( region.getNumberOfParticles< ParticleRegion >() );
  std::vector< int > cellType;
  cellType.reserve( region.getNumberOfParticles< ParticleRegion >() );

  vtkIdType nodeIndex = 0;

  region.forParticleSubRegions< ParticleSubRegion >( [&]( ParticleSubRegion const & subRegion )
  {
    std::vector< int > vtkOrdering = getVtkToGeosxNodeOrdering( subRegion.getParticleType() );
    std::vector< vtkIdType > connectivity( vtkOrdering.size() );
    int vtkCellType = toVTKCellType( subRegion.getParticleType() );
    for( localIndex c = 0; c < subRegion.size(); c++ )
    {
      for( std::size_t i = 0; i < connectivity.size(); i++ )
      {

        connectivity[i] = vtkOrdering.size()*nodeIndex + vtkOrdering[i];
      }
      nodeIndex++;
      cellType.push_back( vtkCellType );
      cellsArray->InsertNextCell( vtkOrdering.size(), connectivity.data() );
    }
  } );
  return std::make_pair( cellType, cellsArray );
}

/**
 * @brief Writes timestamp information required by VisIt
 * @param[in] ug the VTK unstructured grid.
 * @param[in] time the current time-step
 */
static void
writeTimestamp( vtkUnstructuredGrid * ug,
                real64 const time )
{
  auto t = vtkSmartPointer< vtkDoubleArray >::New();
  t->SetName( "TIME" );
  t->SetNumberOfTuples( 1 );
  t->SetTuple1( 0, time );
  ug->GetFieldData()->AddArray( t );
}

/**
 * @brief Writes a field from @p wrapper.
 * @param[in] wrapper a wrapper around the field to be written
 * @param[in] offset the cell index offset at which to start writing data (in case of multiple subregions)
 * @param[in,out] data a VTK data container, must be a vtkAOSDataArrayTemplate of the correct value type
 */
static void
writeField( WrapperBase const & wrapper,
            localIndex const offset,
            vtkDataArray * data )
{
  types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    using T = typename ArrayType::ValueType;
    vtkAOSDataArrayTemplate< T > * typedData = vtkAOSDataArrayTemplate< T >::FastDownCast( data );
    auto const sourceArray = Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();

    forAll< parallelHostPolicy >( sourceArray.size( 0 ), [sourceArray, offset, typedData]( localIndex const i )
    {
      LvArray::forValuesInSlice( sourceArray[i], [&, compIndex = 0]( T const & value ) mutable
      {
        typedData->SetTypedComponent( offset + i, compIndex++, value );
      } );
    } );
  }, wrapper );
}

/**
 * @brief Writes a field from @p wrapper using an index list.
 * @param[in] wrapper a wrapper around the field to be written
 * @param[in] indices a list of indices into @p wrapper that will be written
 * @param[in] offset the cell index offset at which to start writing data (in case of multiple subregions)
 * @param[in,out] data a VTK data container, must be a vtkAOSDataArrayTemplate of the correct value type
 */
static void
writeField( WrapperBase const & wrapper,
            arrayView1d< localIndex const > const & indices,
            localIndex const offset,
            vtkDataArray * data )
{
  types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    using T = typename ArrayType::ValueType;
    vtkAOSDataArrayTemplate< T > * typedData = vtkAOSDataArrayTemplate< T >::FastDownCast( data );
    auto const sourceArray = Wrapper< ArrayType >::cast( wrapper ).reference().toViewConst();

    forAll< parallelHostPolicy >( indices.size(), [=]( localIndex const i )
    {
      LvArray::forValuesInSlice( sourceArray[indices[i]], [&, compIndex = 0]( T const & value ) mutable
      {
        typedData->SetTypedComponent( offset + i, compIndex++, value );
      } );
    } );
  }, wrapper );
}

/**
 * @brief Build/expand a list of strings used as default component labels on demand.
 * @param size number of labels requested
 * @return a span over range of strings (stored permanently in memory)
 */
static Span< string const >
getDefaultLabels( localIndex const size )
{
  static std::vector< string > labels;
  localIndex oldSize = LvArray::integerConversion< localIndex >( labels.size() );
  std::generate_n( std::back_inserter( labels ), size - oldSize, [&] { return std::to_string( oldSize++ ); } );
  return { labels.begin(), labels.begin() + size };
}

/**
 * @brief Checks consistency of user-provided per-dimension labels (must match the size of array).
 * @param wrapper the array wrapper
 * @param dim dimension index to check
 */
template< typename T, int NDIM, typename PERM >
void checkLabels( Wrapper< Array< T, NDIM, PERM > > const & wrapper, int const dim )
{
  GEOS_ERROR_IF_NE_MSG( LvArray::integerConversion< localIndex >( wrapper.getDimLabels( dim ).size() ),
                        wrapper.reference().size( dim ),
                        "VTK writer: component names are set, but don't match the array size.\n"
                        "This is likely a bug in physics module (solver or constitutive model)." );
}

/**
 * @brief Get a list of component labels for a particular dimension of an array.
 * @param wrapper array wrapper
 * @param dim dimension index
 * @return a span over range of strings representing labels
 */
template< typename T, int NDIM, typename PERM >
static Span< string const >
getDimLabels( Wrapper< Array< T, NDIM, PERM > > const & wrapper,
              int const dim )
{
  Span< string const > const labels = wrapper.getDimLabels( dim );
  if( labels.empty() )
  {
    return getDefaultLabels( wrapper.reference().size( dim ) );
  }
  checkLabels( wrapper, dim );
  return labels;
}

/**
 * @brief Build a multidimensional component name out of distinct dimension-wise labels.
 * @tparam Ts types of indices
 * @tparam Is dummy template argument required for positional expansion of the pack
 * @param dimLabels per-dimension component labels
 * @param indices per-dimension component indices
 * @return combined component name
 */
template< typename ... Ts, integer ... Is >
static string
makeComponentName( std::vector< string >(&dimLabels)[sizeof...( Ts )],
                   std::integer_sequence< integer, Is... >,
                   Ts const & ... indices )
{
  return stringutilities::concat( '/', dimLabels[Is][indices] ... );
}

/**
 * @brief Specialized component metadata handler for 1D arrays.
 * @param data VTK typed data array
 */
template< typename T, typename PERM >
static void
setComponentMetadata( Wrapper< Array< T, 1, PERM > > const &,
                      vtkAOSDataArrayTemplate< T > * data )
{
  data->SetNumberOfComponents( 1 );
}

/**
 * @brief Specialized component metadata handler for 2D arrays.
 * @param wrapper GEOSX typed wrapper over source array
 * @param data VTK typed data array
 *
 * This exists because we want to keep default VTK handling for unlabeled components
 * (i.e. X/Y/Z for 1-3 components, numeric indices for higher) for the time being.
 * This function can be removed if we force each physics package to always set its labels.
 */
template< typename T, typename PERM >
static void
setComponentMetadata( Wrapper< Array< T, 2, PERM > > const & wrapper,
                      vtkAOSDataArrayTemplate< T > * data )
{
  auto const view = wrapper.referenceAsView();
  data->SetNumberOfComponents( view.size( 1 ) );

  Span< string const > const labels = wrapper.getDimLabels( 1 );
  if( !labels.empty() )
  {
    checkLabels( wrapper, 1 );
    for( localIndex i = 0; i < view.size( 1 ); ++i )
    {
      data->SetComponentName( i, labels[i].c_str() );
    }
  }
}

/**
 * @brief Produces a temporary array slice from a view that can be looped over.
 * @param view the source view
 * @return a fake slice that does not point to real data but has correct dims/strides.
 * @note The slice is only valid as long as the @p view is in scope.
 *       Values in the slice may be uninitialized and should not be used.
 */
template< typename T, int NDIM, int USD >
static ArraySlice< T const, NDIM - 1, USD - 1 >
makeTemporarySlice( ArrayView< T const, NDIM, USD > const & view )
{
  // The following works in all compilers, but technically invokes undefined behavior:
  // return ArraySlice< T, NDIM - 1, USD - 1 >( nullptr, view.dims() + 1, view.strides() + 1 );
  localIndex const numComp = LvArray::indexing::multiplyAll< NDIM - 1 >( view.dims() + 1 );
  static array1d< T > arr;
  arr.template resizeWithoutInitializationOrDestruction( numComp );
  return ArraySlice< T const, NDIM - 1, USD - 1 >( arr.data(), view.dims() + 1, view.strides() + 1 );
}

/**
 * @brief Generic component metadata handler for multidimensional arrays.
 * @param wrapper GEOSX typed wrapper over source array
 * @param data VTK typed data array
 */
template< typename T, int NDIM, typename PERM >
static void
setComponentMetadata( Wrapper< Array< T, NDIM, PERM > > const & wrapper,
                      vtkAOSDataArrayTemplate< T > * data )
{
  data->SetNumberOfComponents( wrapper.numArrayComp() );

  std::vector< string > labels[NDIM-1];
  for( integer dim = 1; dim < NDIM; ++dim )
  {
    Span< string const > dimLabels = getDimLabels( wrapper, dim );
    labels[dim-1].assign( dimLabels.begin(), dimLabels.end() );
  }

  auto const view = wrapper.referenceAsView();
  auto const slice = view.size( 0 ) > 0 ? view[0] : makeTemporarySlice( view );

  integer compIndex = 0;
  LvArray::forValuesInSliceWithIndices( slice, [&]( T const &, auto const ... indices )
  {
    using idx_seq = std::make_integer_sequence< integer, sizeof...(indices) >;
    data->SetComponentName( compIndex++, makeComponentName( labels, idx_seq{}, indices ... ).c_str() );
  } );
}

template< class SUBREGION = Group >
static void
writeElementField( Group const & subRegions,
                   string const & field,
                   vtkCellData * cellData )
{
  // instantiate vtk array of the correct type
  vtkSmartPointer< vtkDataArray > data;
  localIndex numElements = 0;
  bool first = true;
  int numDims = 0;
  subRegions.forSubGroups< SUBREGION >( [&]( SUBREGION const & subRegion )
  {
    numElements += subRegion.size();
    WrapperBase const & wrapper = subRegion.getWrapperBase( field );
    if( first )
    {
      types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
      {
        using ArrayType = camp::first< decltype( tupleOfTypes ) >;
        using T = typename ArrayType::ValueType;
        auto typedData = vtkAOSDataArrayTemplate< T >::New();
        data.TakeReference( typedData );
        setComponentMetadata( Wrapper< ArrayType >::cast( wrapper ), typedData );
      }, wrapper );
      first = false;
      numDims = wrapper.numArrayDims();
    }
    else
    {
      // Sanity check
      GEOS_ERROR_IF_NE_MSG( wrapper.numArrayDims(), numDims,
                            "VTK writer: sanity check failed for " << field << " (inconsistent array dimensions)" );
      GEOS_ERROR_IF_NE_MSG( wrapper.numArrayComp(), data->GetNumberOfComponents(),
                            "VTK writer: sanity check failed for " << field << " (inconsistent array sizes)" );
    }
  } );

  data->SetNumberOfTuples( numElements );
  data->SetName( field.c_str() );

  // write each subregion in turn, keeping track of element offset
  localIndex offset = 0;
  subRegions.forSubGroups< SUBREGION >( [&]( SUBREGION const & subRegion )
  {
    WrapperBase const & wrapper = subRegion.getWrapperBase( field );
    writeField( wrapper, offset, data.GetPointer() );
    offset += subRegion.size();
  } );
  cellData->AddArray( data );
}

void VTKPolyDataWriterInterface::writeParticleFields( ParticleRegionBase const & region,
                                                      vtkCellData * cellData ) const
{
  std::unordered_set< string > materialFields;
  conduit::Node fakeRoot;
  Group materialData( "materialData", fakeRoot );
  region.forParticleSubRegions( [&]( ParticleSubRegionBase const & subRegion )
  {
    // Register a dummy group for each subregion
    Group & subReg = materialData.registerGroup( subRegion.getName() );
    subReg.resize( subRegion.size() );

    // Collect a list of plotted constitutive fields and create wrappers containing averaged data
    subRegion.getConstitutiveModels().forSubGroups( [&]( Group const & material )
    {
      material.forWrappers( [&]( WrapperBase const & wrapper )
      {
        string const fieldName = constitutive::ConstitutiveBase::makeFieldName( material.getName(), wrapper.getName() );
        if( outputUtilities::isFieldPlotEnabled( wrapper.getPlotLevel(), m_plotLevel, fieldName, m_fieldNames, m_onlyPlotSpecifiedFieldNames ) )
        {
          subReg.registerWrapper( wrapper.averageOverSecondDim( fieldName, subReg ) );
          materialFields.insert( fieldName );
        }
      } );
    } );
  } );

  // Write averaged material data
  for( string const & field : materialFields )
  {
    writeElementField( materialData, field, cellData );
  }

  // Collect a list of regular fields (filter out material field wrappers)
  // TODO: this can be removed if we stop hanging constitutive wrappers on the mesh
  std::unordered_set< string > regularFields;
  region.forParticleSubRegions( [&]( ParticleSubRegionBase const & subRegion )
  {
    for( auto const & wrapperIter : subRegion.wrappers() )
    {
      if( isFieldPlotEnabled( *wrapperIter.second ) && materialFields.count( wrapperIter.first ) == 0 )
      {
        regularFields.insert( wrapperIter.first );
      }
    }
  } );

  // Write regular fields
  for( string const & field : regularFields )
  {
    writeElementField( region.getGroup( ParticleRegionBase::viewKeyStruct::particleSubRegions() ), field, cellData );
  }
}

void VTKPolyDataWriterInterface::writeNodeFields( NodeManager const & nodeManager,
                                                  arrayView1d< localIndex const > const & nodeIndices,
                                                  vtkPointData * pointData ) const
{
  for( auto const & wrapperIter : nodeManager.wrappers() )
  {
    auto const & wrapper = *wrapperIter.second;
    if( isFieldPlotEnabled( wrapper ) )
    {
      vtkSmartPointer< vtkDataArray > data;
      types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto tupleOfTypes )
      {
        using ArrayType = camp::first< decltype( tupleOfTypes ) >;
        using T = typename ArrayType::ValueType;
        auto typedData = vtkAOSDataArrayTemplate< T >::New();
        data.TakeReference( typedData );
        setComponentMetadata( Wrapper< ArrayType >::cast( wrapper ), typedData );
      }, wrapper );

      data->SetNumberOfTuples( nodeIndices.size() );
      data->SetName( wrapper.getName().c_str() );

      writeField( wrapper, nodeIndices, 0, data.GetPointer() );
      pointData->AddArray( data );
    }
  }
}

void VTKPolyDataWriterInterface::writeElementFields( ElementRegionBase const & region,
                                                     vtkCellData * cellData ) const
{
  std::unordered_set< string > materialFields;
  conduit::Node fakeRoot;
  Group materialData( "averagedMaterialData", fakeRoot );
  region.forElementSubRegions( [&]( ElementSubRegionBase const & subRegion )
  {
    // Register a dummy group for each subregion
    Group & subReg = materialData.registerGroup( subRegion.getName() );
    subReg.resize( subRegion.size() );

    // Collect a list of plotted constitutive fields and create wrappers containing averaged data
    subRegion.getConstitutiveModels().forSubGroups( [&]( Group const & material )
    {
      material.forWrappers( [&]( WrapperBase const & wrapper )
      {
        string const fieldName = constitutive::ConstitutiveBase::makeFieldName( material.getName(), wrapper.getName() );
        if( outputUtilities::isFieldPlotEnabled( wrapper.getPlotLevel(), m_plotLevel, fieldName, m_fieldNames, m_onlyPlotSpecifiedFieldNames ) )
        {
          subReg.registerWrapper( wrapper.averageOverSecondDim( fieldName, subReg ) );
          materialFields.insert( fieldName );
        }
      } );
    } );
  } );

  // Write averaged material data
  for( string const & field : materialFields )
  {
    writeElementField( materialData, field, cellData );
  }

  // Collect a list of regular fields (filter out material field wrappers)
  // TODO: this can be removed if we stop hanging constitutive wrappers on the mesh
  std::unordered_set< string > regularFields;
  region.forElementSubRegions( [&]( ElementSubRegionBase const & subRegion )
  {
    for( auto const & wrapperIter : subRegion.wrappers() )
    {
      if( isFieldPlotEnabled( *wrapperIter.second ) && materialFields.count( wrapperIter.first ) == 0 )
      {
        regularFields.insert( wrapperIter.first );
      }
    }
  } );

  // Write regular fields
  for( string const & field : regularFields )
  {
    writeElementField( region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() ), field, cellData );
  }
}

void VTKPolyDataWriterInterface::writeCellElementRegions( real64 const time,
                                                          ElementRegionManager const & elemManager,
                                                          NodeManager const & nodeManager,
                                                          string const & path ) const
{
  elemManager.forElementRegions< CellElementRegion >( [&]( CellElementRegion const & region )
  {
    CellData VTKCells = getVtkCells( region, nodeManager.size() );
    vtkSmartPointer< vtkPoints > const VTKPoints = getVtkPoints( nodeManager, VTKCells.nodes );

    auto const ug = vtkSmartPointer< vtkUnstructuredGrid >::New();
    ug->SetCells( VTKCells.cellTypes.data(), VTKCells.cells );
    ug->SetPoints( VTKPoints );

    writeTimestamp( ug.GetPointer(), time );
    writeElementFields( region, ug->GetCellData() );
    writeNodeFields( nodeManager, VTKCells.nodes, ug->GetPointData() );

    string const regionDir = joinPath( path, region.getName() );
    writeUnstructuredGrid( regionDir, ug.GetPointer() );
  } );
}

void VTKPolyDataWriterInterface::writeParticleRegions( real64 const time,
                                                       ParticleManager const & particleManager,
                                                       string const & path ) const
{
  particleManager.forParticleRegions< ParticleRegion >( [&]( ParticleRegion const & region )
  {
    auto VTKCells = getVtkCells( region );
    auto VTKPoints = getVtkPoints( region );

    auto const ug = vtkSmartPointer< vtkUnstructuredGrid >::New();
    ug->SetPoints( VTKPoints );
    ug->SetCells( VTKCells.first.data(), VTKCells.second );

    writeTimestamp( ug.GetPointer(), time );
    writeParticleFields( region, ug->GetCellData() );

    string const regionDir = joinPath( path, region.getName() );
    writeUnstructuredGrid( regionDir, ug.GetPointer() );
  } );
}

void VTKPolyDataWriterInterface::writeWellElementRegions( real64 const time,
                                                          ElementRegionManager const & elemManager,
                                                          NodeManager const & nodeManager,
                                                          string const & path ) const
{
  elemManager.forElementRegions< WellElementRegion >( [&]( WellElementRegion const & region )
  {
    auto const & subRegion = region.getSubRegion< WellElementSubRegion >( 0 );
    ElementData well = getWell( subRegion, nodeManager );

    auto const ug = vtkSmartPointer< vtkUnstructuredGrid >::New();
    ug->SetPoints( well.points );
    ug->SetCells( well.cellTypes.data(), well.cells );

    writeTimestamp( ug.GetPointer(), time );
    writeElementFields( region, ug->GetCellData() );

    string const regionDir = joinPath( path, region.getName() );
    writeUnstructuredGrid( regionDir, ug.GetPointer() );
  } );
}

void VTKPolyDataWriterInterface::writeSurfaceElementRegions( real64 const time,
                                                             ElementRegionManager const & elemManager,
                                                             NodeManager const & nodeManager,
                                                             EmbeddedSurfaceNodeManager const & embSurfNodeManager,
                                                             FaceManager const & faceManager,
                                                             string const & path ) const
{
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & region )
  {
    auto const ug = vtkSmartPointer< vtkUnstructuredGrid >::New();
    ElementData surface = [&]()
    {
      switch( region.subRegionType() )
      {
        case SurfaceElementRegion::SurfaceSubRegionType::embeddedElement:
          {
            auto const & subRegion = region.getUniqueSubRegion< EmbeddedSurfaceSubRegion >();
            return getEmbeddedSurface( subRegion, embSurfNodeManager );
          }
        case SurfaceElementRegion::SurfaceSubRegionType::faceElement:
          {
            auto const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
            return getSurface( subRegion, nodeManager, faceManager, m_writeFaceElementsAs3D );
          }
        default:
          {
            return ElementData{};
          }
      }
    }();

    ug->SetPoints( surface.points );
    ug->SetCells( surface.cellTypes.data(), surface.cells );

    writeTimestamp( ug.GetPointer(), time );
    writeElementFields( region, ug->GetCellData() );

    string const regionDir = joinPath( path, region.getName() );
    writeUnstructuredGrid( regionDir, ug.GetPointer() );
  } );
}

static string getCycleSubFolder( integer const cycle )
{
  return GEOS_FMT( "{:06d}", cycle );
}

static string getRankFileName( integer const rank )
{
  int const width = static_cast< int >( std::log10( MpiWrapper::commSize() ) ) + 1;
  return GEOS_FMT( "rank_{0:0>{1}}", rank, width );
}

void VTKPolyDataWriterInterface::writeVtmFile( integer const cycle,
                                               DomainPartition const & domain,
                                               VTKVTMWriter const & vtmWriter ) const
{
  GEOS_ASSERT_EQ_MSG( MpiWrapper::commRank(), 0, "Must only be called on rank 0" );

  // loop over mesh bodies - use domain to get element regions
  domain.forMeshBodies( [&]( MeshBody const & meshBody )
  {
    meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
    {

      if( meshLevel.isShallowCopy() )
        return;

      string const & meshLevelName = meshLevel.getName();

      if( !m_levelNames.empty())
      {
        if( m_levelNames.find( meshLevelName ) == m_levelNames.end())
          return;
      }

      string const & meshBodyName = meshBody.getName();

      ElementRegionManager const & elemManager = meshLevel.getElemManager();

      ParticleManager const & particleManager = meshLevel.getParticleManager();

      string const meshPath = joinPath( getCycleSubFolder( cycle ), meshBodyName, meshLevelName );

      int const mpiSize = MpiWrapper::commSize();

      auto addElementRegion = [&]( ElementRegionBase const & region )
      {
        std::vector< string > const blockPath{ meshBody.getName(), meshLevel.getName(), region.getCatalogName(), region.getName() };
        string const regionPath = joinPath( meshPath, region.getName() );
        for( int i = 0; i < mpiSize; i++ )
        {
          string const dataSetName = getRankFileName( i );
          string const dataSetFile = joinPath( regionPath, dataSetName + ".vtu" );
          vtmWriter.addDataSet( blockPath, dataSetName, dataSetFile );
        }
      };

      auto addParticleRegion = [&]( ParticleRegionBase const & region )
      {
        string const & regionName = region.getName();
        std::vector< string > const blockPath{ meshBodyName, meshLevelName, region.getCatalogName(), regionName };
        string const regionPath = joinPath( meshPath, regionName );
        for( int i = 0; i < mpiSize; i++ )
        {
          string const dataSetName = getRankFileName( i );
          string const dataSetFile = joinPath( regionPath, dataSetName + ".vtu" );
          vtmWriter.addDataSet( blockPath, dataSetName, dataSetFile );
        }
      };

      // Output each of the region types
      if( m_outputRegionType == VTKRegionTypes::CELL || m_outputRegionType == VTKRegionTypes::ALL )
      {
        elemManager.forElementRegions< CellElementRegion >( addElementRegion );
      }

      if( m_outputRegionType == VTKRegionTypes::WELL || m_outputRegionType == VTKRegionTypes::ALL )
      {
        elemManager.forElementRegions< WellElementRegion >( addElementRegion );
      }

      if( m_outputRegionType == VTKRegionTypes::SURFACE || m_outputRegionType == VTKRegionTypes::ALL )
      {
        elemManager.forElementRegions< SurfaceElementRegion >( addElementRegion );
      }

      if( m_outputRegionType == VTKRegionTypes::PARTICLE || m_outputRegionType == VTKRegionTypes::ALL )
      {
        particleManager.forParticleRegions< ParticleRegion >( addParticleRegion );
      }
    } );
  } );

  vtmWriter.write();
}

int toVtkOutputMode( VTKOutputMode const mode )
{
  switch( mode )
  {
    case VTKOutputMode::ASCII: return vtkXMLWriterBase::Ascii;
    case VTKOutputMode::BINARY: return vtkXMLWriterBase::Binary;
    default:
    {
      GEOS_ERROR( "Unsupported VTK output mode" );
      return -1;
    }
  }
}

void VTKPolyDataWriterInterface::writeUnstructuredGrid( string const & path,
                                                        vtkUnstructuredGrid * ug ) const
{
  vtkSmartPointer< vtkAlgorithm > filter;

  // If we want to get rid of the ghost ranks, we use the appropriate `vtkThreshold` filter.
  // If we don't, to keep the symetry in the code, we use a `vtkPassThrough`
  // that will allow a more generic code down the line.
  if( !m_writeGhostCells && ug->GetCellData()->HasArray( ObjectManagerBase::viewKeyStruct::ghostRankString() ) )
  {
    auto threshold = vtkSmartPointer< vtkThreshold >::New();
    // Ghost ranks values are integers, and negative values mean that the cell is owned by another rank.
    // Removing the cells with negative ghost ranks remove duplicated cells in the vtk output.
    threshold->SetUpperThreshold( -0.5 );
    threshold->SetInputArrayToProcess( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, ObjectManagerBase::viewKeyStruct::ghostRankString() );

    filter = threshold;
  }
  else
  {
    filter = vtkSmartPointer< vtkPassThrough >::New();
  }

  filter->SetInputDataObject( ug );
  filter->Update();

  makeDirectory( path );
  string const vtuFilePath = joinPath( path, getRankFileName( MpiWrapper::commRank() ) + ".vtu" );
  auto const vtuWriter = vtkSmartPointer< vtkXMLUnstructuredGridWriter >::New();
  vtuWriter->SetInputData( filter->GetOutputDataObject( 0 ) );
  vtuWriter->SetFileName( vtuFilePath.c_str() );
  vtuWriter->SetDataMode( toVtkOutputMode( m_outputMode ) );
  vtuWriter->Write();
}

void VTKPolyDataWriterInterface::write( real64 const time,
                                        integer const cycle,
                                        DomainPartition const & domain )
{
  // This guard prevents crashes due to a floating point exception (SIGFPE)
  // triggered inside VTK by a progress indicator
  LvArray::system::FloatingPointExceptionGuard guard;

  string const stepSubDir = joinPath( m_outputName, getCycleSubFolder( cycle ) );
  string const stepSubDirFull = joinPath( m_outputDir, stepSubDir );

  int const rank = MpiWrapper::commRank();
  if( rank == 0 )
  {
    makeDirsForPath( stepSubDirFull );
  }
  MpiWrapper::barrier( MPI_COMM_GEOS );

  // loop over all mesh levels and mesh bodies
  domain.forMeshBodies( [&]( MeshBody const & meshBody )
  {
    meshBody.forMeshLevels( [&]( MeshLevel const & meshLevel )
    {

      if( meshLevel.isShallowCopy() )
        return;

      string const & meshLevelName = meshLevel.getName();

      if( !m_levelNames.empty())
      {
        if( m_levelNames.find( meshLevelName ) == m_levelNames.end())
          return;
      }

      ElementRegionManager const & elemManager = meshLevel.getElemManager();
      ParticleManager const & particleManager = meshLevel.getParticleManager();
      NodeManager const & nodeManager = meshLevel.getNodeManager();
      FaceManager const & faceManager = meshLevel.getFaceManager();
      EmbeddedSurfaceNodeManager const & embSurfNodeManager = meshLevel.getEmbSurfNodeManager();
      string const & meshBodyName = meshBody.getName();

      if( m_requireFieldRegistrationCheck && !m_fieldNames.empty() )
      {
        outputUtilities::checkFieldRegistration( elemManager,
                                                 nodeManager,
                                                 m_fieldNames,
                                                 "VTKOutput" );
        m_requireFieldRegistrationCheck = false;
      }

      string const meshDir = joinPath( stepSubDirFull, meshBodyName, meshLevelName );
      makeDirsForPath( meshDir );

      if( m_outputRegionType == VTKRegionTypes::CELL || m_outputRegionType == VTKRegionTypes::ALL )
      {
        writeCellElementRegions( time, elemManager, nodeManager, meshDir );
      }
      if( m_outputRegionType == VTKRegionTypes::WELL || m_outputRegionType == VTKRegionTypes::ALL )
      {
        writeWellElementRegions( time, elemManager, nodeManager, meshDir );
      }
      if( m_outputRegionType == VTKRegionTypes::SURFACE || m_outputRegionType == VTKRegionTypes::ALL )
      {
        writeSurfaceElementRegions( time, elemManager, nodeManager, embSurfNodeManager, faceManager, meshDir );
      }
      if( m_outputRegionType == VTKRegionTypes::PARTICLE || m_outputRegionType == VTKRegionTypes::ALL )
      {
        writeParticleRegions( time, particleManager, meshDir );
      }
    } );
  } );

  if( rank == 0 )
  {
    string const vtmName = stepSubDir + ".vtm";
    VTKVTMWriter vtmWriter( joinPath( m_outputDir, vtmName ) );
    writeVtmFile( cycle, domain, vtmWriter );

    if( cycle != m_previousCycle )
    {
      m_pvd.addData( time, vtmName );
      m_pvd.save();
    }
  }

  m_previousCycle = cycle;
}

void VTKPolyDataWriterInterface::clearData()
{
  m_pvd.reinitData();
}

bool VTKPolyDataWriterInterface::isFieldPlotEnabled( dataRepository::WrapperBase const & wrapper ) const
{
  return outputUtilities::isFieldPlotEnabled( wrapper.getPlotLevel(),
                                              m_plotLevel,
                                              wrapper.getName(),
                                              m_fieldNames,
                                              m_onlyPlotSpecifiedFieldNames );
}

} // namespace vtk
} // namespace geos
