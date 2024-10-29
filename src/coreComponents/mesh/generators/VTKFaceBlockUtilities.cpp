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

#include "VTKFaceBlockUtilities.hpp"

#include "mesh/generators/VTKUtilities.hpp"
#include "mesh/generators/CollocatedNodes.hpp"

#include "dataRepository/Group.hpp"

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkExtractEdges.h>
#include <vtkGeometryFilter.h>
#include <vtkPointData.h>

#include <algorithm>

namespace geos::vtk
{

namespace internal
{

/**
 * @brief This classes builds the mapping from global (vtk) element indices to the indices local to the cell blocks.
 * Then the global element to global faces can be be also computed.
 */
class ElementToFace
{
public:
  /**
   * @brief Builds from the cell blocks.
   * @param[in] cellBlocks The cell blocks group.
   */
  explicit ElementToFace( dataRepository::Group const & cellBlocks )
  {
    localIndex const numCellBlocks = cellBlocks.numSubGroups();
    for( int c = 0; c < numCellBlocks; ++c )
    {
      CellBlock const & cb = cellBlocks.getGroup< CellBlock >( c );
      auto const & l2g = cb.localToGlobalMapConstView();

      std::map< globalIndex, localIndex > g2l;

      for( auto l = 0; l < l2g.size(); ++l )
      {
        globalIndex const & g = l2g[l];
        m_elementToCellBlock[g] = c;
        g2l[g] = l;
      }

      m_cbe[c] = g2l;
      m_cbf[c] = cb.getElemToFacesConstView();
    }
  }

  /**
   * @brief Given the global (vtk) id @p ei, returns the cell block index it belongs to.
   * @param[in] ei global cell index.
   * @return The cell block index.
   */
  int getCellBlockIndex( vtkIdType const & ei ) const
  {
    return m_elementToCellBlock.at( ei );
  }

  /**
   * @brief From global index @p ei, returns the local element index in its cell block.
   * @param ei The global vtk element index.
   * @return The integer.
   */
  localIndex getElementIndexInCellBlock( vtkIdType const & ei ) const
  {
    localIndex const & cbi = getCellBlockIndex( ei );
    return m_cbe.at( cbi ).at( ei );
  }

  /**
   * @brief Returns the faces indices touching the global (vtk) element id.
   * @param[in] ei global cell index.
   * @return A container for the face indices of element @p ei.
   */
  auto operator[]( vtkIdType const & ei ) const
  {
    localIndex const & cbi = getCellBlockIndex( ei );
    localIndex const & e = m_cbe.at( cbi ).at( ei );
    arrayView2d< localIndex const > const & e2f = m_cbf.at( cbi );
    return e2f[e];
  }

private:
  /// global element index to the local cell block index
  std::map< globalIndex, localIndex > m_elementToCellBlock;

  /// Cell block index to a mapping from global element index to the local (to the cell block) element index.
  std::map< localIndex, std::map< globalIndex, localIndex > > m_cbe;

  /// Cell block index to a mapping from global element index to the faces indices.
  std::map< localIndex, arrayView2d< localIndex const > > m_cbf;
};

} // end of namespace internal


/**
 * @brief Organize the collocated nodes information as an @p LvArray::ArrayOfArrays.
 * @param cns The collocated nodes information.
 * @return An iterable of arrays. Each array containing the global indices of the nodes which are collocated of each others.
 */
ArrayOfArrays< globalIndex > buildCollocatedNodesMap( CollocatedNodes const & cns )
{
  ArrayOfArrays< globalIndex > result;

  std::vector< int > sizes( cns.size() );
  for( std::size_t i = 0; i < cns.size(); ++i )
  {
    sizes[i] = cns[i].size();
  }
  result.resizeFromCapacities< serialPolicy >( sizes.size(), sizes.data() );

  for( std::size_t i = 0; i < cns.size(); ++i )
  {
    for( globalIndex const & g: cns[i] )
    {
      result.emplaceBack( i, g );
    }
  }

  return result;
}


/**
 * @brief Build the mapping between the elements and their nodes.
 * @param mesh The mesh for which to compute the mapping.
 * @return The mapping (that allows heterogeneous meshes).
 */
ArrayOfArrays< localIndex > build2dElemTo2dNodes( vtkSmartPointer< vtkDataSet > mesh )
{
  // First allocate...
  vtkIdType const numCells = mesh->GetNumberOfCells();
  std::vector< localIndex > sizes( numCells );
  for( auto i = 0; i < numCells; ++i )
  {
    sizes[i] = mesh->GetCell( i )->GetNumberOfPoints();
  }
  ArrayOfArrays< localIndex > result;
  result.resizeFromCapacities< geos::serialPolicy >( sizes.size(), sizes.data() );
  // ... then fill with data.
  for( auto i = 0; i < numCells; ++i )
  {
    vtkIdList * const pointIds = mesh->GetCell( i )->GetPointIds();
    vtkIdType const numPoints = pointIds->GetNumberOfIds();
    for( int j = 0; j < numPoints; ++j )
    {
      result.emplaceBack( i, pointIds->GetId( j ) );
    }
  }

  return result;
}


ArrayOfArrays< array1d< globalIndex > > buildCollocatedNodesBucketsOf2dElemsMap( ArrayOfArrays< localIndex > const & elem2dTo2dNodes,
                                                                                 ArrayOfArrays< globalIndex > const & nodes2dToCollocatedNodes )
{
  localIndex const num2dElems = elem2dTo2dNodes.size();
  ArrayOfArrays< array1d< globalIndex > > result;
  // Allocation...
  result.resize( num2dElems );
  for( localIndex e2d = 0; e2d < num2dElems; ++e2d )
  {
    result.resizeArray( e2d, elem2dTo2dNodes[e2d].size() );
  }
  for( localIndex e2d = 0; e2d < num2dElems; ++e2d )
  {
    auto const numNodes = elem2dTo2dNodes[e2d].size();
    for( integer ni = 0; ni < numNodes; ++ni )
    {
      auto & dest = result( e2d, ni );
      localIndex const node = elem2dTo2dNodes[e2d][ni];
      auto const src = nodes2dToCollocatedNodes[node];
      dest.reserve( src.size() );
      // ...Definition
      for( globalIndex const s: src )
      {
        dest.emplace_back( s );
      }
    }
  }

  return result;
}


/**
 * @brief Some hash function for pairs of integers.
 */
struct pairHashComputer
{
  /**
   * @brief Computes the hash of a pair of integers.
   * @param p The input pair.
   * @return The hash.
   */
  size_t operator()( const std::pair< vtkIdType, vtkIdType > & p ) const
  {
    auto hash1 = std::hash< vtkIdType >{} ( p.first );
    auto hash2 = std::hash< vtkIdType >{} ( p.second );

    // If hash1 == hash2, their XOR is zero.
    return hash1 != hash2 ? hash1 ^ hash2: hash1;
  }
};


/**
 * @brief Compute the 2d faces (geometrical edges in 2d) to 2d elems (geometrical faces in 3d) mapping.
 * @param edges The edges of the fracture mesh.
 * @param faceMesh The fracture mesh.
 * @return The mapping
 */
ArrayOfArrays< localIndex > buildFace2dToElems2d( vtkPolyData * edges,
                                                  vtkSmartPointer< vtkDataSet > faceMesh )
{
  // Each edge is first associated to an hash and to an id.
  // Then we loop over all the edges of each cell and compute its hash to recover the associated id.
  std::unordered_map< std::pair< vtkIdType, vtkIdType >, int, pairHashComputer > face2dIds;
  for( int i = 0; i < edges->GetNumberOfCells(); ++i )
  {
    vtkCell * c = edges->GetCell( i );
    vtkIdList * edgePointIds = c->GetPointIds();
    GEOS_ASSERT( edgePointIds->GetNumberOfIds() == 2 );
    std::pair< vtkIdType, vtkIdType > const minMax = std::minmax( edgePointIds->GetId( 0 ), edgePointIds->GetId( 1 ) );
    face2dIds[minMax] = i;
  }

  ArrayOfArrays< localIndex > face2dToElems2d( LvArray::integerConversion< localIndex >( edges->GetNumberOfCells() ) );
  for( vtkIdType i = 0; i < faceMesh->GetNumberOfCells(); ++i )
  {
    vtkCell * c = faceMesh->GetCell( i );
    for( int j = 0; j < c->GetNumberOfEdges(); ++j )
    {
      vtkCell * e = c->GetEdge( j );
      vtkIdList * edgePointIds = e->GetPointIds();
      GEOS_ASSERT( edgePointIds->GetNumberOfIds() == 2 );
      std::pair< vtkIdType, vtkIdType > const minMax = std::minmax( edgePointIds->GetId( 0 ), edgePointIds->GetId( 1 ) );
      face2dToElems2d.emplaceBack( LvArray::integerConversion< localIndex >( face2dIds.at( minMax ) ), i );
    }
  }

  return face2dToElems2d;
}


/**
 * @brief For each 2d face (a segment in 3d), returns one single overlapping 3d edge (and not 2).
 * @param globalPtIds[in] The global point ids.
 * @param edges[in] The edges as computed by vtk.
 * @param collocatedNodes[in] The collocated nodes information.
 * @param nodeToEdges[in] The node to edges mapping.
 * @return The 2d face to 3d edge mapping. In the case where the face is at the boundary of the MPI domain,
 * then the edge index will be set to @e -1 for further actions.
 */
array1d< localIndex > buildFace2dToEdge( vtkIdTypeArray const * globalPtIds,
                                         vtkPolyData * edges,
                                         CollocatedNodes const & collocatedNodes,
                                         ArrayOfArraysView< localIndex const > nodeToEdges )
{
  std::map< globalIndex, std::vector< localIndex > > n2e;
  for( auto i = 0; i < nodeToEdges.size(); ++i )
  {
    std::vector< localIndex > es;
    for( auto j = 0; j < nodeToEdges[i].size(); ++j )
    {
      es.push_back( nodeToEdges[i][j] );
    }
    n2e[globalPtIds->GetValue( i )] = es;
  }

  auto const comp = []( std::pair< vtkIdType, int > const & l, std::pair< vtkIdType, int > const & r ) { return l.second < r.second; };
  array1d< localIndex > face2dToEdge( edges->GetNumberOfCells() );
  // We loop over all the (duplicated) nodes of each edge.
  // Then, thanks to the node to edges mapping,
  // we can find all (ie 2) the edges that share (at max) 2 duplicated nodes.
  // Eventually we select one of these edges.
  for( int i = 0; i < edges->GetNumberOfCells(); ++i )
  {
    std::vector< vtkIdType > allDuplicatedNodesOfEdge;
    vtkCell * edge = edges->GetCell( i );
    for( int j = 0; j < edge->GetNumberOfPoints(); ++j )
    {
      for( auto const & d: collocatedNodes[ edge->GetPointId( j ) ] )
      {
        allDuplicatedNodesOfEdge.emplace_back( d );
      }
    }
    std::map< vtkIdType, int > edgeCount;
    for( vtkIdType const & d: allDuplicatedNodesOfEdge )
    {
      localIndex const dd = LvArray::integerConversion< localIndex >( d );
      for( localIndex const & val: n2e[dd] )
      {
        edgeCount[val]++;
      }
    }
    auto const res = std::max_element( edgeCount.cbegin(), edgeCount.cend(), comp );
    // If we're in a case where there aren't two edges sharing two nodes,
    // then it means that we're in a corner case where the 2d element is on the boundary of the MPI domain,
    // and maybe some nodes are missing for the 2d element to be properly and consistently defines.
    // In this case, we explicitly set the edge index at `-1`, so we can get back on it later.
    face2dToEdge[i] = res->second < 2 ? -1: LvArray::integerConversion< localIndex >( res->first );
  }

  return face2dToEdge;
}

/**
 * @brief Compute the 2d elements to 2d faces mapping (mainly by inverting @p face2dToElems2d).
 * @param num2dElements Number of (2d) elements in the fracture.
 * @param face2dToElems2d Which 2d faces touch which 2d elements.
 * @return The inverted mapping.
 */
ArrayOfArrays< localIndex > buildElem2dToFace2d( vtkIdType num2dElements,
                                                 ArrayOfArraysView< localIndex const > face2dToElems2d )
{
  // Inversion of the input mapping.
  ArrayOfArrays< localIndex > elem2dToFace2d( LvArray::integerConversion< localIndex >( num2dElements ) );
  for( localIndex face2dIndex = 0; face2dIndex < face2dToElems2d.size(); ++face2dIndex )
  {
    for( localIndex const & elem2dIndex: face2dToElems2d[face2dIndex] )
    {
      elem2dToFace2d.emplaceBack( elem2dIndex, face2dIndex );
    }
  }

  return elem2dToFace2d;
}

/**
 * @brief Builds the 2d element to edges mapping.
 * @param num2dElements Number of (2d) elements in the fracture.
 * @param face2dToEdge The 2d face to 3d edge mapping.
 * @param elem2dToFace2d The 2d element to 2d face mapping.
 * @return The computed mapping.
 *
 * @see FaceBlockABC::get2dElemToEdges for more information.
 */
ArrayOfArrays< localIndex > buildElem2dToEdges( vtkIdType num2dElements,
                                                arrayView1d< localIndex const > face2dToEdge,
                                                ArrayOfArraysView< localIndex const > elem2dToFace2d )
{
  ArrayOfArrays< localIndex > elem2dToEdges( LvArray::integerConversion< localIndex >( num2dElements ) );
  // Let's compose elem2dToFace2d with face2dToEdge to create elem2dToEdges.
  for( localIndex elemIndex = 0; elemIndex < elem2dToFace2d.size(); ++elemIndex )
  {
    for( auto const & face2dIndex: elem2dToFace2d[elemIndex] )
    {
      localIndex const & e = face2dToEdge[face2dIndex];
      if( e > -1 )
      {
        elem2dToEdges.emplaceBack( elemIndex, e );
      }
    }
  }

  return elem2dToEdges;
}

/**
 * @brief Utility structure which connects 2d element to their 3d element, faces and nodes neighbors.
 */
struct Elem2dTo3dInfo
{
  ToCellRelation< ArrayOfArrays< localIndex > > elem2dToElem3d;
  ArrayOfArrays< localIndex > elem2dToFaces;
  /**
   * @brief All the neighboring points of the 2d element @e available on the current MPI rank.
   * @note During the MPI partitioning, a 2d element may not already have its neighbors (face, 3d element).
   * But this does not mean that the 2d element has no neighboring points.
   * @note This container will not keep the order of the nodes.
   * @note This container may not contain all the nodes of the 2d element.
   * The missing nodes must therefore be added as part of the ghosting process.
   */
  ArrayOfArrays< localIndex > elem2dToNodes;

  Elem2dTo3dInfo( ToCellRelation< ArrayOfArrays< localIndex > > && elem2dToElem3d_,
                  ArrayOfArrays< localIndex > && elem2dToFaces_,
                  ArrayOfArrays< localIndex > && elem2dToNodes_ )
    : elem2dToElem3d( elem2dToElem3d_ ),
    elem2dToFaces( elem2dToFaces_ ),
    elem2dToNodes( elem2dToNodes_ )
  { }
};

/**
 * @brief Computes the mappings that link the 2d elements to their 3d element and faces neighbors.
 * @param faceMesh The face mesh.
 * @param mesh The 3d mesh.
 * @param collocatedNodes The collocated nodes information.
 * @param faceToNodes The (3d) face to nodes mapping.
 * @param elemToFaces The element to faces information.
 * @return All the information gathered into a single instance.
 */
Elem2dTo3dInfo buildElem2dTo3dElemAndFaces( vtkSmartPointer< vtkDataSet > faceMesh,
                                            vtkSmartPointer< vtkDataSet > mesh,
                                            CollocatedNodes const & collocatedNodes,
                                            ArrayOfArraysView< localIndex const > faceToNodes,
                                            vtk::internal::ElementToFace const & elemToFaces )
{
  // First, we'll only consider the boundary cells,
  // since only boundary cells can be involved in this kind of computations.
  const char * const boundaryPointsName = "boundary points";
  const char * const boundaryCellsName = "boundary cells";

  auto boundaryExtractor = vtkSmartPointer< vtkGeometryFilter >::New();
  boundaryExtractor->PassThroughPointIdsOn();
  boundaryExtractor->PassThroughCellIdsOn();
  boundaryExtractor->FastModeOff();
  boundaryExtractor->SetOriginalPointIdsName( boundaryPointsName );
  boundaryExtractor->SetOriginalCellIdsName( boundaryCellsName );
  boundaryExtractor->SetInputData( mesh );
  boundaryExtractor->Update();
  vtkPolyData * boundary = boundaryExtractor->GetOutput();

  vtkIdTypeArray const * boundaryPoints = vtkIdTypeArray::FastDownCast( boundary->GetPointData()->GetArray( boundaryPointsName ) );
  vtkIdTypeArray const * boundaryCells = vtkIdTypeArray::FastDownCast( boundary->GetCellData()->GetArray( boundaryCellsName ) );

  vtkIdTypeArray const * globalPtIds = vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );
  vtkIdTypeArray const * globalCellIds = vtkIdTypeArray::FastDownCast( mesh->GetCellData()->GetGlobalIds() );

  std::map< vtkIdType, localIndex > ng2l;  // global to local mapping for nodes.
  for( vtkIdType i = 0; i < globalPtIds->GetNumberOfValues(); ++i )
  {
    ng2l[globalPtIds->GetValue( i )] = i;
  }

  // Let's build the elem2d to elem3d mapping.
  // We need to find the 3d elements (and only the 3d elements, so we can safely ignore the others).
  // First we compute the mapping from all the boundary nodes to the 3d elements that rely on those nodes.
  std::map< vtkIdType, std::vector< vtkIdType > > nodesToCellsFull;
  for( vtkIdType i = 0; i < boundary->GetNumberOfCells(); ++i )
  {
    vtkIdType const cellId = boundaryCells->GetValue( i );
    if( mesh->GetCell( cellId )->GetCellDimension() != 3 )
    {
      continue;
    }
    vtkIdList * pointIds = boundary->GetCell( i )->GetPointIds();
    for( int j = 0; j < pointIds->GetNumberOfIds(); ++j )
    {
      vtkIdType const pointId = boundaryPoints->GetValue( pointIds->GetId( j ) );
      nodesToCellsFull[globalPtIds->GetValue( pointId )].emplace_back( globalCellIds->GetValue( cellId ) );
    }
  }

  // Then we only keep the duplicated nodes. It's only for optimisation purpose.
  std::map< vtkIdType, std::set< vtkIdType > > nodesToCells;
  { // scope reduction
    std::set< vtkIdType > allDuplicatedNodes;
    for( std::size_t i = 0; i < collocatedNodes.size(); ++i )
    {
      std::vector< vtkIdType > const & ns = collocatedNodes[ i ];
      allDuplicatedNodes.insert( ns.cbegin(), ns.cend() );
    }

    for( vtkIdType const & n: allDuplicatedNodes )
    {
      auto const it = nodesToCellsFull.find( n );
      if( it != nodesToCellsFull.cend() )
      {
        std::vector< vtkIdType > const & tmp = it->second;
        std::set< vtkIdType > const cells{ tmp.cbegin(), tmp.cend() };
        nodesToCells[n] = cells;
      }
    }
  }

  vtkIdType const num2dElements = faceMesh->GetNumberOfCells();

  ArrayOfArrays< localIndex > elem2dToElem3d( num2dElements, 2 );
  ArrayOfArrays< localIndex > elem2dToCellBlock( num2dElements, 2 );
  ArrayOfArrays< localIndex > elem2dToFaces( num2dElements, 2 );
  ArrayOfArrays< localIndex > elem2dToNodes( num2dElements, 10 );

  // Now we loop on all the 2d elements.
  for( int e2d = 0; e2d < num2dElements; ++e2d )
  {
    // We collect all the duplicated points that are involved for each 2d element.
    vtkIdList * pointIds = faceMesh->GetCell( e2d )->GetPointIds();
    std::size_t const elem2dNumPoints = pointIds->GetNumberOfIds();
    // All the duplicated points of the 2d element. Note that we lose the collocation of the duplicated nodes.
    std::set< vtkIdType > duplicatedPointOfElem2d;
    for( vtkIdType j = 0; j < pointIds->GetNumberOfIds(); ++j )
    {
      std::vector< vtkIdType > const & ns = collocatedNodes[ pointIds->GetId( j ) ];
      duplicatedPointOfElem2d.insert( ns.cbegin(), ns.cend() );
    }

    for( vtkIdType const & gni: duplicatedPointOfElem2d )
    {
      auto it = ng2l.find( gni );
      if( it != ng2l.cend() )  // If the node is not on this rank, we want to ignore this entry.
      {
        // The node lists in `elem2dToNodes` may be in any order.
        // Anyway, there will be a specific step to reset the appropriate order.
        // Note that due to ghosting, there will also always be some cases where
        // some nodes will be missing in `elem2dToNodes` _before_ the ghosting.
        // After the ghosting, the whole information will be gathered,
        // but a reordering step will always be compulsory.
        elem2dToNodes.emplaceBack( e2d, it->second );
      }
    }

    // Here, we collect all the 3d elements that are concerned by at least one of those duplicated elements.
    std::map< vtkIdType, std::set< vtkIdType > > elem3dToDuplicatedNodes;
    for( vtkIdType const & n: duplicatedPointOfElem2d )
    {
      auto const ncs = nodesToCells.find( n );
      if( ncs != nodesToCells.cend() )
      {
        for( vtkIdType const & c: ncs->second )
        {
          elem3dToDuplicatedNodes[c].insert( n );
        }
      }
    }
    // Last we extract which of those candidate 3d elements are the ones actually neighboring the 2d element.
    for( auto const & e2n: elem3dToDuplicatedNodes )
    {
      // If the face of the element 3d has the same number of nodes than the elem 2d, it should be a successful (the mesh is conformal).
      if( e2n.second.size() == elem2dNumPoints )
      {
        // Now we know that the element 3d has a face that touches the element 2d. Let's find which one.
        elem2dToElem3d.emplaceBack( e2d, elemToFaces.getElementIndexInCellBlock( e2n.first ) );
        // Computing the elem2dToFaces mapping.
        auto faces = elemToFaces[e2n.first];
        for( int j = 0; j < faces.size( 0 ); ++j )
        {
          localIndex const faceIndex = faces[j];
          auto nodes = faceToNodes[faceIndex];
          std::set< vtkIdType > globalNodes;
          for( auto const & n: nodes )
          {
            globalNodes.insert( globalPtIds->GetValue( n ) );
          }
          if( globalNodes == e2n.second )
          {
            elem2dToFaces.emplaceBack( e2d, faceIndex );
            elem2dToCellBlock.emplaceBack( e2d, elemToFaces.getCellBlockIndex( e2n.first ) );
            break;
          }
        }
      }
    }
  }

  auto cellRelation = ToCellRelation< ArrayOfArrays< localIndex > >( std::move( elem2dToCellBlock ), std::move( elem2dToElem3d ) );
  return Elem2dTo3dInfo( std::move( cellRelation ), std::move( elem2dToFaces ), std::move( elem2dToNodes ) );
}


/**
 * @brief Computes the local to global mappings for the 2d elements of the face mesh.
 * @param faceMeshCellGlobalIds The cell global ids for the face mesh.
 * @param meshCellGlobalIds The cell global ids for the volumic mesh.
 * @return The mapping as an array.
 * @details The vtk global ids of the elements of the @p faceMesh are used,
 * but they are shifted by the maximum element global id of the @p mesh,
 * to avoid any collision.
 */
array1d< globalIndex > buildLocalToGlobal( vtkIdTypeArray const * faceMeshCellGlobalIds,
                                           vtkIdTypeArray const * meshCellGlobalIds )
{
  vtkIdType const numCells = faceMeshCellGlobalIds ? faceMeshCellGlobalIds->GetNumberOfTuples() : 0;
  array1d< globalIndex > l2g( numCells );

  // In order to avoid any cell global id collision, we gather the max cell global id over all the ranks.
  // Then we use this maximum as on offset.
  // TODO This does not take into account multiple fractures.
  vtkIdType const maxLocalCellId = meshCellGlobalIds->GetMaxId();
  vtkIdType const maxGlobalCellId = MpiWrapper::max( maxLocalCellId );
  vtkIdType const cellGlobalOffset = maxGlobalCellId + 1;

  for( auto i = 0; i < l2g.size(); ++i )
  {
    // Note that `l2g.size()` is zero if `faceMeshGlobalIds` is 0 too.
    // This prevents from calling a member on a null instance.
    l2g[i] = faceMeshCellGlobalIds->GetValue( i ) + cellGlobalOffset;
  }

  return l2g;
}


void importFractureNetwork( string const & faceBlockName,
                            vtkSmartPointer< vtkDataSet > faceMesh,
                            vtkSmartPointer< vtkDataSet > mesh,
                            CellBlockManager & cellBlockManager )
{
  ArrayOfArrays< localIndex > const faceToNodes = cellBlockManager.getFaceToNodes();
  vtk::internal::ElementToFace const elemToFaces( cellBlockManager.getCellBlocks() );
  ArrayOfArrays< localIndex > const nodeToEdges = cellBlockManager.getNodeToEdges();

  CollocatedNodes const collocatedNodes( faceBlockName, faceMesh );
  // Add the appropriate validations (only 2d cells...)

  auto edgesExtractor = vtkSmartPointer< vtkExtractEdges >::New();
  edgesExtractor->SetInputData( faceMesh );
  edgesExtractor->UseAllPointsOn();  // Important: we want to prevent any node renumbering.
  edgesExtractor->Update();
  vtkPolyData * edges = edgesExtractor->GetOutput();

  vtkIdType const num2dFaces = edges->GetNumberOfCells();
  vtkIdType const num2dElements = faceMesh->GetNumberOfCells();
  // Now let's build the elem2dTo* mappings.
  Elem2dTo3dInfo elem2dTo3d = buildElem2dTo3dElemAndFaces( faceMesh, mesh, collocatedNodes, faceToNodes.toViewConst(), elemToFaces );

  ArrayOfArrays< localIndex > face2dToElems2d = buildFace2dToElems2d( edges, faceMesh );
  array1d< localIndex > face2dToEdge = buildFace2dToEdge( vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() ), edges, collocatedNodes, nodeToEdges.toViewConst() );
  ArrayOfArrays< localIndex > const elem2dToFace2d = buildElem2dToFace2d( num2dElements, face2dToElems2d.toViewConst() );
  ArrayOfArrays< localIndex > elem2dToEdges = buildElem2dToEdges( num2dElements, face2dToEdge.toViewConst(), elem2dToFace2d.toViewConst() );

  // Mappings are now computed. Just create the face block by value.
  FaceBlock & faceBlock = cellBlockManager.registerFaceBlock( faceBlockName );

  faceBlock.setNum2dElements( num2dElements );
  faceBlock.setNum2dFaces( num2dFaces );
  faceBlock.set2dElemToNodes( std::move( elem2dTo3d.elem2dToNodes ) );
  faceBlock.set2dElemToEdges( std::move( elem2dToEdges ) );
  faceBlock.set2dFaceToEdge( std::move( face2dToEdge ) );
  faceBlock.set2dFaceTo2dElems( std::move( face2dToElems2d ) );

  faceBlock.set2dElemToFaces( std::move( elem2dTo3d.elem2dToFaces ) );
  faceBlock.set2dElemToElems( std::move( elem2dTo3d.elem2dToElem3d ) );

  faceBlock.setLocalToGlobalMap(
    buildLocalToGlobal( vtkIdTypeArray::FastDownCast( faceMesh->GetCellData()->GetGlobalIds() ),
                        vtkIdTypeArray::FastDownCast( mesh->GetCellData()->GetGlobalIds() ) )
    );

  faceBlock.set2dElemsToCollocatedNodesBuckets(
    buildCollocatedNodesBucketsOf2dElemsMap( build2dElemTo2dNodes( faceMesh ),
                                             buildCollocatedNodesMap( collocatedNodes ) )
    );
}

} // end of namespace
