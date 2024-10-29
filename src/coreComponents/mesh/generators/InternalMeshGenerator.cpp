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

#include "InternalMeshGenerator.hpp"
#include "CellBlockManager.hpp"

#include "common/MpiWrapper.hpp"

#include "common/DataTypes.hpp"

#include "LvArray/src/tensorOps.hpp"

#include <cmath>

namespace geos
{
using namespace dataRepository;

InternalMeshGenerator::InternalMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_dim( 3 ),
  m_min(),
  m_max(),
  m_coordinatePrecision( 1e-10 ),
  m_vertices{},
  m_nElems{},
  m_nElemBias{},
  m_setCoords{},
  m_periodic( 3 )
{
  registerWrapper( viewKeyStruct::xCoordsString(), &(m_vertices[0]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "x-coordinates of each mesh block vertex" );

  registerWrapper( viewKeyStruct::yCoordsString(), &(m_vertices[1]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "y-coordinates of each mesh block vertex" );

  registerWrapper( viewKeyStruct::zCoordsString(), &(m_vertices[2]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "z-coordinates of each mesh block vertex" );

  registerWrapper( viewKeyStruct::xElemsString(), &(m_nElems[0]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Number of elements in the x-direction within each mesh block" );

  registerWrapper( viewKeyStruct::yElemsString(), &(m_nElems[1]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Number of elements in the y-direction within each mesh block" );

  registerWrapper( viewKeyStruct::zElemsString(), &(m_nElems[2]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Number of elements in the z-direction within each mesh block" );

  registerWrapper( viewKeyStruct::xBiasString(), &(m_nElemBias[0]) ).
    setApplyDefaultValue( 1.0 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the x-direction within each mesh block (dx_left=(1+b)*L/N, dx_right=(1-b)*L/N)" );

  registerWrapper( viewKeyStruct::yBiasString(), &(m_nElemBias[1]) ).
    setApplyDefaultValue( 1.0 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the y-direction within each mesh block (dy_left=(1+b)*L/N, dx_right=(1-b)*L/N)" );

  registerWrapper( viewKeyStruct::zBiasString(), &(m_nElemBias[2]) ).
    setApplyDefaultValue( 1.0 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the z-direction within each mesh block (dz_left=(1+b)*L/N, dz_right=(1-b)*L/N)" );

  registerWrapper( viewKeyStruct::cellBlockNamesString(), &m_regionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of each mesh block" );

  registerWrapper( viewKeyStruct::elementTypesString(), &m_elementType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Element types of each mesh block" );

  registerWrapper( viewKeyStruct::trianglePatternString(), &m_trianglePattern ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Pattern by which to decompose the hex mesh into wedges" );

  registerWrapper( viewKeyStruct::positionToleranceString(), &m_coordinatePrecision ).
    setApplyDefaultValue( 1e-10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "A position tolerance to verify if a node belong to a nodeset" );

  registerWrapper( viewKeyStruct::periodicString(), &m_periodic ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Flag for periodicity in each direction" );
}

static int getNumElemPerBox( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Triangle:      return 2;
    case ElementType::Quadrilateral: return 1;
    case ElementType::Tetrahedron:   return 6;
    case ElementType::Wedge:         return 2;
    case ElementType::Pyramid:       return 6;
    case ElementType::Hexahedron:    return 1;
    default:
    {
      GEOS_THROW( "Unsupported element type " << elementType, InputError );
      return 0;
    }
  }
}

void InternalMeshGenerator::postInputInitialization()
{
  m_dim = getElementDim( EnumStrings< ElementType >::fromString( m_elementType[0] ) );

  // Should this have the size of the problem as enforced by m_dim?
  if( m_periodic.size() == 0 ){
    m_periodic.resize( 3 );
    LvArray::tensorOps::fill< 3 >(m_periodic, 0);
  }
  else
  {
    GEOS_ERROR_IF( m_periodic.size() !=3, "Periodic flags must have size 3" );
  }

  {
    // Check for vertex/element matching
    bool failFlag = false;
    for( int i=0; i<m_dim; ++i )
    {
      failFlag += ( m_nElems[i].size() != m_vertices[i].size()-1 );
    }
    if( failFlag )
    {
      GEOS_ERROR( getDataContext() << ": vertex/element mismatch.\n" << generalMeshErrorAdvice );
    }

    // If specified, check to make sure bias values have the correct length
    for( int i=0; i<m_dim; ++i )
    {
      if( m_nElemBias[i].size() > 0 )
      {
        failFlag += ( m_nElems[i].size() != m_nElemBias[i].size() );
      }
    }
    if( failFlag )
    {
      GEOS_ERROR( getDataContext() << ": element/bias mismatch.\n" << generalMeshErrorAdvice );
    }
  }

  m_numElePerBox.resize( m_nElems[0].size() * m_nElems[1].size() * m_nElems[2].size());

  if( m_elementType.size() != m_numElePerBox.size() )
  {
    if( m_elementType.size() == 1 )
    {
      string const elementType = m_elementType[0];
      m_elementType.resizeDefault( m_numElePerBox.size(), elementType );
    }
    else
    {
      GEOS_ERROR( getDataContext() << ": InternalMeshGenerator: The number of element types is inconsistent" <<
                  " with the number of total cell blocks.\n" << generalMeshErrorAdvice );
    }
  }

  for( localIndex i = 0; i < LvArray::integerConversion< localIndex >( m_elementType.size() ); ++i )
  {
    try
    {
      m_numElePerBox[i] = getNumElemPerBox( EnumStrings< ElementType >::fromString( m_elementType[i] ) );
    } catch( InputError const & e )
    {
      WrapperBase const & wrapper = getWrapperBase( viewKeyStruct::elementTypesString() );
      throw InputError( e, "InternalMesh " + wrapper.getDataContext().toString() +
                        ", element index = " + std::to_string( i ) + ": " );
    }
  }

  {
    localIndex numBlocks = 1;
    for( int i=0; i<m_dim; ++i )
    {
      numBlocks *= m_nElems[i].size();
    }
    if( numBlocks != m_regionNames.size() )
    {
      if( m_regionNames.size() == 1 )
      {
        string const regionName = m_regionNames[0];
        m_regionNames.resizeDefault( numBlocks, regionName );
      }
      else
      {
        GEOS_ERROR( getDataContext() << ": Incorrect number of regionLayout entries specified." );
      }
    }
  }

  for( int i=0; i<3; ++i )
  {
    m_min[i] = m_vertices[i].front();
    m_max[i] = m_vertices[i].back();
  }

  for( int dir=0; dir<3; ++dir )
  {
    m_firstElemIndexForBlock[dir].resize( m_nElems[dir].size() );
    m_lastElemIndexForBlock[dir].resize( m_nElems[dir].size() );
    m_firstElemIndexForBlock[dir][0] = 0;
    m_lastElemIndexForBlock[dir][0] = m_nElems[dir][0]-1;
    for( int block=1; block<m_nElems[dir].size(); ++block )
    {
      m_firstElemIndexForBlock[dir][block] = m_lastElemIndexForBlock[dir][block-1] + 1;
      m_lastElemIndexForBlock[dir][block] = m_firstElemIndexForBlock[dir][block] + m_nElems[dir][block]-1;
    }
  }

  m_fPerturb = 0.0;
}

/**
 * @brief Get the label mapping of element vertices indexes onto node indexes for a type of element.
 * @param[in] elementType the element type
 * @param[in] index ndim-spatialized Element index.
 * @param[in] iEle the index of Element begin processed
 * @param[out] nodeIDInBox array to map element vertices index to node indexes
 * @param[in] size the number of node on the element
 */
static void getElemToNodesRelationInBox( ElementType const elementType,
                                         integer const trianglePattern,
                                         int const (&index)[3],
                                         int const & iEle,
                                         int (& nodeIDInBox)[8],
                                         int const node_size )

{
  switch( elementType )
  {
    case ElementType::Hexahedron:
    {
      nodeIDInBox[0] = 0;
      nodeIDInBox[1] = 1;
      nodeIDInBox[2] = 3;
      nodeIDInBox[3] = 2;
      nodeIDInBox[4] = 4;
      nodeIDInBox[5] = 5;
      nodeIDInBox[6] = 7;
      nodeIDInBox[7] = 6;
      break;
    }
    case ElementType::Wedge:
    {
      if( trianglePattern == 0 )
      {
        if( ( index[0] + index[1] ) % 2 == 1 )
        {
          if( iEle % 2 == 0 )
          {
            nodeIDInBox[0] = 6;
            nodeIDInBox[1] = 2;
            nodeIDInBox[2] = 4;
            nodeIDInBox[3] = 0;
            nodeIDInBox[4] = 7;
            nodeIDInBox[5] = 3;
            nodeIDInBox[6] = 7;
            nodeIDInBox[7] = 3;
          }
          else
          {
            nodeIDInBox[0] = 5;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 4;
            nodeIDInBox[3] = 0;
            nodeIDInBox[4] = 6;
            nodeIDInBox[5] = 2;
            nodeIDInBox[6] = 6;
            nodeIDInBox[7] = 2;
          }
        }
        else
        {
          if( iEle % 2 == 0 )
          {
            nodeIDInBox[0] = 5;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 4;
            nodeIDInBox[3] = 0;
            nodeIDInBox[4] = 7;
            nodeIDInBox[5] = 3;
            nodeIDInBox[6] = 7;
            nodeIDInBox[7] = 3;
          }
          else
          {
            nodeIDInBox[0] = 5;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 7;
            nodeIDInBox[3] = 3;
            nodeIDInBox[4] = 6;
            nodeIDInBox[5] = 2;
            nodeIDInBox[6] = 6;
            nodeIDInBox[7] = 2;
          }
        }
      }
      else if( trianglePattern == 1 )
      {
        if( index[1] % 2 == 0 )
        {
          if( iEle % 2 == 0 )
          {
            nodeIDInBox[0] = 5;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 4;
            nodeIDInBox[3] = 0;
            nodeIDInBox[4] = 6;
            nodeIDInBox[5] = 2;
            nodeIDInBox[6] = 6;
            nodeIDInBox[7] = 2;
          }
          else
          {
            nodeIDInBox[0] = 6;
            nodeIDInBox[1] = 2;
            nodeIDInBox[2] = 4;
            nodeIDInBox[3] = 0;
            nodeIDInBox[4] = 7;
            nodeIDInBox[5] = 3;
            nodeIDInBox[6] = 7;
            nodeIDInBox[7] = 3;
          }
        }
        else
        {
          if( iEle % 2 == 0 )
          {
            nodeIDInBox[0] = 5;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 4;
            nodeIDInBox[3] = 0;
            nodeIDInBox[4] = 7;
            nodeIDInBox[5] = 3;
            nodeIDInBox[6] = 7;
            nodeIDInBox[7] = 3;
          }
          else
          {
            nodeIDInBox[0] = 5;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 7;
            nodeIDInBox[3] = 3;
            nodeIDInBox[4] = 6;
            nodeIDInBox[5] = 2;
            nodeIDInBox[6] = 6;
            nodeIDInBox[7] = 2;
          }
        }
      }
      break;
    }
    case ElementType::Quadrilateral:
    {
      nodeIDInBox[0] = 0;
      nodeIDInBox[1] = 1;
      nodeIDInBox[2] = 3;
      nodeIDInBox[3] = 2;
      break;
    }
    case ElementType::Triangle:
    {
      if( trianglePattern == 0 )
      {
        if( ( index[0] + index[1] ) % 2 == 1 )
        {
          if( iEle % 2 == 0 )
          {
            nodeIDInBox[0] = 0;
            nodeIDInBox[1] = 2;
            nodeIDInBox[2] = 3;
          }
          else
          {
            nodeIDInBox[0] = 0;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 2;
          }
        }
        else
        {
          if( iEle % 2 == 0 )
          {
            nodeIDInBox[0] = 0;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 3;
          }
          else
          {
            nodeIDInBox[0] = 1;
            nodeIDInBox[1] = 2;
            nodeIDInBox[2] = 3;
          }
        }
      }
      else if( trianglePattern == 1 )
      {
        if( index[1] % 2 == 0 )
        {
          if( iEle % 2 == 0 )
          {
            nodeIDInBox[0] = 0;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 2;
          }
          else
          {
            nodeIDInBox[0] = 2;
            nodeIDInBox[1] = 3;
            nodeIDInBox[2] = 0;
          }
        }
        else
        {
          if( iEle % 2 == 0 )
          {
            nodeIDInBox[0] = 0;
            nodeIDInBox[1] = 1;
            nodeIDInBox[2] = 3;
          }
          else
          {
            nodeIDInBox[0] = 1;
            nodeIDInBox[1] = 2;
            nodeIDInBox[2] = 3;
          }
        }
      }
      break;
    }
    case ElementType::Tetrahedron:
    {
      int mapBoxTet[8][6][4] =
      {
        {
          { 0, 3, 7, 6 },
          { 0, 7, 4, 6 },
          { 0, 5, 1, 6 },
          { 0, 4, 5, 6 },
          { 0, 1, 2, 6 },
          { 0, 2, 3, 6 }
        },
        {
          { 1, 5, 6, 7 },
          { 1, 6, 2, 7 },
          { 0, 4, 1, 7 },
          { 1, 4, 5, 7 },
          { 0, 1, 3, 7 },
          { 1, 2, 3, 7 }
        },
        {
          { 0, 3, 4, 2 },
          { 3, 7, 4, 2 },
          { 0, 4, 1, 2 },
          { 1, 4, 5, 2 },
          { 4, 7, 6, 2 },
          { 4, 6, 5, 2 }
        },
        {
          { 1, 5, 2, 3 },
          { 2, 5, 6, 3 },
          { 0, 5, 1, 3 },
          { 0, 4, 5, 3 },
          { 4, 7, 5, 3 },
          { 5, 7, 6, 3 }
        },
        {
          { 0, 3, 4, 5 },
          { 3, 7, 4, 5 },
          { 3, 6, 7, 5 },
          { 3, 2, 6, 5 },
          { 3, 0, 1, 5 },
          { 1, 2, 3, 5 }
        },
        {
          { 1, 5, 2, 4 },
          { 2, 5, 6, 4 },
          { 2, 6, 7, 4 },
          { 2, 7, 3, 4 },
          { 0, 2, 3, 4 },
          { 0, 1, 2, 4 }
        },
        {
          { 0, 7, 4, 1 },
          { 0, 3, 7, 1 },
          { 2, 7, 3, 1 },
          { 2, 6, 7, 1 },
          { 4, 7, 5, 1 },
          { 5, 7, 6, 1 }
        },
        {
          { 1, 5, 6, 0 },
          { 1, 6, 2, 0 },
          { 2, 6, 3, 0 },
          { 3, 6, 7, 0 },
          { 4, 6, 5, 0 },
          { 4, 7, 6, 0 }
        }
      };

      int mapBoxType[2][2][2];
      mapBoxType[0][0][0] = 0;
      mapBoxType[1][0][0] = 1;
      mapBoxType[0][0][1] = 2;
      mapBoxType[1][0][1] = 3;
      mapBoxType[0][1][0] = 4;
      mapBoxType[1][1][0] = 5;
      mapBoxType[0][1][1] = 6;
      mapBoxType[1][1][1] = 7;

      int boxType = mapBoxType[index[0] % 2][index[1] % 2][index[2] % 2];
      for( int i = 0; i < node_size; ++i )
      {
        nodeIDInBox[i] = mapBoxTet[boxType][iEle][i];
      }
      break;
    }
    default:
    {
      GEOS_ERROR( "InternalMeshGenerator: unsupported element type " << elementType << ".\n" << generalMeshErrorAdvice );
    }
  }
}

void InternalMeshGenerator::fillCellBlockManager( CellBlockManager & cellBlockManager, array1d< int > const & partition )
{
  GEOS_MARK_FUNCTION;
  
  //Needs to be set before addNeighbors call that occurs in SpatialPartition::setSizes
  m_partition.setPeriodic( m_periodic );

  m_partition.setPartitions( partition );
  // Partition based on even spacing to get load balance
  // Partition geometrical boundaries will be corrected in the end.
  {
    m_min[0] = m_vertices[0].front();
    m_min[1] = m_vertices[1].front();
    m_min[2] = m_vertices[2].front();

    m_max[0] = m_vertices[0].back();
    m_max[1] = m_vertices[1].back();
    m_max[2] = m_vertices[2].back();

    m_partition.setSizes( m_min, m_max );
  }

  // Make sure that the node manager fields are initialized
  auto & nodeSets = cellBlockManager.getNodeSets();

  real64 size[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_max );
  LvArray::tensorOps::subtract< 3 >( size, m_min );
  cellBlockManager.setGlobalLength( LvArray::tensorOps::l2Norm< 3 >( size ) );

//  bool isRadialWithOneThetaPartition = false;

  // This should probably handled elsewhere:
  int aa = 0;
  for( auto & cellBlockName : m_regionNames )
  {
    CellBlock & cellBlock = cellBlockManager.registerCellBlock( cellBlockName );
    cellBlock.setElementType( EnumStrings< ElementType >::fromString( m_elementType[aa++] ) );
  }

  SortedArray< localIndex > & xnegNodes = nodeSets["xneg"];
  SortedArray< localIndex > & xposNodes = nodeSets["xpos"];
  SortedArray< localIndex > & ynegNodes = nodeSets["yneg"];
  SortedArray< localIndex > & yposNodes = nodeSets["ypos"];
  SortedArray< localIndex > & znegNodes = nodeSets["zneg"];
  SortedArray< localIndex > & zposNodes = nodeSets["zpos"];
  SortedArray< localIndex > & allNodes = nodeSets["all"];

  // Find elemCenters for even uniform element sizes
  array1d< array1d< real64 > > elemCenterCoords( 3 );
  for( int dim = 0; dim < m_dim; ++dim )
  {
    m_numElemsTotal[dim] = 0;
    for( int block = 0; block < m_nElems[dim].size(); ++block )
    {
      m_numElemsTotal[dim] += m_nElems[dim][block];
    }
    array1d< int > const parts = m_partition.getPartitions();
    GEOS_ERROR_IF( parts[dim] > m_numElemsTotal[dim], "Number of partitions in a direction should not exceed the number of elements in that direction" );

    elemCenterCoords[dim].resize( m_numElemsTotal[dim] );
    array1d< real64 > elemCenterCoordsLocal( m_numElemsTotal[dim] );
    for( integer k = 0; k < m_numElemsTotal[dim]; ++k )
    {
      elemCenterCoordsLocal[k] = m_min[dim] + ( m_max[dim] - m_min[dim] ) * ( k + 0.5 ) / m_numElemsTotal[dim];
    }
    MpiWrapper::allReduce( elemCenterCoordsLocal.data(),
                           elemCenterCoords[dim].data(),
                           m_numElemsTotal[dim],
                           MPI_MAX,
                           MPI_COMM_GEOSX );
  }

  // Find starting/ending index
  // Get the first and last indices in this partition each direction
  integer firstElemIndexInPartition[3] = { -1, -1, -1 };
  integer lastElemIndexInPartition[3] = { -2, -2, -2 };

  for( int dim = 0; dim < m_dim; ++dim )
  {
    //    firstElemIndexInPartition[i] = -1;
    //    lastElemIndexInPartition[i] = -2;
    for( int k = 0; k < m_numElemsTotal[dim]; ++k )
    {
      if( m_partition.isCoordInPartition( elemCenterCoords[dim][k], dim ) )
      {
        firstElemIndexInPartition[dim] = k;
        break;
      }
    }

    if( firstElemIndexInPartition[dim] > -1 )
    {
      for( int k = firstElemIndexInPartition[dim]; k < m_numElemsTotal[dim]; ++k )
      {
        if( m_partition.isCoordInPartition( elemCenterCoords[dim][k], dim ) )
        {
          lastElemIndexInPartition[dim] = k;
        }
      }
    }
  }

  // Calculate number of elements in this partition from each region, and the
  // total number of nodes

  std::map< string, int > numElemsInRegions;
  std::map< string, ElementType > elemTypeInRegions;

  array1d< integer > firstElemIndexForBlockInPartition[3];
  array1d< integer > lastElemIndexForBlockInPartition[3];

  for( int dim = 0; dim < 3; ++dim )
  {
    firstElemIndexForBlockInPartition[dim] = m_firstElemIndexForBlock[dim];
    lastElemIndexForBlockInPartition[dim] = m_lastElemIndexForBlock[dim];

    for( integer block = 0; block < m_nElems[dim].size(); ++block )
    {
      if( firstElemIndexForBlockInPartition[dim][block] > lastElemIndexInPartition[dim] ||
          lastElemIndexForBlockInPartition[dim][block] < firstElemIndexInPartition[dim] )
      {
        firstElemIndexForBlockInPartition[dim][block] = -1;
        lastElemIndexForBlockInPartition[dim][block] = -2;
      }
      else
      {
        if( firstElemIndexForBlockInPartition[dim][block] < firstElemIndexInPartition[dim] )
        {
          firstElemIndexForBlockInPartition[dim][block] = firstElemIndexInPartition[dim];
        }
        if( lastElemIndexForBlockInPartition[dim][block] > lastElemIndexInPartition[dim] )
        {
          lastElemIndexForBlockInPartition[dim][block] = lastElemIndexInPartition[dim];
        }
      }
    }
  }

  // TODO This needs to be rewritten for dimensions lower than 3.
  localIndex regionOffset = 0;
  for( integer kblock = 0; kblock < m_nElems[2].size(); ++kblock )
  {
    for( integer jblock = 0; jblock < m_nElems[1].size(); ++jblock )
    {
      for( integer iblock = 0; iblock < m_nElems[0].size(); ++iblock, ++regionOffset )
      {
        numElemsInRegions[ m_regionNames[ regionOffset ] ] = 0;
        elemTypeInRegions[ m_regionNames[ regionOffset ] ] = ElementType::Quadrilateral;
      }
    }
  }

  regionOffset = 0;
  {
    localIndex iR = 0;
    for( integer kblock = 0; kblock < m_nElems[2].size(); ++kblock )
    {
      for( integer jblock = 0; jblock < m_nElems[1].size(); ++jblock )
      {
        for( integer iblock = 0; iblock < m_nElems[0].size(); ++iblock, ++regionOffset, ++iR )
        {
          integer numElemsInRegion = 1;
          numElemsInRegion *= lastElemIndexForBlockInPartition[0][iblock] - firstElemIndexForBlockInPartition[0][iblock] + 1;

          if( m_dim > 1 )
          {
            numElemsInRegion *= lastElemIndexForBlockInPartition[1][jblock] - firstElemIndexForBlockInPartition[1][jblock] + 1;
          }
          if( m_dim > 2 )
          {
            numElemsInRegion *= lastElemIndexForBlockInPartition[2][kblock] - firstElemIndexForBlockInPartition[2][kblock] + 1;
          }

          numElemsInRegion *= m_numElePerBox[iR];
          numElemsInRegions[ m_regionNames[ regionOffset ] ] += numElemsInRegion;
          elemTypeInRegions[ m_regionNames[ regionOffset ] ] = EnumStrings< ElementType >::fromString( m_elementType[iR] );
        }
      }
    }
  }

  localIndex numNodes = 1;
  integer numNodesInDir[3] = { 1, 1, 1 };
  for( int dim = 0; dim < m_dim; ++dim )
  {
    numNodesInDir[dim] = lastElemIndexInPartition[dim] - firstElemIndexInPartition[dim] + 2;
  }
  reduceNumNodesForPeriodicBoundary( m_partition, numNodesInDir );
  numNodes = numNodesInDir[0] * numNodesInDir[1] * numNodesInDir[2];

  cellBlockManager.setNumNodes( numNodes );

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > X = cellBlockManager.getNodePositions();

  arrayView1d< globalIndex > const nodeLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();

  {
    localIndex localNodeIndex = 0;
    for( integer k = 0; k < numNodesInDir[2]; ++k )
    {
      for( integer j = 0; j < numNodesInDir[1]; ++j )
      {
        for( integer i = 0; i < numNodesInDir[0]; ++i )
        {
          integer globalIJK[3] = { i, j, k };

          for( int dim = 0; dim < m_dim; ++dim )
          {
            globalIJK[dim] += firstElemIndexInPartition[dim];
          }

          getNodePosition( globalIJK, m_trianglePattern, X[localNodeIndex] );

          // Alter global node map for radial mesh
          setNodeGlobalIndicesOnPeriodicBoundary( globalIJK );

          nodeLocalToGlobal[localNodeIndex] = nodeGlobalIndex( globalIJK );

          // Cartesian-specific nodesets
          if( isCartesian() )
          {
            if( isEqual( X( localNodeIndex, 0 ), m_min[0], m_coordinatePrecision ) )
            {
              xnegNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 0 ), m_max[0], m_coordinatePrecision ) )
            {
              xposNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 1 ), m_min[1], m_coordinatePrecision ) )
            {
              ynegNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 1 ), m_max[1], m_coordinatePrecision ) )
            {
              yposNodes.insert( localNodeIndex );
            }
          }

          // General nodesets
          if( isEqual( X( localNodeIndex, 2 ), m_min[2], m_coordinatePrecision ) )
          {
            znegNodes.insert( localNodeIndex );
          }
          if( isEqual( X( localNodeIndex, 2 ), m_max[2], m_coordinatePrecision ) )
          {
            zposNodes.insert( localNodeIndex );
          }
          allNodes.insert( localNodeIndex );

          ++localNodeIndex;
        }
      }
    }
  }

  {
    array1d< integer > numElements;
    array1d< string > elementRegionNames;
    std::map< string, localIndex > localElemIndexInRegion;

    for( auto const & numElemsInRegion : numElemsInRegions )
    {
      numElements.emplace_back( numElemsInRegion.second );
      elementRegionNames.emplace_back( numElemsInRegion.first );
      localElemIndexInRegion[numElemsInRegion.first] = 0;
    }

    cellBlockManager.resize( numElements, elementRegionNames );

    // Assign global numbers to elementsl
    regionOffset = 0;
    SortedArray< string > processedRegionNames;
    localIndex iR = 0;

    // Reset the number of nodes in each dimension in case of periodic BCs so the element firstNodeIndex
    //  calculation is correct? Not actually needed in parallel since we still have ghost nodes in that case and
    //  the count has not been altered due to periodicity.
    auto const periodic = m_partition.getPeriodic();
    if( std::any_of( periodic.begin(), periodic.end(), []( int & dimPeriodic ) { return dimPeriodic == 1; } ) )
    {
      for( int i = 0; i < m_dim; ++i )
      {
        numNodesInDir[i] = lastElemIndexInPartition[i] - firstElemIndexInPartition[i] + 2;
      }
      numNodes = numNodesInDir[0] * numNodesInDir[1] * numNodesInDir[2];
    }

    for( integer kblock = 0; kblock < m_nElems[2].size(); ++kblock )
    {
      for( integer jblock = 0; jblock < m_nElems[1].size(); ++jblock )
      {
        for( integer iblock = 0; iblock < m_nElems[0].size(); ++iblock, ++regionOffset, ++iR )
        {
          ElementType const elementType = EnumStrings< ElementType >::fromString( m_elementType[iR] );

          CellBlock & cellBlock = cellBlockManager.getCellBlock( m_regionNames[regionOffset] );
          int const numNodesPerElem = LvArray::integerConversion< int >( cellBlock.numNodesPerElement());
          integer nodeIDInBox[ 8 ];

          arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodes = cellBlock.getElemToNode();
          arrayView1d< globalIndex > const & elemLocalToGlobal = cellBlock.localToGlobalMap();

          integer numElemsInDirForBlock[3] =
          { lastElemIndexForBlockInPartition[0][iblock] - firstElemIndexForBlockInPartition[0][iblock] + 1,
            lastElemIndexForBlockInPartition[1][jblock] - firstElemIndexForBlockInPartition[1][jblock] + 1,
            lastElemIndexForBlockInPartition[2][kblock] - firstElemIndexForBlockInPartition[2][kblock] + 1 };

          for( integer k = 0; k < numElemsInDirForBlock[2]; ++k )
          {
            for( integer j = 0; j < numElemsInDirForBlock[1]; ++j )
            {
              for( integer i = 0; i < numElemsInDirForBlock[0]; ++i )
              {
                integer globalIJK[3] =
                { i + firstElemIndexForBlockInPartition[0][iblock],
                  j + firstElemIndexForBlockInPartition[1][jblock],
                  k + firstElemIndexForBlockInPartition[2][kblock] };

                localIndex const firstNodeIndex = ( globalIJK[0] - firstElemIndexInPartition[0] )
                                                  + numNodesInDir[0] * ( globalIJK[1] - firstElemIndexInPartition[1] )
                                                  + numNodesInDir[0] * numNodesInDir[1] * ( globalIJK[2] - firstElemIndexInPartition[2] );
                localIndex nodeOfBox[8];

                if( elementType == ElementType::Quadrilateral || elementType == ElementType::Triangle )
                {
                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;
                }
                else
                {
                  localIndex const stride[3] = { 1, numNodesInDir[0], numNodesInDir[0] * numNodesInDir[1] };

                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] = nodeOfBox[0] + stride[0];
                  nodeOfBox[2] = nodeOfBox[1] + stride[1];
                  nodeOfBox[3] = nodeOfBox[0] + stride[1];

                  nodeOfBox[4] = nodeOfBox[0] + stride[2];
                  nodeOfBox[5] = nodeOfBox[1] + stride[2];
                  nodeOfBox[6] = nodeOfBox[2] + stride[2];
                  nodeOfBox[7] = nodeOfBox[3] + stride[2];

                  //               7___________________ 6
                  //               /                   /|
                  //              /                   / |
                  //             /                   /  |
                  //           4/__________________5/   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |   3               |   /2        z
                  //            |                   |  /          |   y
                  //            |                   | /           |  /
                  //            |___________________|/            | /
                  //            0                   1             |/____ x

                }


                // Fix local connectivity for single theta (y) partition (radial meshes only)

                setConnectivityForPeriodicBoundaries( globalIJK,
                                                      numNodesInDir,
                                                      firstElemIndexInPartition,
                                                      nodeOfBox );

                for( int iEle = 0; iEle < m_numElePerBox[iR]; ++iEle )
                {
                  localIndex & localElemIndex = localElemIndexInRegion[ m_regionNames[ regionOffset ] ];
                  elemLocalToGlobal[localElemIndex] = elemGlobalIndex( globalIJK ) * m_numElePerBox[iR] + iEle;

                  getElemToNodesRelationInBox( elementType,
                                               m_trianglePattern,
                                               globalIJK,
                                               iEle,
                                               nodeIDInBox,
                                               numNodesPerElem );

                  for( localIndex iN = 0; iN < numNodesPerElem; ++iN )
                  {
                    elemsToNodes[localElemIndex][iN] = nodeOfBox[nodeIDInBox[iN]];
                  }
                  ++localElemIndex;
                }
              }
            }
          }
        }
      }
    }
  }

  // Node perturbation
  if( m_fPerturb > 0 )
  {

    for( localIndex iN = 0; iN != numNodes; ++iN )
    {

      for( int i = 0; i < m_dim; ++i )
      {
        if( X[iN][i] > m_min[i] && X[iN][i] < m_max[i] )
        {
          // This ensures that the perturbation pattern is unaffected by domain
          srand( LvArray::integerConversion< int >( nodeLocalToGlobal[iN] ) + m_randSeed + i );

          X[iN][i] += ( ( m_max[i] - m_min[i] ) / m_numElemsTotal[i] ) * ( ( rand() * 1.0 ) / RAND_MAX - 0.5 ) * 2 * m_fPerturb;
        }
      }
    }
  }

  if( std::fabs( m_skewAngle ) > 0.0 )
  {
    for( localIndex iN = 0; iN != numNodes; ++iN )
    {
      X[iN][0] -= ( X[iN][1] - m_skewCenter[1] ) * std::tan( m_skewAngle );
    }
  }

  coordinateTransformation( X, nodeSets );

  cellBlockManager.buildMaps();

  GEOS_LOG_RANK_0( GEOS_FMT( "{}: total number of nodes = {}", getName(),
                             ( m_numElemsTotal[0] + 1 ) * ( m_numElemsTotal[1] + 1 ) * ( m_numElemsTotal[2] + 1 ) ) );
  GEOS_LOG_RANK_0( GEOS_FMT( "{}: total number of elems = {}", getName(),
                             m_numElemsTotal[0] * m_numElemsTotal[1] * m_numElemsTotal[2] ) );
}

void
InternalMeshGenerator::
  setConnectivityForPeriodicBoundary( int const component,
                                      int const (&globalIJK)[3],
                                      integer const (&numNodesInDir)[3],
                                      int const (&firstElemIndexInPartition)[3],
                                      localIndex (& nodeOfBox)[8] )
{
  GEOS_ERROR_IF( component != 1, "Connectivity for periodic boundary implemented for the 1-component only" );

  // Condition is:
  // 1) element is last index in component direction
  // 2) first local element in component partition is zero
  if( ( globalIJK[component] == m_numElemsTotal[component] - 1 ) &&
      ( firstElemIndexInPartition[component] == 0) )
  {
    // Last set of nodes
    int modGlobalIJK[3] = { globalIJK[0], globalIJK[1], globalIJK[2] };
    modGlobalIJK[component] = 0;
    localIndex const firstNodeIndex = ( modGlobalIJK[0] - firstElemIndexInPartition[0] )
                                      + numNodesInDir[0] * ( modGlobalIJK[1] - firstElemIndexInPartition[1] )
                                      + numNodesInDir[0] * numNodesInDir[1] * ( modGlobalIJK[2] - firstElemIndexInPartition[2] );

    nodeOfBox[3] = firstNodeIndex;
    nodeOfBox[2] = nodeOfBox[3] + 1;
    nodeOfBox[7] = nodeOfBox[3] + numNodesInDir[0] * numNodesInDir[1];
    nodeOfBox[6] = nodeOfBox[7] + 1;
  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, InternalMeshGenerator, string const &, Group * const )
} /* namespace geos */
