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
 * @file LAIHelperFunctions.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_
#define GEOS_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geos
{
namespace LAIHelperFunctions
{

/**
 * @brief Create an identity matrix.
 * @tparam MATRIX type of matrix
 * @param n local size of the square identity matrix
 * @param comm MPI communicator
 * @param mat the output matrix
 */
template< typename MATRIX >
void makeIdentity( localIndex const n,
                   MPI_Comm const & comm,
                   MATRIX & mat )
{
  mat.createWithLocalSize( n, 1, comm );
  mat.open();
  for( globalIndex i = mat.ilower(); i < mat.iupper(); ++i )
  {
    mat.insert( i, i, 1.0 );
  }
  mat.close();
}

/**
 * @brief Create a permutation matrix for a given nodal variable.
 * @tparam     MATRIX            the parallel matrix type
 * @param[in]  nodeManager       the node manager
 * @param[in]  nDofPerNode       number of degrees-of-freedom per node
 * @param[in]  dofKey            DofManager key used to access dof index array
 * @param[out] permutationMatrix the target matrix
 */
template< typename MATRIX >
void createPermutationMatrix( NodeManager const & nodeManager,
                              int const nDofPerNode,
                              string const & dofKey,
                              MATRIX & permutationMatrix )
{
  /* Crates a permutation matrix for a given nodal variable specified by the DofKey.
   * It consider that nDofPerNode dofs are associated with each node (e.g., nDofPerNode = 3 for the displacement).
   *
   * The permutation matrix maps from the dofs ordering set by the DOFManager to the ordering based on the global
   * indexes of the nodes.
   */

  localIndex const numLocalRows = nodeManager.getNumberOfLocalIndices() * nDofPerNode;
  permutationMatrix.createWithLocalSize( numLocalRows, numLocalRows, 1, MPI_COMM_GEOSX );

  arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
  arrayView1d< globalIndex const > const & localToGlobal = nodeManager.localToGlobalMap();

  permutationMatrix.open();
  for( localIndex a = 0; a < nodeManager.size(); ++a )
  {
    if( dofNumber[a] >= 0 )
    {
      for( int d = 0; d < nDofPerNode; ++d )
      {
        globalIndex const rowIndex    = localToGlobal[a] * nDofPerNode + d;
        globalIndex const columnIndex = dofNumber[a] + d;
        permutationMatrix.insert( rowIndex, columnIndex, 1.0 );
      }
    }
  }
  permutationMatrix.close();
  permutationMatrix.set( 1.0 );
}

/**
 * @brief Create a permutation matrix for a given cell-centered variable.
 * @tparam     MATRIX            the parallel matrix type
 * @param[in]  elemManager       the element region manager
 * @param[in]  nDofPerCell       number of degrees-of-freedom per node
 * @param[in]  dofKey            DofManager key used to access dof index array
 * @param[out] permutationMatrix the target matrix
 */
template< typename MATRIX >
void createPermutationMatrix( ElementRegionManager const & elemManager,
                              int const nDofPerCell,
                              string const & dofKey,
                              MATRIX & permutationMatrix )
{
  /* Crates a permutation matrix for a given cell centered variable specified by the DofKey.
   * It consider that nDofPerCell dofs are associated with each cell (e.g., nDofPerCell = 3 for the displacement).
   *
   * The permutation matrix maps from the dofs ordering set by the DOFManager to the ordering based on the global
   * indexes of the cells.
   */

  // Create permutation matrix
  localIndex numLocalRows = 0;
  elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & elementSubRegion )
  {
    if( elementSubRegion.hasWrapper( dofKey ) )
    {
      numLocalRows += elementSubRegion.getNumberOfLocalIndices() * nDofPerCell;
    }
  } );
  permutationMatrix.createWithLocalSize( numLocalRows, numLocalRows, 1, MPI_COMM_GEOSX );

  permutationMatrix.open();
  elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();
    arrayView1d< globalIndex const > const & dofNumber = elementSubRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< globalIndex const > const & localToGlobal = elementSubRegion.localToGlobalMap();

    for( localIndex k = 0; k < numElems; ++k )
    {
      if( dofNumber[k] >= 0 )
      {
        for( int d = 0; d < nDofPerCell; ++d )
        {
          globalIndex const rowIndex    = localToGlobal[k] * nDofPerCell + d;
          globalIndex const columnIndex = dofNumber[k] + d;
          permutationMatrix.insert( rowIndex, columnIndex, 1.0 );
        }
      }
    }
  } );
  permutationMatrix.close();
  permutationMatrix.set( 1.0 );
}

/**
 * @brief Permute a vector.
 * @tparam    VECTOR            the parallel vector type
 * @tparam    MATRIX            the parallel matrix type
 * @param[in] vector            the source vector
 * @param[in] permutationMatrix the permutation matrix
 * @return the permuted vector
 */
template< typename VECTOR, typename MATRIX >
VECTOR permuteVector( VECTOR const & vector,
                      MATRIX const & permutationMatrix )
{
  VECTOR permutedVector;
  permutedVector.create( vector.localSize(), permutationMatrix.comm() );
  permutationMatrix.apply( vector, permutedVector );
  return permutedVector;
}

/**
 * @brief Permute rows and columns of a square matrix.
 * @param[in] matrix            the source matrix
 * @param[in] permutationMatrix permutation matrix
 * @return the permuted matrix
 */
template< typename MATRIX >
MATRIX permuteMatrix( MATRIX const & matrix,
                      MATRIX const & permutationMatrix )
{
  MATRIX permutedMatrix;
  matrix.multiplyRARt( permutationMatrix, permutedMatrix );
  return matrix;
}

/**
 * Permute rows and columns of a rectangular matrix
 * @param[in] matrix                 the source matrix
 * @param[in] permutationMatrixLeft  left permutation matrix
 * @param[in] permutationMatrixRight right permutation matrix
 * @return the permuted matrix
 */
template< typename MATRIX >
MATRIX permuteMatrix( MATRIX const & matrix,
                      MATRIX const & permutationMatrixLeft,
                      MATRIX const & permutationMatrixRight )
{
  MATRIX permutedMatrix;
  matrix.multiplyRAP( permutationMatrixLeft, permutationMatrixRight, permutedMatrix );
  return matrix;
}

/**
 * @brief Computes rigid body modes
 * @tparam VECTOR output vector type
 * @param mesh the mesh
 * @param dofManager the degree-of-freedom manager
 * @param selection list of field names
 * @param rigidBodyModes the output array of linear algebra vectors containing RBMs
 */
template< typename VECTOR >
void computeRigidBodyModes( MeshLevel const & mesh,
                            DofManager const & dofManager,
                            std::vector< string > const & selection,
                            array1d< VECTOR > & rigidBodyModes )
{
  NodeManager const & nodeManager = mesh.getNodeManager();

  localIndex numComponents = 0;
  array1d< localIndex > globalNodeList;
  for( localIndex k = 0; k < LvArray::integerConversion< localIndex >( selection.size() ); ++k )
  {
    if( dofManager.location( selection[k] ) ==  FieldLocation::Node )
    {
      string const & dispDofKey = dofManager.getKey( selection[k] );
      arrayView1d< globalIndex const > const & dofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
      localIndex const numComponentsField = dofManager.numComponents( selection[k] );
      numComponents = numComponents > 0 ? numComponents : numComponentsField;
      GEOS_ERROR_IF( numComponents != numComponentsField, "Rigid body modes called with different number of components." );
      globalIndex const globalOffset = dofManager.globalOffset( selection[k] );
      globalIndex const numLocalDofs = LvArray::integerConversion< globalIndex >( dofManager.numLocalDofs( selection[k] ) );
      for( globalIndex i = 0; i < dofNumber.size(); ++i )
      {
        if( dofNumber[i] >= globalOffset && ( dofNumber[i] - globalOffset ) < numLocalDofs )
        {
          globalNodeList.emplace_back( LvArray::integerConversion< localIndex >( dofNumber[i] - globalOffset ) / numComponentsField );
        }
      }
    }
  }
  localIndex const numNodes = globalNodeList.size();
  arrayView1d< localIndex const > globalNodeListView = globalNodeList.toViewConst();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  localIndex const numRidigBodyModes = numComponents * ( numComponents + 1 ) / 2;
  rigidBodyModes.resize( numRidigBodyModes );
  for( localIndex k = 0; k < numComponents; ++k )
  {
    rigidBodyModes[k].create( numNodes * numComponents, MPI_COMM_GEOSX );
    arrayView1d< real64 > const values = rigidBodyModes[k].open();
    forAll< parallelHostPolicy >( numNodes, [=]( localIndex const i )
    {
      values[numComponents * i + k] = 1.0;
    } );
    rigidBodyModes[k].close();
    rigidBodyModes[k].scale( 1.0 / rigidBodyModes[k].norm2() );
  }
  switch( numComponents )
  {
    case 2:
    {
      localIndex const k = 2;
      rigidBodyModes[k].create( numNodes*numComponents, MPI_COMM_GEOSX );
      {
        arrayView1d< real64 > const values = rigidBodyModes[k].open();
        forAll< parallelHostPolicy >( numNodes, [=]( localIndex const i )
        {
          values[numComponents * i + 0] = -nodePosition[globalNodeListView[i]][1];
          values[numComponents * i + 1] = +nodePosition[globalNodeListView[i]][0];
        } );
        rigidBodyModes[k].close();
      }

      for( localIndex j = 0; j < k; ++j )
      {
        rigidBodyModes[k].axpy( -rigidBodyModes[k].dot( rigidBodyModes[j] ), rigidBodyModes[j] );
      }
      rigidBodyModes[k].scale( 1.0 / rigidBodyModes[k].norm2() );
      break;
    }
    case 3:
    {
      localIndex k = 3;
      rigidBodyModes[k].create( numNodes*numComponents, MPI_COMM_GEOSX );
      {
        arrayView1d< real64 > const values = rigidBodyModes[k].open();
        forAll< parallelHostPolicy >( numNodes, [=]( localIndex const i )
        {
          values[numComponents * i + 0] = +nodePosition[globalNodeListView[i]][1];
          values[numComponents * i + 1] = -nodePosition[globalNodeListView[i]][0];
        } );
        rigidBodyModes[k].close();
      }

      for( localIndex j = 0; j < k; ++j )
      {
        rigidBodyModes[k].axpy( -rigidBodyModes[k].dot( rigidBodyModes[j] ), rigidBodyModes[j] );
      }
      rigidBodyModes[k].scale( 1.0 / rigidBodyModes[k].norm2() );

      ++k;
      rigidBodyModes[k].create( numNodes*numComponents, MPI_COMM_GEOSX );
      {
        arrayView1d< real64 > const values = rigidBodyModes[k].open();
        forAll< parallelHostPolicy >( numNodes, [=]( localIndex const i )
        {
          values[numComponents * i + 1] = -nodePosition[globalNodeListView[i]][2];
          values[numComponents * i + 2] = +nodePosition[globalNodeListView[i]][1];
        } );
        rigidBodyModes[k].close();
      }

      for( localIndex j = 0; j < k; ++j )
      {
        rigidBodyModes[k].axpy( -rigidBodyModes[k].dot( rigidBodyModes[j] ), rigidBodyModes[j] );
      }
      rigidBodyModes[k].scale( 1.0 / rigidBodyModes[k].norm2() );

      ++k;
      rigidBodyModes[k].create( numNodes*numComponents, MPI_COMM_GEOSX );
      {
        arrayView1d< real64 > const values = rigidBodyModes[k].open();
        forAll< parallelHostPolicy >( numNodes, [=]( localIndex const i )
        {
          values[numComponents * i + 0] = +nodePosition[globalNodeListView[i]][2];
          values[numComponents * i + 2] = -nodePosition[globalNodeListView[i]][0];
        } );
        rigidBodyModes[k].close();
      }

      for( localIndex j = 0; j < k; ++j )
      {
        rigidBodyModes[k].axpy( -rigidBodyModes[k].dot( rigidBodyModes[j] ), rigidBodyModes[j] );
      }
      rigidBodyModes[k].scale( 1.0 / rigidBodyModes[k].norm2() );
      break;
    }
    default:
    {
      GEOS_ERROR( "Rigid body modes computation unsupported for " << numComponents << " components." );
    }
  }
}

} // LAIHelperFunctions namespace

} // geos namespace

#endif /*GEOS_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_*/
