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
 * @file testLinearAlgebraUtils.hpp
 */

#ifndef GEOS_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP_
#define GEOS_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP_

#include "common/DataTypes.hpp"
#include "common/initializeEnvironment.hpp"
#include "common/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshManager.hpp"


#include <gtest/gtest.h>

namespace geos
{
namespace testing
{

/**
 * @brief Simple scope-based manager for initialization/finalization in tests.
 */
class LinearAlgebraTestScope
{
public:

  LinearAlgebraTestScope( int argc, char * * argv )
  {
    ::testing::InitGoogleTest( &argc, argv );
    geos::setupEnvironment( argc, argv );
    setupLAI();
  }

  ~LinearAlgebraTestScope()
  {
    finalizeLAI();
    geos::cleanupEnvironment();
  }
};

/**
 * @name Utility functions for linear algebra unit tests.
 * @brief Functions used to construct useful matrices in the test files.
 */
///@{

/**
 * @brief Compute an identity matrix
 * @tparam MATRIX type of matrix
 * @param comm MPI communicator
 * @param N global size of the square identity matrix
 * @param I the output matrix
 *
 * This function computes the identity matrix. It can be used to generate a dummy
 * preconditioner.
 */
template< typename MATRIX >
void computeIdentity( MPI_Comm comm,
                      globalIndex const N,
                      MATRIX & I )
{
  // Construct a local CRS matrix
  int const rank = MpiWrapper::commRank( comm );
  int const nproc = MpiWrapper::commSize( comm );

  localIndex const localRowSize = LvArray::integerConversion< localIndex >( N / nproc );
  localIndex const rowResidual = LvArray::integerConversion< localIndex >( N % nproc );

  CRSMatrix< real64, globalIndex > matrix( rank == 0 ? localRowSize + rowResidual : localRowSize, N, 1 );
  CRSMatrixView< real64, globalIndex > matrixView = matrix.toView();

  globalIndex const ilower = rank * localRowSize + ( rank == 0 ? 0 : rowResidual );
  globalIndex const iupper = ilower + localRowSize + ( rank == 0 ? rowResidual : 0 );

  // Loop over rows to fill the matrix
  forAll< parallelDevicePolicy<> >( LvArray::integerConversion< localIndex >( iupper - ilower ),
                                    [=] GEOS_DEVICE ( localIndex const localRow )
  {
    // Set the values for global row i
    globalIndex const i = localRow + ilower;
    matrixView.insertNonZero( localRow, i, 1 );
  } );

  // Construct the 2D Laplace matrix
  I.create( matrix.toViewConst(), matrix.numRows(), comm );

}

/**
 * @brief Compute the 2D Laplace operator
 * @tparam MATRIX type of matrix
 * @param comm      MPI communicator.
 * @param n         size of the nxn mesh for the square 2D Laplace operator matrix. Matrix size will be N=n^2.
 * @param laplace2D the output matrix
 *
 * This function computes the matrix corresponding to a 2D Laplace operator. These
 * matrices arise from a classical finite volume formulation on a cartesian mesh
 * (5-point stencil).
 */
template< typename MATRIX >
void compute2DLaplaceOperator( MPI_Comm comm,
                               globalIndex const n,
                               MATRIX & laplace2D )
{
  // total dofs = n^2
  globalIndex N = n * n;

  // Construct a local CRSMatrix
  int const rank = MpiWrapper::commRank( comm );
  int const nproc = MpiWrapper::commSize( comm );

  localIndex const localRowSize = LvArray::integerConversion< localIndex >( N / nproc );
  localIndex const rowResidual = LvArray::integerConversion< localIndex >( N % nproc );

  CRSMatrix< real64, globalIndex > matrix( rank == 0 ? localRowSize + rowResidual : localRowSize, N, 5 );
  CRSMatrixView< real64, globalIndex > matrixView = matrix.toView();

  globalIndex const ilower = rank * localRowSize + ( rank == 0 ? 0 : rowResidual );
  globalIndex const iupper = ilower + localRowSize + ( rank == 0 ? rowResidual : 0 );

  // Loop over rows to fill the matrix
  forAll< parallelDevicePolicy<> >( LvArray::integerConversion< localIndex >( iupper - ilower ),
                                    [matrixView, n, N, ilower] GEOS_DEVICE ( localIndex const localRow )
  {

    // Allocate arrays to fill the matrix (values and columns)
    real64 values[5];
    globalIndex cols[5];

    globalIndex const globalRow = localRow + ilower;

    // Re-set the number of non-zeros for row globalRow to 0.
    localIndex nnz = 0;

    // The left -n: position globalRow-n
    if( globalRow - n >= 0 )
    {
      cols[nnz] = globalRow - n;
      values[nnz] = -1.0;
      ++nnz;
    }

    // The left -1: position globalRow-1
    if( globalRow - 1 >= 0 )
    {
      cols[nnz] = globalRow - 1;
      values[nnz] = -1.0;
      ++nnz;
    }

    // Set the diagonal: position globalRow
    cols[nnz] = globalRow;
    values[nnz] = 4.0;
    ++nnz;

    // The right +1: position globalRow+1
    if( globalRow + 1 < N )
    {
      cols[nnz] = globalRow + 1;
      values[nnz] = -1.0;
      ++nnz;
    }

    // The right +n: position globalRow+n
    if( globalRow + n < N )
    {
      cols[nnz] = globalRow + n;
      values[nnz] = -1.0;
      ++nnz;
    }

    // Set the values for row globalRow
    matrixView.insertNonZeros( localRow, cols, values, nnz );
  } );

  // Construct the 2D Laplace matrix
  laplace2D.create( matrix.toViewConst(), matrix.numRows(), comm );
}

/**
 * @brief Compute a 1st order FEM local stiffness matrix for a quad element.
 * @param hx element width
 * @param hy element height
 * @param E Young's modulus
 * @param nu Poisson ratio
 * @return the quad element stiffness matrix
 */
GEOS_HOST_DEVICE
inline stackArray2d< real64, 8*8 > Q12d_local( real64 const & hx,
                                               real64 const & hy,
                                               real64 const & E,
                                               real64 const & nu )
{
  real64 fac = E / ( 1. - 2. * nu ) / ( 1. + nu );

  // Populate stiffness matrix
  stackArray2d< real64, 8*8 > Ke( 8, 8 );

  // --- Fill diagonal entries
  real64 Dxx = ( fac * hx * ( 1. - 2. * nu ) ) / ( 6. * hy )
               - ( fac * hy * ( -1. + nu ) ) / ( 3. * hx );
  real64 Dyy = ( fac * hy * ( 1. - 2. * nu ) ) / ( 6. * hx )
               - ( fac * hx * ( -1. + nu ) ) / ( 3. * hy );
  for( localIndex i = 0; i < 8; i += 2 )
  {
    Ke[i][i] = Dxx;
    Ke[i+1][i+1] = Dyy;
  }

  // --- Fill upper triangular part
  // --- --- Ke( 0, 1:7 )
  Ke[0][1] = fac / 8.;
  Ke[0][2] = ( fac * hx * ( 1. - 2. * nu ) ) / ( 12. * hy )
             + ( fac * hy * ( -1. + nu ) ) / ( 3. * hx );
  Ke[0][3] = ( fac * ( -1 + 4. * nu ) ) / 8.;
  Ke[0][4] = -( fac * hy * ( -1. + nu ) ) / ( 6. * hx )
             + ( fac * hx * ( -1. + 2. * nu ) ) / ( 6. * hy );
  Ke[0][5] = -( fac * ( -1. + 4. * nu ) ) / 8.;
  Ke[0][6] = ( fac * hy * ( -1. + nu ) ) / ( 6. * hx )
             + ( fac * hx * ( -1. + 2. * nu ) ) / ( 12. * hy );
  Ke[0][7] = -Ke[0][1];

  // --- --- Ke( 1, 2:7 )
  Ke[1][2] = Ke[0][5];
  Ke[1][3] = -( fac * ( hy * hy * ( 1. - 2. * nu ) + hx * hx * ( -1. + nu ) ) ) / ( 6. * hx * hy );
  Ke[1][4] = Ke[0][3];
  Ke[1][5] = ( fac * hy * ( 1. - 2. * nu ) ) / ( 12. * hx )
             + ( fac * hx * ( -1. + nu ) ) / ( 3. * hy );
  Ke[1][6] = Ke[0][7];
  Ke[1][7] = ( fac * hx * ( -1. + nu ) ) / ( 6. * hy ) + ( fac * hy * ( -1. + 2. * nu ) ) / ( 12. * hx );

  // --- --- Ke( 2, 3:7 )
  Ke[2][3] = Ke[0][7];
  Ke[2][6] = Ke[0][4];
  Ke[2][7] = Ke[1][4];
  Ke[2][4] = ( fac * hy * ( -1 + nu ) ) / ( 6. * hx ) + ( fac * hx * ( -1. + 2. * nu ) ) / ( 12. * hy );
  Ke[2][5] = Ke[0][1];

  // --- --- Ke( 3, 4:7 )
  Ke[3][4] = Ke[0][1];
  Ke[3][5] = Ke[1][7];
  Ke[3][6] = Ke[1][2];
  Ke[3][7] = Ke[1][5];

  // --- --- Ke( 4, 5:7 )
  Ke[4][5] = Ke[0][7];
  Ke[4][6] = Ke[0][2];
  Ke[4][7] = Ke[0][5];

  // --- --- Ke( 5, 6:7 )
  Ke[5][6] = Ke[0][3];
  Ke[5][7] = Ke[1][3];

  // --- --- Ke( 6, 7 )
  Ke[6][7] = Ke[0][1];

  // --- Fill lower triangular part
  for( localIndex i = 1; i < 8; ++i )
  {
    for( localIndex j = 0; j < i; ++j )
    {
      Ke[i][j] = Ke[j][i];
    }
  }

  return Ke;
}

/**
 * @brief Compute cell dof indices for a quad element.
 * @param iCell global cell index
 * @param nCellsX number of cells in the X-direction
 * @param localDofIndex indices of local degrees of freedom
 */
GEOS_HOST_DEVICE
inline void computeQuadElementDofIndices( globalIndex const & iCell,
                                          globalIndex const & nCellsX,
                                          globalIndex (& localDofIndex)[8] )
{
  globalIndex cellNodes[4];
  cellNodes[0] = iCell / nCellsX + iCell;
  cellNodes[1] = cellNodes[0] + 1;
  cellNodes[2] = cellNodes[1] + nCellsX;
  cellNodes[3] = cellNodes[2] + 1;

  for( integer i = 0; i < 4; ++i )
  {
    localDofIndex[2*i] = cellNodes[i] * 2;
    localDofIndex[2*i+1] = localDofIndex[2*i] + 1;
  }
}

/**
 * @brief Compute the 2D elasticity (plane strain) operator.
 * @tparam MATRIX type of output matrix
 * @param comm MPI communicator
 * @param domainSizeX domain size in the X-direction
 * @param domainSizeY domain size in the Y-direction
 * @param nCellsX number of cells in the X-direction
 * @param nCellsY number of cells in the Y-direction
 * @param youngModulus Young's modulus value (same for all cells)
 * @param poissonRatio Poisson's ratio value (same for all cells)
 * @param elasticity2D the output matrix
 *
 * This function computes the matrix corresponding to a 2D elasticity operator,
 * assuming plane strain conditions, based on a Q1 finite element discretization.
 * The discretized domain has dimensions domainSizeX * domainSizeY.  A regular grid
 * consting of nCellsX * nCellsY cells is used. The medium is characterized by
 * homogeneous Young's modulus and Poisson's ratio. The assembled matrix is
 * singular, meaning that Dirichlet boundary conditions have not been enforced.
 */
template< typename MATRIX >
void compute2DElasticityOperator( MPI_Comm const comm,
                                  real64 const domainSizeX,
                                  real64 const domainSizeY,
                                  globalIndex const nCellsX,
                                  globalIndex const nCellsY,
                                  real64 const youngModulus,
                                  real64 const poissonRatio,
                                  MATRIX & elasticity2D )
{
  int const rank = MpiWrapper::commRank( comm );
  int const nproc = MpiWrapper::commSize( comm );

  GEOS_ERROR_IF( nCellsY < nproc, "Less than one cell row per processor is not supported" );

  real64 const hx = domainSizeX / nCellsX;
  real64 const hy = domainSizeY / nCellsY;

  // Compute cell partitioning
  globalIndex const nCells = nCellsX * nCellsY;
  localIndex const nLocalCells = LvArray::integerConversion< localIndex >( nCellsY / nproc * nCellsX );
  localIndex const nExtraCells = LvArray::integerConversion< localIndex >( nCells ) - nLocalCells * nproc;
  globalIndex const iCellLower = rank * nLocalCells + ( rank == 0 ? 0 : nExtraCells );
  globalIndex const iCellUpper = iCellLower + nLocalCells + ( rank == 0 ? nExtraCells : 0 );

  // Compute node partitioning
  globalIndex const iNodeLower = iCellLower / nCellsX + iCellLower;
  globalIndex const iNodeUpper = iCellUpper / nCellsX + iCellUpper + ( rank == nproc - 1 ? nCellsX + 1 : 0 );
  localIndex const nLocalNodes = iNodeUpper - iNodeLower;

  // Construct a local CRSMatrix
  SparsityPattern< globalIndex > sparsity( nLocalNodes * 2,
                                           ( nCellsX + 1 ) * ( nCellsY + 1 ) * 2,
                                           18 );

  // Construct local stiffness matrix (same for all cells)
  stackArray2d< real64, 8*8 > const Ke =  Q12d_local( hx, hy, youngModulus, poissonRatio );

  globalIndex const iStart = LvArray::math::max( iCellLower - nCellsX, globalIndex( 0 ) );
  globalIndex const iEnd   = LvArray::math::min( iCellUpper + nCellsX, nCells );

  globalIndex const minGlobalDof = iNodeLower*2;
  globalIndex const maxGlobalDof = iNodeUpper*2 - 1;

#if 1
  for( globalIndex iCell = iStart; iCell < iEnd; ++iCell )
  {
    // Loop over grid cells
    globalIndex dofIndex[8];

    computeQuadElementDofIndices( iCell, nCellsX, dofIndex );

    // Element matrix assembly into the local matrix

    for( integer i = 0; i < 8; ++i )
    {
      if( ( minGlobalDof <= dofIndex[i] ) && ( dofIndex[i] <= maxGlobalDof ) )
      {
        sparsity.insertNonZeros( LvArray::integerConversion< localIndex >( dofIndex[i] - minGlobalDof ),
                                 dofIndex,
                                 dofIndex + 8 );
      }
    }
  }
#else // not working on Lassen
  SparsityPatternView< globalIndex > sparsityView = sparsity.toView();
  forAll< parallelDevicePolicy<> >( iEnd-iStart, [=] GEOS_DEVICE ( localIndex const iCell )
  {
    // Loop over grid cells
    globalIndex localDofIndex[8];

    computeQuadElementDofIndices( iCell + iStart, nCellsX, localDofIndex );

    // Element matrix assembly into the local matrix

    for( integer i = 0; i < 8; ++i )
    {
      if( ( minGlobalDof <= dofIndex[i] ) && ( localDofIndex[i] <= maxGlobalDof )  )
      {
        sparsityView.insertNonZeros( LvArray::integerConversion< localIndex >( dofIndex[i] - minGlobalDof ),
                                     localDofIndex,
                                     localDofIndex + 8 );
      }
    }
  } );
#endif

  sparsity.compress();

  CRSMatrix< real64, globalIndex > matrix;
  matrix.assimilate< parallelDevicePolicy<> >( std::move( sparsity ) );
  CRSMatrixView< real64, globalIndex > matrixView = matrix.toView();

  forAll< parallelDevicePolicy<> >( LvArray::integerConversion< localIndex >( iEnd - iStart ),
                                    [=] GEOS_DEVICE ( localIndex const iCell )
  {
    // Loop over grid cells
    globalIndex dofIndex[8];
    computeQuadElementDofIndices( iCell + iStart, nCellsX, dofIndex );

    // Element matrix assembly into the local matrix
    for( integer i = 0; i < 8; ++i )
    {
      if( ( minGlobalDof <= dofIndex[i] ) && ( dofIndex[i] <= maxGlobalDof )  )
      {
        matrixView.addToRowBinarySearch< parallelDeviceAtomic >( LvArray::integerConversion< localIndex >( dofIndex[i] - minGlobalDof ),
                                                                 dofIndex,
                                                                 Ke[i],
                                                                 8 );
      }
    }
  } );

  // Impose Dirichlet boundary conditions: fix domain bottom (first 2*(nCellsX + 1) rows of matrix)
  if( rank == 0 )
  {
    forAll< parallelDevicePolicy<> >( 2 * (nCellsX + 1), [=] GEOS_DEVICE ( localIndex const localRow )
    {
      arraySlice1d< globalIndex const > const columns = matrixView.getColumns( localRow );
      arraySlice1d< real64 > const entries = matrixView.getEntries( localRow );
      localIndex const numEntries = matrixView.numNonZeros( localRow );

      for( localIndex j = 0; j < numEntries; ++j )
      {
        if( columns[ j ] != localRow )
        {
          entries[ j ] = 0;
        }
      }
    } );
  }

  // Construct the 2D elasticity (plane strain) operator
  elasticity2D.create( matrix.toViewConst(), matrix.numRows(), comm );
}

///@}

} // namespace testing
} // namespace geos

#endif //GEOS_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP_
