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
 * @file testLAIHelperFunctions.cpp
 */

#include "common/DataTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mainInterface/GeosxState.hpp"
#include "unitTests/linearAlgebraTests/testDofManagerUtils.hpp"

#include <gtest/gtest.h>

using namespace geos;

char const * xmlInput =
  R"xml(
  <Problem>
    <Mesh>
      <InternalMesh name="mesh1"
                    elementTypes="{C3D8}"
                    xCoords="{0, 1}"
                    yCoords="{0, 1}"
                    zCoords="{0, 1}"
                    nx="{6}"
                    ny="{9}"
                    nz="{5}"
                    cellBlockNames="{block1}"/>
    </Mesh>
    <ElementRegions>
      <CellElementRegion name="region1" cellBlocks="{block1}" materialList="{dummy}"/>
    </ElementRegions>
  </Problem>
  )xml";

template< typename LAI >
class LAIHelperFunctionsTest : public ::testing::Test
{
protected:

  using Base = ::testing::Test;

  LAIHelperFunctionsTest():
    Base(),
    state( std::make_unique< CommandLineOptions >() )
  {
    geos::testing::setupProblemFromXML( &state.getProblemManager(), xmlInput );
    mesh = &state.getProblemManager().getDomainPartition().getMeshBody( 0 ).getBaseDiscretization();
  }

  GeosxState state;
  MeshLevel * mesh;
};

TYPED_TEST_SUITE_P( LAIHelperFunctionsTest );

void assembleGlobalIndexVector( arrayView1d< globalIndex const > const & localToGlobal,
                                arrayView1d< integer const > const & ghostRank,
                                arrayView1d< globalIndex const > const & dofNumber,
                                globalIndex const rankOffset,
                                integer const numDofPerPoint,
                                arrayView1d< real64 > const & values )
{
  forAll< geos::parallelDevicePolicy<> >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    if( dofNumber[k] >= 0 && ghostRank[k] < 0 )
    {
      for( localIndex c = 0; c < numDofPerPoint; ++c )
      {
        values[dofNumber[k] - rankOffset + c] = static_cast< real64 >( localToGlobal[k] * numDofPerPoint + c );
      }
    }
  } );
}

TYPED_TEST_P( LAIHelperFunctionsTest, nodalVectorPermutation )
{
  using Matrix = typename TypeParam::ParallelMatrix;
  using Vector = typename TypeParam::ParallelVector;

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getBaseDiscretization();
  NodeManager const & nodeManager = meshLevel.getNodeManager();

  string const fieldName = "nodalVariable";
  integer constexpr numDofPerNode = 3;

  DofManager dofManager( "test" );
  dofManager.setDomain( domain );

  std::vector< DofManager::FieldSupport > regions;
  DofManager::FieldSupport region = { "mesh1", "Level0", {"region1"} };
  regions.emplace_back( region );

  dofManager.addField( "nodalVariable", FieldLocation::Node, 3, regions );
  dofManager.addCoupling( "nodalVariable", "nodalVariable", DofManager::Connector::Elem );
  dofManager.reorderByRank();

  Vector nodalVariable;
  nodalVariable.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );
  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< real64 > const nodalVariableView = nodalVariable.open();
  assembleGlobalIndexVector( nodeManager.localToGlobalMap(),
                             nodeManager.ghostRank(),
                             nodeManager.getReference< globalIndex_array >( dofManager.getKey( fieldName ) ),
                             rankOffset,
                             numDofPerNode,
                             nodalVariableView );
  nodalVariable.close();

  Matrix permutationMatrix;
  LAIHelperFunctions::createPermutationMatrix( nodeManager,
                                               numDofPerNode,
                                               dofManager.getKey( fieldName ),
                                               permutationMatrix );
  Vector permutedVector = LAIHelperFunctions::permuteVector( nodalVariable, permutationMatrix );

  // After permutation, the vector should become an increasing sequence of numbers
  arrayView1d< real64 const > permutedValues = permutedVector.values();
  forAll< serialPolicy >( permutedValues.size(), [=]( localIndex const i )
  {
    EXPECT_EQ( permutedValues[i], rankOffset + i );
  } );
}

TYPED_TEST_P( LAIHelperFunctionsTest, cellCenteredVectorPermutation )
{
  using Matrix = typename TypeParam::ParallelMatrix;
  using Vector = typename TypeParam::ParallelVector;

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getBaseDiscretization();
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  string const fieldName = "cellCenteredVariable";
  integer constexpr numDofPerCell = 3;

  DofManager dofManager( "test" );
  dofManager.setDomain( domain );

  std::vector< DofManager::FieldSupport > regions;
  DofManager::FieldSupport region = { "mesh1", "Level0", {"region1"} };
  regions.emplace_back( region );

  dofManager.addField( fieldName, FieldLocation::Elem, numDofPerCell, regions );
  dofManager.addCoupling( fieldName, fieldName, DofManager::Connector::Face );
  dofManager.reorderByRank();

  Vector cellCenteredVariable;
  cellCenteredVariable.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );
  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< real64 > const cellCenteredVariableView = cellCenteredVariable.open();
  elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & elementSubRegion )
  {
    assembleGlobalIndexVector( elementSubRegion.localToGlobalMap(),
                               elementSubRegion.ghostRank(),
                               elementSubRegion.getReference< globalIndex_array >( dofManager.getKey( fieldName ) ),
                               rankOffset,
                               numDofPerCell,
                               cellCenteredVariableView );
  } );
  cellCenteredVariable.close();

  Matrix permutationMatrix;
  LAIHelperFunctions::createPermutationMatrix( elemManager,
                                               numDofPerCell,
                                               dofManager.getKey( fieldName ),
                                               permutationMatrix );
  Vector permutedVector = LAIHelperFunctions::permuteVector( cellCenteredVariable, permutationMatrix );

  // After permutation, the vector should become an increasing sequence of numbers
  arrayView1d< real64 const > permutedValues = permutedVector.values();
  forAll< serialPolicy >( permutedValues.size(), [=]( localIndex const i )
  {
    EXPECT_EQ( permutedValues[i], rankOffset + i );
  } );
}

REGISTER_TYPED_TEST_SUITE_P( LAIHelperFunctionsTest,
                             nodalVectorPermutation,
                             cellCenteredVectorPermutation );

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, LAIHelperFunctionsTest, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, LAIHelperFunctionsTest, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, LAIHelperFunctionsTest, PetscInterface, );
#endif

/**
 * @function main
 * @brief Main function to setup the GEOSX environment, read the xml file and run all cases.
 */
int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
