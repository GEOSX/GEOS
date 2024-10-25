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
 * @file testDofManager.cpp
 * @brief This test file is part of the ctest suite and tests the DofManager functionality.
 */

#include "codingUtilities/UnitTestUtilities.hpp"
#include "common/DataTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mainInterface/GeosxState.hpp"
#include "unitTests/linearAlgebraTests/testDofManagerUtils.hpp"

#include "gtest/gtest.h"

#include <memory>

using namespace geos;
using namespace geos::testing;
using namespace geos::dataRepository;

char const * xmlInput =
  R"xml(
  <Problem>
    <Mesh>
      <InternalMesh name="mesh"
                    elementTypes="{C3D8}"
                    xCoords="{0, 1, 2, 3, 4}"
                    yCoords="{0, 1}"
                    zCoords="{0, 1}"
                    nx="{4, 4, 4, 4}"
                    ny="{4}"
                    nz="{5}"
                    cellBlockNames="{block1, block2, block3, block4}"/>
    </Mesh>
    <ElementRegions>
      <CellElementRegion name="region1" cellBlocks="{block1}" materialList="{}" />
      <CellElementRegion name="region2" cellBlocks="{block2}" materialList="{}" />
      <CellElementRegion name="region3" cellBlocks="{block3}" materialList="{}" />
      <CellElementRegion name="region4" cellBlocks="{block4}" materialList="{}" />
    </ElementRegions>
  </Problem>
  )xml";

/**
 * @brief Base class for all DofManager test fixtures.
 */
class DofManagerTestBase : public ::testing::Test
{
protected:

  using Base = ::testing::Test;

  DofManagerTestBase():
    Base(),
    state( std::make_unique< CommandLineOptions >() ),
    domain( state.getProblemManager().getDomainPartition() ),
    dofManager( "test" )
  {
    geos::testing::setupProblemFromXML( &state.getProblemManager(), xmlInput );
    dofManager.setDomain( domain );
  }

  GeosxState state;
  DomainPartition & domain;
  DofManager dofManager;
};

/**
 * @brief Check that all local objects have dofs assigned. Collect all local dofs into a list.
 * @param [in] domain the domain partition
 * @param [in] dofIndexKey key used to look up dof index array
 * @param [in] support list of target region names
 * @param [out] dofNumbers the list of dof indices
 */
template< FieldLocation LOC >
void collectLocalDofNumbers( DomainPartition const & domain,
                             string const & dofIndexKey,
                             std::vector< DofManager::FieldSupport > const & support,
                             array1d< globalIndex > & dofNumbers )
{
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & meshLevel = meshBody.getMeshLevel( regions.meshLevelName );

    ObjectManagerBase const & manager = meshLevel.getGroup< ObjectManagerBase >
                                          ( geos::testing::internal::testMeshHelper< LOC >::managerKey() );
    arrayView1d< globalIndex const > dofIndex = manager.getReference< array1d< globalIndex > >( dofIndexKey );

    forLocalObjects< LOC >( meshLevel, regions.regionNames, [&]( localIndex const idx )
    {
      SCOPED_TRACE( "idx = " + std::to_string( idx ) );
      EXPECT_GE( dofIndex[idx], 0 );
      dofNumbers.emplace_back( dofIndex[idx] );
    } );
  }
}

/**
 * @brief Collect all local dofs into a list.
 * @param [in] domain the domain partition
 * @param [in] dofIndexKey key used to look up dof index array
 * @param [in] support list of target region names
 * @param [out] dofNumbers the list of dof indices
 */
template<>
void collectLocalDofNumbers< FieldLocation::Elem >( DomainPartition const & domain,
                                                    string const & dofIndexKey,
                                                    std::vector< DofManager::FieldSupport > const & support,
                                                    array1d< globalIndex > & dofNumbers )
{
  for( DofManager::FieldSupport const & regions : support )
  {
    MeshBody const & meshBody = domain.getMeshBody( regions.meshBodyName );
    MeshLevel const & meshLevel = meshBody.getMeshLevel( regions.meshLevelName );

    // make a list of regions
    ElementRegionManager const & elemManager = meshLevel.getElemManager();
    auto const dofNumber = elemManager.constructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofIndexKey );

    forLocalObjects< FieldLocation::Elem >( meshLevel, regions.regionNames, [&]( auto const idx )
    {
      globalIndex const dofIndex = dofNumber[idx[0]][idx[1]][idx[2]];
      EXPECT_GE( dofIndex, 0 );
      dofNumbers.emplace_back( dofIndex );
    } );
  }
}

/**
 * @brief Checks that dof numbers are unique
 * @param [in] dofNumbers list of dof numbers, must be sorted
 */
void checkUniqueness( arrayView1d< globalIndex const > const & dofNumbers )
{
  EXPECT_TRUE( std::adjacent_find( dofNumbers.begin(), dofNumbers.end() ) == dofNumbers.end() );
}

/**
 * @brief Checks that dof numbers have a given stride between them.
 * @param [in] dofNumbers list of dof numbers, must be sorted
 * @param [in] stride expected stride between dof numbers
 */
void checkStride( arrayView1d< globalIndex const > const & dofNumbers, globalIndex const stride )
{
  // Check dof stride is equal to numComp
  for( localIndex i = 1; i < dofNumbers.size(); ++i )
  {
    SCOPED_TRACE( "i = " + std::to_string( i ) );
    EXPECT_EQ( dofNumbers[i] - dofNumbers[i - 1], stride );
  }
}

/**
 * @brief Checks that dof numbers are globally consistent across ranks.
 * @param [in] dofNumbers list of dof numbers, must be sorted
 */
void checkGlobalOrdering( arrayView1d< globalIndex const > const & dofNumbers )
{
  array1d< globalIndex > localDofBounds( 2 );
  localDofBounds[0] = dofNumbers.empty() ? -1 : dofNumbers.front();
  localDofBounds[1] = dofNumbers.empty() ? -1 : dofNumbers.back();

  array1d< globalIndex > globalDofBounds;
  MpiWrapper::allGather( localDofBounds.toViewConst(), globalDofBounds );

  if( MpiWrapper::commRank() == 0 )
  {
    localIndex const newSize =
      std::distance( globalDofBounds.begin(), std::remove( globalDofBounds.begin(), globalDofBounds.end(), -1 ));
    globalDofBounds.resize( newSize );

    // Check that lower ranks have strictly lower dof numbers
    EXPECT_TRUE( std::is_sorted( globalDofBounds.begin(), globalDofBounds.end() ) );
    checkUniqueness( globalDofBounds );
  }
}

/**
 * @brief Test fixture for all non-typed (LAI independent) DofManager tests.
 */
class DofManagerIndicesTest : public DofManagerTestBase
{
protected:

  struct FieldDesc
  {
    string name;
    FieldLocation location;
    integer components;
    std::vector< DofManager::FieldSupport > regions{};
  };

  void test( std::vector< FieldDesc > const & fields );
};

/**
 * @brief General test function for single-field dof indices.
 * @param location location of the field
 * @param fieldName name of the field
 * @param numComp number of components
 * @param numDofIndicesExpected expected global number of dof indices
 * @param regions list of support regions (empty = whole domain)
 */
void DofManagerIndicesTest::test( std::vector< FieldDesc > const & fields )
{
  for( FieldDesc const & f : fields )
  {
    dofManager.addField( f.name, f.location, f.components, f.regions );
  }
  dofManager.reorderByRank();

  array1d< globalIndex > allDofNumbers;
  localIndex lastNumComp = -1;

  for( FieldDesc const & f : fields )
  {
    array1d< globalIndex > dofNumbers;
    string const key = dofManager.getKey( f.name );
    switch( f.location )
    {
      case FieldLocation::Elem:
      {
        collectLocalDofNumbers< FieldLocation::Elem >( domain, key, getRegions( domain, f.regions ), dofNumbers );
        break;
      }
      case FieldLocation::Face:
      {
        collectLocalDofNumbers< FieldLocation::Face >( domain, key, getRegions( domain, f.regions ), dofNumbers );
        break;
      }
      case FieldLocation::Node:
      {
        collectLocalDofNumbers< FieldLocation::Node >( domain, key, getRegions( domain, f.regions ), dofNumbers );
        break;
      }
      default:
      {
        GEOS_ERROR( "Unsupported" );
      }
    }
    std::sort( dofNumbers.begin(), dofNumbers.end() );

    if( !dofNumbers.empty() )
    {
      checkUniqueness( dofNumbers );
      checkStride( dofNumbers, f.components );
      if( lastNumComp >= 0 )
      {
        EXPECT_EQ( dofNumbers.front(), allDofNumbers.back() + lastNumComp );
      }
      allDofNumbers.insert( allDofNumbers.size(), dofNumbers.begin(), dofNumbers.end() );
    }
    lastNumComp = f.components;
  }
  checkGlobalOrdering( allDofNumbers );
}

/**
 * @brief Check dof indexing for a node-based field with 3 dof/node.
 */
TEST_F( DofManagerIndicesTest, Node_Full )
{
  test(
  {
    {
      "displacement",
      FieldLocation::Node,
      3,
      { { "mesh", "Level0", {} } }
    }
  } );
}

/**
 * @brief Check dof indexing for a node-based field with 3 dof/node with partial domain support.
 */
TEST_F( DofManagerIndicesTest, Node_Partial )
{
  test(
  {
    {
      "displacement",
      FieldLocation::Node,
      3,
      { { "mesh", "Level0", { "region1", "region3", "region4" } } }
    }
  } );
}

/**
 * @brief Check dof indexing for a cell-based field with 2 dof/cell.
 */
TEST_F( DofManagerIndicesTest, Elem_Full )
{
  test(
  {
    {
      "pressure",
      FieldLocation::Elem,
      2,
      { { "mesh", "Level0", {} } }
    }
  } );
}

/**
 * @brief Check dof indexing for a cell-based field with 2 dof/cell with partial domain support.
 */
TEST_F( DofManagerIndicesTest, Elem_Partial )
{
  test(
  {
    {
      "pressure",
      FieldLocation::Elem,
      2,
      { { "mesh", "Level0", { "region1", "region3", "region4" } } }
    }
  } );
}

/**
 * @brief Check dof indexing for a typical FVM system with 2 dof/cell.
 */
TEST_F( DofManagerIndicesTest, Face_Full )
{
  test(
  {
    {
      "flux",
      FieldLocation::Face,
      2,
      { { "mesh", "Level0", {} } }
    }
  } );
}

/**
 * @brief Check dof indexing for a typical FVM system with 2 dof/cell.
 */
TEST_F( DofManagerIndicesTest, Face_Partial )
{
  test(
  {
    {
      "flux",
      FieldLocation::Face,
      2,
      { { "mesh", "Level0", { "region1", "region3", "region4" } } }
    }
  } );
}

/**
 * @brief Check dof indexing for a typical FEM/FVM system.
 */
TEST_F( DofManagerIndicesTest, Node_Elem_Full )
{
  test(
  {
    {
      "displacement",
      FieldLocation::Node,
      3,
      { { "mesh", "Level0", {} } }
    },
    {
      "pressure",
      FieldLocation::Elem,
      2,
      { { "mesh", "Level0", {} } }
    }
  } );
}

/**
 * @brief Check dof indexing for a typical FEM/FVM system with partial domain support.
 */
TEST_F( DofManagerIndicesTest, Node_Elem_Partial )
{
  test(
  {
    {
      "displacement",
      FieldLocation::Node,
      3,
      { { "mesh", "Level0", { "region1", "region3", "region4" } } }
    },
    {
      "pressure",
      FieldLocation::Elem,
      2,
      { { "mesh", "Level0", { "region1", "region2", "region4" } } }
    }
  } );
}

/**
 * @brief Check dof indexing for a typical Hybrid FVM system.
 */
TEST_F( DofManagerIndicesTest, Face_Elem_Full )
{
  test(
  {
    {
      "flux",
      FieldLocation::Face,
      2,
      { { "mesh", "Level0", {} } }
    },
    {
      "pressure",
      FieldLocation::Elem,
      2,
      { { "mesh", "Level0", {} } }
    }
  } );
}

/**
 * @brief Check dof indexing for a typical Hybrid FVM system with partial domain support.
 */
TEST_F( DofManagerIndicesTest, Face_Elem_Partial )
{
  test(
  {
    {
      "flux",
      FieldLocation::Face,
      2,
      { { "mesh", "Level0", { "region1", "region3", "region4" } } }
    },
    {
      "pressure",
      FieldLocation::Elem,
      2,
      { { "mesh", "Level0", { "region1", "region3", "region4" } } }
    }
  } );
}

/**
 * @brief Test fixture for all typed (LAI dependent) DofManager tests.
 * @tparam LAI linear algebra interface type
 */
template< typename LAI >
class DofManagerMatrixTest : public DofManagerTestBase
{
protected:

  using Matrix = typename LAI::ParallelMatrix;

  using PatternFunc = void ( * )( DomainPartition const & mesh,
                                  string const & dofIndexKey,
                                  std::vector< DofManager::FieldSupport > const & regions,
                                  globalIndex const rankOffset,
                                  localIndex const numComp,
                                  CRSMatrix< real64 > & sparsity );

  using CoupledPatternFunc = void ( * )( DomainPartition const & mesh,
                                         string const & dofIndexKey1,
                                         string const & dofIndexKey2,
                                         std::vector< DofManager::FieldSupport > const & regions,
                                         globalIndex const rankOffset,
                                         localIndex const numComp1,
                                         localIndex const numComp2,
                                         CRSMatrix< real64 > & sparsity );

  struct FieldDesc
  {
    string name;
    FieldLocation location;
    DofManager::Connector connectivity;
    localIndex components;
    PatternFunc makePattern;
    std::vector< DofManager::FieldSupport > regions = {};
  };

  struct CouplingDesc
  {
    DofManager::Connector connectivity;
    CoupledPatternFunc makeCouplingPattern;
    bool symmetric = true;
    std::vector< DofManager::FieldSupport > regions = {};
  };

  void addFields( std::vector< FieldDesc > fields,
                  std::map< std::pair< string, string >, CouplingDesc > couplings = {} )
  {
    for( FieldDesc const & f : fields )
    {
      std::vector< DofManager::FieldSupport > const regions = getRegions( domain, f.regions );
      dofManager.addField( f.name, f.location, f.components, regions );
      dofManager.addCoupling( f.name, f.name, f.connectivity );
    }
    for( auto const & entry : couplings )
    {
      std::pair< string, string > const & fieldNames = entry.first;
      CouplingDesc const & c = entry.second;
      std::vector< DofManager::FieldSupport > const regions = getRegions( domain, c.regions );
      dofManager.addCoupling( fieldNames.first, fieldNames.second, c.connectivity, regions, c.symmetric );
    }
    dofManager.reorderByRank();
  }
};

template< typename LAI >
class DofManagerSparsityTest : public DofManagerMatrixTest< LAI >
{
protected:

  using Base = DofManagerMatrixTest< LAI >;
  using Matrix = typename Base::Matrix;
  using FieldDesc = typename Base::FieldDesc;
  using CouplingDesc = typename Base::CouplingDesc;

  using Base::domain;
  using Base::dofManager;
  using Base::addFields;

  void test( std::vector< FieldDesc > fields,
             std::map< std::pair< string, string >, CouplingDesc > couplings = {} );
};

TYPED_TEST_SUITE_P( DofManagerSparsityTest );

template< typename LAI >
void DofManagerSparsityTest< LAI >::test( std::vector< FieldDesc > fields,
                                          std::map< std::pair< string, string >, CouplingDesc > couplings )
{
  addFields( fields, couplings );

  // Create a sparsity pattern via regular face loop
  localIndex numLocalDof = 0;
  localIndex numCompTotal = 0;
  for( FieldDesc const & f : fields )
  {
    std::vector< DofManager::FieldSupport > const regions = getRegions( domain, f.regions );
    localIndex numLocalObj = 0;
    switch( f.location )
    {
      case FieldLocation::Elem:
      {
        numLocalObj = countLocalObjects< FieldLocation::Elem >( domain, regions );
      }
      break;
      case FieldLocation::Face:
      {
        numLocalObj = countLocalObjects< FieldLocation::Face >( domain, regions );
      }
      break;
      case FieldLocation::Node:
      {
        numLocalObj = countLocalObjects< FieldLocation::Node >( domain, regions );
      }
      break;
      default:
        GEOS_ERROR( "Unsupported" );
    }
    numLocalDof += numLocalObj * f.components;
    numCompTotal += f.components;
  }

  Matrix pattern;
  {
    // Create a sparsity pattern via DofManager
    SparsityPattern< globalIndex > localPattern;
    dofManager.setSparsityPattern( localPattern );
    CRSMatrix< real64, globalIndex > localMatrix;
    localMatrix.assimilate< parallelHostPolicy >( std::move( localPattern ) );
    pattern.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
    pattern.set( 1.0 );
  }

  CRSMatrix< real64, globalIndex > localPatternExpected( numLocalDof,
                                                         dofManager.numGlobalDofs(),
                                                         27 * numCompTotal );
  Matrix patternExpected;

  for( FieldDesc const & f : fields )
  {
    GEOS_LOG_RANK( "rankOffset = "<<dofManager.rankOffset() );
    f.makePattern( domain,
                   dofManager.getKey( f.name ),
                   getRegions( domain, f.regions ),
                   dofManager.rankOffset(),
                   f.components,
                   localPatternExpected );
  }
  for( auto const & entry : couplings )
  {
    std::pair< string, string > const & names = entry.first;
    CouplingDesc const & c = entry.second;

    FieldDesc const & f1 = *std::find_if( fields.begin(), fields.end(),
                                          [&]( auto const & f ){ return f.name == names.first; } );
    FieldDesc const & f2 = *std::find_if( fields.begin(), fields.end(),
                                          [&]( auto const & f ){ return f.name == names.second; } );

    c.makeCouplingPattern( domain,
                           dofManager.getKey( f1.name ),
                           dofManager.getKey( f2.name ),
                           getRegions( domain, c.regions ),
                           dofManager.rankOffset(),
                           f1.components,
                           f2.components,
                           localPatternExpected );
  }
  patternExpected.create( localPatternExpected.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );

  // Compare the sparsity patterns
  pattern.set( 1.0 );
  patternExpected.set( 1.0 );
  compareMatrices( pattern, patternExpected, 0.0, 0.0 );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, TPFA_Full )
{
  TestFixture::test( {
    { "pressure",
      FieldLocation::Elem,
      DofManager::Connector::Face,
      2, makeSparsityTPFA }
  } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, TPFA_Partial )
{
  TestFixture::test( {
    { "pressure",
      FieldLocation::Elem,
      DofManager::Connector::Face,
      2, makeSparsityTPFA,
      { {"mesh", "Level0", { "region1", "region3", "region4" } } }
    }
  } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, FEM_Full )
{
  TestFixture::test( {
    { "displacement",
      FieldLocation::Node,
      DofManager::Connector::Elem,
      3, makeSparsityFEM }
  } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, FEM_Partial )
{
  TestFixture::test( {
    { "displacement",
      FieldLocation::Node,
      DofManager::Connector::Elem,
      3, makeSparsityFEM,
      { {"mesh", "Level0", { "region1", "region3", "region4" } } }
    }
  } );
}

/**
 * @brief Compare mass matrix sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, Mass_Full )
{
  TestFixture::test( {
    { "mass",
      FieldLocation::Elem,
      DofManager::Connector::None,
      2, makeSparsityMass }
  } );
}

/**
 * @brief Compare mass matrix sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, Mass_Partial )
{
  TestFixture::test( {
    { "mass",
      FieldLocation::Elem,
      DofManager::Connector::None,
      2, makeSparsityMass,
      { {"mesh", "Level0", { "region1", "region3", "region4" } } }
    }
  } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, Flux_Full )
{
  TestFixture::test( {
    { "flux",
      FieldLocation::Face,
      DofManager::Connector::Elem,
      2, makeSparsityFlux }
  } );
}

/**
 * @brief Compare TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, Flux_Partial )
{
  TestFixture::test( {
    { "flux",
      FieldLocation::Face,
      DofManager::Connector::Elem,
      2, makeSparsityFlux,
      { {"mesh", "Level0", { "region1", "region3", "region4" } } }
    }
  } );
}

/**
 * @brief Compare a mixed FEM/TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop.
 */
TYPED_TEST_P( DofManagerSparsityTest, FEM_TPFA_Full )
{
  TestFixture::test( {
    { "displacement",
      FieldLocation::Node,
      DofManager::Connector::Elem,
      3, makeSparsityFEM },
    { "pressure",
      FieldLocation::Elem,
      DofManager::Connector::Face,
      2, makeSparsityTPFA }
  },
  {
    {
      { "displacement", "pressure" },
      { DofManager::Connector::Elem,
        makeSparsityFEM_FVM,
        true }
    }
  } );
}

/**
 * @brief Compare a mixed FEM/TPFA sparsity pattern produced by DofManager against one
 *        created with a direct assembly loop, with partial domain support.
 */
TYPED_TEST_P( DofManagerSparsityTest, FEM_TPFA_Partial )
{
  TestFixture::test( {
    { "displacement",
      FieldLocation::Node,
      DofManager::Connector::Elem,
      3, makeSparsityFEM,
      { {"mesh", "Level0", { "region1", "region3", "region4" } } }
    },
    { "pressure",
      FieldLocation::Elem,
      DofManager::Connector::Face,
      2, makeSparsityTPFA,
      { {"mesh", "Level0", { "region1", "region2", "region4" } } }
    }
  },
  {
    {
      { "displacement", "pressure" },
      { DofManager::Connector::Elem,
        makeSparsityFEM_FVM,
        true,
        { {"mesh", "Level0", { "region4" } } }
      }
    }
  } );
}

REGISTER_TYPED_TEST_SUITE_P( DofManagerSparsityTest,
                             TPFA_Full,
                             TPFA_Partial,
                             FEM_Full,
                             FEM_Partial,
                             Mass_Full,
                             Mass_Partial,
                             Flux_Full,
                             Flux_Partial,
                             FEM_TPFA_Full,
                             FEM_TPFA_Partial );

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, DofManagerSparsityTest, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, DofManagerSparsityTest, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, DofManagerSparsityTest, PetscInterface, );
#endif

/**
 * @brief Test fixture for all typed (LAI dependent) DofManager tests.
 * @tparam LAI linear algebra interface type
 */
template< typename LAI >
class DofManagerRestrictorTest : public DofManagerMatrixTest< LAI >
{
protected:

  using Base = DofManagerMatrixTest< LAI >;
  using Matrix = typename Base::Matrix;
  using FieldDesc = typename Base::FieldDesc;
  using CouplingDesc = typename Base::CouplingDesc;

  using Base::domain;
  using Base::dofManager;
  using Base::addFields;

  void test( std::vector< FieldDesc > fields,
             std::vector< DofManager::SubComponent > selection,
             std::map< std::pair< string, string >, CouplingDesc > couplings = {} );
};

template< typename LAI >
void DofManagerRestrictorTest< LAI >::test( std::vector< FieldDesc > fields,
                                            std::vector< DofManager::SubComponent > selection,
                                            std::map< std::pair< string, string >, CouplingDesc > couplings )
{
  addFields( fields, couplings );

  // Create and fill the full matrix
  Matrix A;
  {
    SparsityPattern< globalIndex > localPattern;
    dofManager.setSparsityPattern( localPattern );
    CRSMatrix< real64, globalIndex > localMatrix;
    localMatrix.assimilate< parallelHostPolicy >( std::move( localPattern ) );
    A.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
    A.set( 1.0 );
  }

  // Create prolongation and restriction to 2 out of 3 components
  Matrix R, P;
  dofManager.makeRestrictor( selection, A.comm(), false, R );
  dofManager.makeRestrictor( selection, A.comm(), true, P );

  // Compute the sub-matrix via PtAP
  Matrix Asub_PtAP;
  A.multiplyPtAP( P, Asub_PtAP );

  // Compute the sub-matrix via RAP
  Matrix Asub_RAP;
  A.multiplyRAP( R, P, Asub_RAP );

  // Filter the selected fields
  std::vector< FieldDesc > selectedFields( selection.size() );
  for( std::size_t k = 0; k < selection.size(); ++k )
  {
    selectedFields[k] = *std::find_if( fields.begin(), fields.end(),
                                       [&]( FieldDesc const & f ) { return f.name == selection[k].fieldName; } );
    selectedFields[k].components = selection[k].mask.size();
  }

  // Filter the couplings of selected fields
  std::map< std::pair< string, string >, CouplingDesc > couplingsSelected;
  for( auto it = couplings.begin(); it != couplings.end(); ++it )
  {
    std::pair< string, string > const & fieldNames = it->first;
    bool const f1 = std::count_if( selectedFields.begin(), selectedFields.end(),
                                   [&]( FieldDesc const & f ) { return f.name == fieldNames.first; } ) > 0;
    bool const f2 = std::count_if( selectedFields.begin(), selectedFields.end(),
                                   [&]( FieldDesc const & f ) { return f.name == fieldNames.second; } ) > 0;
    if( f1 && f2 )
    {
      couplingsSelected.emplace( *it );
    }
  }

  // Now reset the DofManager and make a field with sub-components only
  dofManager.clear();
  addFields( selectedFields, couplingsSelected );

  // Compute the expected matrix
  Matrix B;
  {
    SparsityPattern< globalIndex > localPattern;
    dofManager.setSparsityPattern( localPattern );
    CRSMatrix< real64, globalIndex > localMatrix;
    localMatrix.assimilate< parallelHostPolicy >( std::move( localPattern ) );
    B.create( localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
    B.set( 1.0 );
  }

  // Check if the matrices match
  {
    SCOPED_TRACE( "PtAP" );
    compareMatrices( Asub_PtAP, B );
  }
  {
    SCOPED_TRACE( "RAP" );
    compareMatrices( Asub_RAP, B );
  }
}

TYPED_TEST_SUITE_P( DofManagerRestrictorTest );

TYPED_TEST_P( DofManagerRestrictorTest, SingleBlock )
{
  TestFixture::test(
  {
    { "pressure",
      FieldLocation::Elem,
      DofManager::Connector::Face,
      3, nullptr,
      { {"mesh", "Level0", {"region1", "region3", "region4"} } }
    }
  },
  {
    { "pressure", { 3, 1, 3 } }
  } );
}

TYPED_TEST_P( DofManagerRestrictorTest, MultiBlock_First )
{
  TestFixture::test(
  {
    { "displacement",
      FieldLocation::Node,
      DofManager::Connector::Elem,
      3, nullptr,
      { {"mesh", "Level0", { "region1", "region3", "region4" } } }
    },
    { "pressure",
      FieldLocation::Elem,
      DofManager::Connector::Face,
      2, nullptr,
      { {"mesh", "Level0", { "region1", "region2", "region4" } } }
    }
  },
  {
    { "displacement", { 3, 1, 3 } }
  },
  {
    {
      { "displacement", "pressure" },
      { DofManager::Connector::Elem,
        nullptr,
        true,
        { {"mesh", "Level0", { "region4" } } }
      }
    }
  } );
}

TYPED_TEST_P( DofManagerRestrictorTest, MultiBlock_Second )
{
  TestFixture::test(
  {
    { "displacement",
      FieldLocation::Node,
      DofManager::Connector::Elem,
      3, nullptr,
      { {"mesh", "Level0", { "region1", "region3", "region4" } } }
    },
    { "pressure",
      FieldLocation::Elem,
      DofManager::Connector::Face,
      2, nullptr,
      { {"mesh", "Level0", {"region1", "region2", "region4"} } }
    }
  },
  {
    { "pressure", { 2, 1, 2 } }
  },
  {
    {
      { "displacement", "pressure" },
      { DofManager::Connector::Elem,
        nullptr,
        true,
        { {"mesh", "Level0", { "region4" } } }
      }
    }
  }
    );
}

TYPED_TEST_P( DofManagerRestrictorTest, MultiBlock_Both )
{
  TestFixture::test(
  {
    { "displacement",
      FieldLocation::Node,
      DofManager::Connector::Elem,
      3, nullptr,
      { { "mesh", "Level0", { "region1", "region3", "region4" } } }
    },
    { "pressure",
      FieldLocation::Elem,
      DofManager::Connector::Face,
      2, nullptr,
      { { "mesh", "Level0", { "region1", "region2", "region4" } } }
    }
  },
  {
    { "displacement", { 3, 1, 3 } }, { "pressure", { 2, 1, 2 } }
  },
  {
    {
      { "displacement", "pressure" },
      { DofManager::Connector::Elem,
        nullptr,
        true,
        { {"mesh", "Level0", { "region4" } } }
      }
    }
  }
    );
}

REGISTER_TYPED_TEST_SUITE_P( DofManagerRestrictorTest,
                             SingleBlock,
                             MultiBlock_First,
                             MultiBlock_Second,
                             MultiBlock_Both );

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, DofManagerRestrictorTest, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, DofManagerRestrictorTest, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, DofManagerRestrictorTest, PetscInterface, );
#endif

TEST( DofManagerRegions, aggregateInitialization )
{
  // The DofManager::FieldSupport and DofManager::SubComponent are sometimes constructed by using aggregate initialization.
  // The danger is that this feature implicitly depends on the order of the parameters of the `struct`.
  // Any reordering of those parameters would result in a bug.
  // This test aims at protecting against any change in this implicit convention.
  // C++20's new feature "designated initializers" would probably make this test useless.

  {
    DofManager::FieldSupport regions0;
    regions0.meshBodyName = "meshBodyName";
    regions0.meshLevelName = "meshLevelName";
    regions0.regionNames = { "regionName0", "regionName1" };

    DofManager::FieldSupport const regions1{ regions0.meshBodyName, regions0.meshLevelName, regions0.regionNames };

    ASSERT_EQ( regions0.meshBodyName, regions1.meshBodyName );
    ASSERT_EQ( regions0.meshLevelName, regions1.meshLevelName );
    ASSERT_EQ( regions0.regionNames, regions1.regionNames );
  }

  {
    DofManager::SubComponent subComponent0;
    subComponent0.fieldName = "fieldName";
    subComponent0.mask = DofManager::CompMask( 7, true );

    DofManager::SubComponent subComponent1{ subComponent0.fieldName, subComponent0.mask };

    ASSERT_EQ( subComponent0.fieldName, subComponent1.fieldName );
    ASSERT_EQ( subComponent0.mask.numComp(), subComponent1.mask.numComp() );
  }
}

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}
