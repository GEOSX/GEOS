/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "common/MpiWrapper.hpp"
#include "common/initializeEnvironment.hpp"

#include <gtest/gtest.h>

using namespace geos;


class MpiTestScope
{
public:

  MpiTestScope( int argc, char * argv[] )
  {
    ::testing::InitGoogleTest( &argc, argv );
    geos::setupEnvironment( argc, argv );
  }

  ~MpiTestScope()
  {
    geos::cleanupEnvironment();
  }

};

template< typename FIRST, typename SECOND >
struct PairTestCase
{
  string_view testName;
  std::vector< std::vector< MpiWrapper::PairType< FIRST, SECOND > > > pairs; // Pair set for each rank
  MpiWrapper::PairType< FIRST, SECOND > expectedMin; // Expected result for Min reduction
  MpiWrapper::PairType< FIRST, SECOND > expectedMax; // Expected result for Max reduction
};

template< typename FIRST, typename SECOND >
string reductionTestFailureMsg( PairTestCase< FIRST, SECOND > const & testCase,
                                MpiWrapper::PairType< FIRST, SECOND > pair,
                                int rankId,
                                string_view opName )
{
  return GEOS_FMT( "Test case '{}', rank {}:\n  Error in {} reduction, incorrect first value in pair ({}, {}}).",
                   testCase.testName, rankId, opName, pair.first, pair.second );
}

/**
 * @name Testing the reduction of single pair per ranks
 */
namespace pairTestCases
{

template< typename FIRST, typename SECOND >
void runAllTestCases( int const rankId )
{
  std::vector< PairTestCase< FIRST, SECOND > > const testCases = {
    {
      "Basic case with distinct values",
      { { { 9, 23 }, { 18, 886 }, { 19, 26 }, { 8, 549 } } }, // pairs
      { 8, 549 }, // expectedMin
      { 19, 26 } // expectedMax
    },
    {
      "Identical second values",
      { { { 23, 12 }, { 886, 12 }, { -26, 12 }, { 549, 12 } } }, // pairs
      { -26, 12 }, // expectedMin
      { 886, 12 } // expectedMax
    },
    {
      "Identical first values",
      { { { 12, 23 }, { 12, 886 }, { 12, -26 }, { 12, 549 } } }, // pairs
      { 12, -26 }, // expectedMin
      { 12, 886 } // expectedMax
    },
  };
  for( auto const & testCase : testCases )
  {
    MpiWrapper::PairType< FIRST, SECOND > minPair = MpiWrapper::minPair( testCase.pairs[0][rankId] );
    EXPECT_EQ( minPair.first, testCase.expectedMin.first ) << reductionTestFailureMsg( testCase, minPair, rankId, "Min" );
    EXPECT_EQ( minPair.second, testCase.expectedMin.second ) << reductionTestFailureMsg( testCase, minPair, rankId, "Min" );
    MpiWrapper::PairType< FIRST, SECOND > maxPair = MpiWrapper::maxPair( testCase.pairs[0][rankId] );
    EXPECT_EQ( maxPair.first, testCase.expectedMax.first ) << reductionTestFailureMsg( testCase, maxPair, rankId, "Max" );
    EXPECT_EQ( maxPair.second, testCase.expectedMax.second ) << reductionTestFailureMsg( testCase, maxPair, rankId, "Max" );
  }
}

TEST( MpiWrapperTesting, MpiPairReductionTest )
{
  int const rankId = MpiWrapper::commRank();
  int const nbRanks = MpiWrapper::commSize();
  if( nbRanks==4 )
  { // here, we only want to test the MIN/MAXLOC behaviour
    // basic MPI types
    runAllTestCases< int32_t, int32_t >( rankId );
    runAllTestCases< double, int32_t >( rankId );

    // custom MPI types
    runAllTestCases< int64_t, int64_t >( rankId );
    runAllTestCases< double, int64_t >( rankId );
    runAllTestCases< double, double >( rankId );
  }
}

} /* namespace pairTestCases */

/**
 * @name Testing the reduction of sets over ranks
 */
namespace pairSetsTestCases
{

template< typename FIRST, typename SECOND >
void runTestCase( PairTestCase< FIRST, SECOND > const & testCase, int rankId, int nbRanks )
{
  // Create local pairs for 1 or n ranks
  std::vector< MpiWrapper::PairType< FIRST, SECOND > > localPairs;
  if( nbRanks == 1 )
  {
    // Cas spécial : tous les rangs sont combinés en un seul
    for( auto const & rankPairs : testCase.pairs )
    {
      localPairs.insert( localPairs.end(), rankPairs.begin(), rankPairs.end());
    }
  }
  else
  {
    // Cas standard : chaque rang récupère ses propres données
    ASSERT_EQ( testCase.pairs.size(), nbRanks )
      << "Test case " << testCase.testName << ": The number of ranks does not match the test case data.";
    localPairs = testCase.pairs[rankId];
  }

  MpiWrapper::PairType< FIRST, SECOND > minPair = MpiWrapper::minPair< FIRST, SECOND >( localPairs );
  EXPECT_EQ( minPair.first, testCase.expectedMin.first ) << reductionTestFailureMsg( testCase, minPair, rankId, "Min" );
  EXPECT_EQ( minPair.second, testCase.expectedMin.second ) << reductionTestFailureMsg( testCase, minPair, rankId, "Min" );

  MpiWrapper::PairType< FIRST, SECOND > maxPair = MpiWrapper::maxPair< FIRST, SECOND >( localPairs );
  EXPECT_EQ( maxPair.first, testCase.expectedMax.first ) << reductionTestFailureMsg( testCase, maxPair, rankId, "Max" );
  EXPECT_EQ( maxPair.second, testCase.expectedMax.second ) << reductionTestFailureMsg( testCase, maxPair, rankId, "Max" );
}

template< typename FIRST, typename SECOND >
void runAllTestCases( int rankId, int nbRanks )
{
  // the set of cases is designed for 4 ranks
  std::vector< PairTestCase< FIRST, SECOND > > const testCases = {
    {
      "Basic case with distinct values",
      { { { 55, 549 }, { 34, 886 } }, { { 68, 23 }, { 10, 6 } }, { { 59, -781 }, { 64, 823 } }, { { 15, 543 }, { 22, 363 } } }, // pairs
      { 10, 6 }, // expectedMin
      { 68, 23 } // expectedMax
    },
    {
      "Distinct values but the first two ranks has no data",
      { { }, { }, { { 59, -781 }, { 64, 823 } }, { { 15, 543 }, { 22, 363 } } }, // pairs
      { 15, 543 }, // expectedMin
      { 64, 823 } // expectedMax
    },
    {
      "Distinct values, but with different data count and the last two ranks has no data",
      { { { 55, 549 }, { 34, 886 }, { 42, 23 }, { 10, 6 }, { 59, -781 }, { 64, 823 } }, { { 68, 23 }, { 22, 363 } }, { }, { } }, // pairs
      { 10, 6 }, // expectedMin
      { 68, 23 } // expectedMax
    },
    {
      "No rank have any data",
      { { }, { }, { }, { } }, // pairs
      { std::numeric_limits< FIRST >::max(), std::numeric_limits< SECOND >::max() }, // expectedMin
      { std::numeric_limits< FIRST >::lowest(), std::numeric_limits< SECOND >::lowest() }, // expectedMax
    },
    {
      "All pairs have the same value but different indices",
      { { { 12, 549 }, { 12, 886 } }, { { 12, 23 }, { 12, -6 } }, { { 12, 781 }, { 12, 823 } }, { { 12, 543 }, { 12, 363 } } }, // pairs
      { 12, -6 }, // expectedMin
      { 12, 886 } // expectedMax
    },
    {
      "All pairs have the same value but different indices, variation",
      { { { 12, -6 }, { 12, 886 } }, { { 12, 23 }, { 12, 549 } }, { { 12, 781 }, { 12, 823 } }, { { 12, 543 }, { 12, 363 } } }, // pairs
      { 12, -6 }, // expectedMin
      { 12, 886 } // expectedMax
    },
    {
      "All pairs have the same value except one",
      { { { 12, 549 }, { 12, 886 } }, { { 3, 23 }, { 12, 6 } }, { { 12, 781 }, { 12, 823 } }, { { 12, 543 }, { 12, 363 } } }, // pairs
      { 3, 23 }, // expectedMin
      { 12, 886 } // expectedMax
    },
    {
      "A unique value for one rank",
      { { { 50, 50 }, { 50, 50 } }, { { 75, -23 }, { 50, 50 } }, { { 50, 50 }, { 50, 50 } }, { { 50, 50 }, { 50, 50 } } }, // pairs
      { 50, 50 }, // expectedMin
      { 75, -23 } // expectedMax
    },
    {
      "All pairs are identical",
      { { { 123, 456 }, { 123, 456 } }, { { 123, 456 }, { 123, 456 } }, { { 123, 456 }, { 123, 456 } }, { { 123, 456 }, { 123, 456 } } }, // pairs
      { 123, 456 }, // expectedMin
      { 123, 456 } // expectedMax
    },
  };
  for( auto const & testCase : testCases )
  {
    runTestCase< FIRST, SECOND >( testCase, rankId, nbRanks );
  }
}

TEST( MpiWrapperTesting, MpiPairSetReductionTest )
{
  int const rankId = MpiWrapper::commRank();
  int const nbRanks = MpiWrapper::commSize();

  // basic MPI types
  runAllTestCases< int32_t, int32_t >( rankId, nbRanks );
  runAllTestCases< double, int32_t >( rankId, nbRanks );

  // custom MPI types
  runAllTestCases< int64_t, int64_t >( rankId, nbRanks );
  runAllTestCases< double, int64_t >( rankId, nbRanks );
  runAllTestCases< double, double >( rankId, nbRanks );
}

} /* namespace pairSetsTestCases */

int main( int argc, char * argv[] )
{
  MpiTestScope testScope{ argc, argv };
  return RUN_ALL_TESTS();
}
