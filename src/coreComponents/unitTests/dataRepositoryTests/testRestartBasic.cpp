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

// Source includes
#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/GeosxState.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "dataRepository/ConduitRestart.hpp"
#include "common/MpiWrapper.hpp"
#include "utils.hpp"

// TPL includes
#include <gtest/gtest.h>

// System includes
#include <random>

namespace geos
{

CommandLineOptions g_commandLineOptions;

namespace dataRepository
{
namespace testing
{

template< typename T >
class SingleWrapperTest : public ::testing::Test
{
public:

  SingleWrapperTest():
    m_node( std::make_unique< conduit::Node >() ),
    m_group( std::make_unique< Group >( m_groupName, *m_node ) ),
    m_groupSize( rand( 0, 100 ) )
  {
    m_group->resize( m_groupSize );
    m_wrapper = &m_group->registerWrapper< T >( m_wrapperName );
    m_wrapperSizedFromParent = rand( 0, 100 );
    m_wrapper->setSizedFromParent( m_wrapperSizedFromParent );
  }

  void test()
  {
    T value;
    fill( value, 100 );

    // Set the value
    m_wrapper->reference() = value;

    // Write out the tree
    m_group->prepareToWrite();
    writeTree( m_fileName, *m_node );
    m_group->finishWriting();

    // Delete geos tree and reset the conduit tree.
    m_group = nullptr;
    m_node = std::make_unique< conduit::Node >();

    // Load in the tree
    loadTree( m_fileName, *m_node );
    m_group = std::make_unique< Group >( m_groupName, *m_node );
    m_wrapper = &m_group->registerWrapper< T >( m_wrapperName );
    m_group->loadFromConduit();

    // Compare metadata
    EXPECT_EQ( m_group->size(), m_groupSize );
    EXPECT_EQ( m_wrapper->sizedFromParent(), m_wrapperSizedFromParent );

    // Compare the values
    compare( value, m_wrapper->reference() );
  }

private:
  string const m_groupName = "root";
  string const m_wrapperName = "wrapper";
  string const m_fileName = "testRestartBasic_SingleWrapperTest";

  std::unique_ptr< conduit::Node > m_node;
  std::unique_ptr< Group > m_group;
  int m_groupSize;
  Wrapper< T > * m_wrapper;
  int m_wrapperSizedFromParent;
};

using TestTypes = ::testing::Types< int,
                                    double,
                                    R1Tensor,
                                    std::pair< int, R1Tensor >, // This should be passed to conduit via an external
                                                                // pointer but currently we're packing it.
                                    std::pair< string, double >,
                                    string,
                                    std::vector< int >,
                                    // std::vector< string > bufferOps currently can't pack this
                                    array1d< double >,
                                    array1d< string >,
                                    array1d< array1d< double > >,
                                    array1d< array1d< string > >,
                                    array2d< double >,
                                    array2d< string >,
                                    array2d< double, RAJA::PERM_JI >,
                                    array2d< string, RAJA::PERM_JI >,
                                    array3d< double >,
                                    array3d< string >,
                                    array3d< double, RAJA::PERM_KJI >,
                                    array3d< string, RAJA::PERM_IKJ >,
                                    SortedArray< int >,
                                    SortedArray< string >,
                                    map< string, int >,
                                    unordered_map< string, int >,
                                    map< long, int >,
                                    unordered_map< long, int >
                                    >;
TYPED_TEST_SUITE( SingleWrapperTest, TestTypes, );

TYPED_TEST( SingleWrapperTest, WriteAndRead )
{
  this->test();
}

} // namespace testing
} // namespace dataRepository
} // namespace geos

int main( int argc, char * argv[] )
{
  testing::InitGoogleTest( &argc, argv );

  geos::g_commandLineOptions = *geos::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
