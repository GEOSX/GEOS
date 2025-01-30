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

#include <iostream>
#include "codingUtilities/RTTypes.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"

#include <gtest/gtest.h>

// TPL includes
#include <conduit.hpp>

// Mock classes to test dynamic casting
#include "BaseClass.hpp"
#include "DerivedFinalClass.hpp"

using namespace geos;
using namespace dataRepository;

// Test for dynamicCast with pointer
TEST( DynamicCastTests, Pointer_Casting_Success )
{
  BaseClass * base = new DerivedFinalClass();
  DerivedFinalClass * derived = geos::dynamicCast< DerivedFinalClass * >( base );
  ASSERT_NE( derived, nullptr ) << "Expected successful cast from Base to Derived.";
  delete base;   // Clean up allocated memory
}

TEST( DynamicCastTests, Pointer_Casting_Failure )
{
  BaseClass * base = new BaseClass();
  DerivedFinalClass * derived = geos::dynamicCast< DerivedFinalClass * >( base );
  ASSERT_EQ( derived, nullptr ) << "Expected nullptr due to failed cast from Base to Derived.";
  delete base;   // Clean up allocated memory
}

// Test for dynamicCast with reference
TEST( DynamicCastTests, Reference_Casting_Success )
{
  DerivedFinalClass derived;
  BaseClass & base_ref = derived;
  DerivedFinalClass & derived_ref = geos::dynamicCast< DerivedFinalClass & >( base_ref );
  ASSERT_EQ( &derived_ref, &derived ) << "Expected successful cast from Base to Derived.";
}

TEST( DynamicCastTests, Reference_Casting_Failure )
{
  BaseClass base;
  BaseClass & base_ref = base;

  BaseClass & derived_base_ref = geos::dynamicCast< BaseClass & >( base_ref );
  ASSERT_EQ( &derived_base_ref, &base ) << "Expected successful cast from Base to Base.";

}


//// Typed test for geos wrapper
//template< typename T >
//class WrapperMock : public ::testing::Test
//{
//public:
//  WrapperMock():
//    m_node(),
//    m_group( "root", m_node ),
//    m_wrapper( "wrapper", m_group ),
//    m_wrapperBase( m_wrapper )
//  {}
//
//  void testDynamicCastWithPointer( )
//  {
//    {
//      WrapperBase * base_pointer = &m_wrapperBase;
//      Wrapper< T > * derived = geos::dynamicCast< Wrapper< T > * >( base_pointer );
//      ASSERT_NE( derived, nullptr ) << "Expected successful cast from Base to Derived.";
//    }
//    {
//      WrapperBase * base_pointer = &m_wrapperBase;
//      WrapperBase * derived = geos::dynamicCast< WrapperBase * >( base_pointer );
//      ASSERT_NE( derived, nullptr ) << "Expected successful cast from Base to Base.";
//    }
//    {
//      Wrapper< T > * defived_pointer = &m_wrapper;
//      Wrapper< T > * derived = geos::dynamicCast< Wrapper< T > * >( defived_pointer );
//      ASSERT_NE( derived, nullptr ) << "Expected successful cast from Derived to Derived.";
//    }
//  }
//
//  void testDynamicCastWithReference( )
//  {
//    {
//      WrapperBase & base_reference = m_wrapperBase;
//      Wrapper< T > & derived = geos::dynamicCast< Wrapper< T > & >( base_reference );
//      ASSERT_EQ( &derived, &base_reference ) << "Expected successful cast from Base to Derived.";
//    }
//    {
//      WrapperBase & base_reference = m_wrapperBase;
//      WrapperBase & derived = geos::dynamicCast< WrapperBase & >( base_reference );
//      ASSERT_EQ( &derived, &base_reference ) << "Expected successful cast from Base to Base.";
//    }
//    {
//      Wrapper< T > & defived_reference = m_wrapper;
//      Wrapper< T > & derived = geos::dynamicCast< Wrapper< T > & >( defived_reference );
//      ASSERT_EQ( &derived, &defived_reference ) << "Expected successful cast from Derived to Derived.";
//    }
//  }
//
//
//private:
//  conduit::Node m_node;
//  Group m_group;
//  Wrapper< T > m_wrapper;
//  WrapperBase & m_wrapperBase;
//};
//
//using WrapperMockTypes = ::testing::Types< int, array1d< real64 >, array1d< array1d< int > >, void *, std::function< void (void) > >;
//TYPED_TEST_SUITE( WrapperMock, WrapperMockTypes, );
//
//TYPED_TEST( WrapperMock, DynamicCastWithPointer )
//{
//  this->testDynamicCastWithPointer( );
//}
//
//TYPED_TEST( WrapperMock, DynamicCastWithReference )
//{
//  this->testDynamicCastWithReference( );
//}
//
//// Test Regex constructor
//TEST( RegexTests, Constructor )
//{
//  geos::Regex regex( "^[0-9]+$", "Input must be a number." );
//  ASSERT_EQ( regex.m_regexStr, "^[0-9]+$" ) << "Regex string is incorrect.";
//  ASSERT_EQ( regex.m_formatDescription, "Input must be a number." ) << "Format description is incorrect.";
//}
//
//TEST( RtTypesTests, GetTypeName )
//{
//  {
//    std::type_index typeIndex( typeid(Base));
//    auto typeName = geos::rtTypes::getTypeName( typeIndex );
//    EXPECT_EQ( typeName, std::string( "Base" ));  // Expected Base
//  }
//  {
//    std::type_index typeIndex( typeid(Derived));
//    auto typeName = geos::rtTypes::getTypeName( typeIndex );
//    EXPECT_EQ( typeName, std::string( "Derived" ));  // Expected Derived
//  }
//}
//
//// Additional tests to validate the functionality of getTypeRegex
//TEST( RtTypesTests, GetTypeRegex_Default )
//{
//  geos::Regex regex = geos::rtTypes::getTypeRegex< int >(); // Assuming int has a default regex defined
//  ASSERT_NE( regex.m_regexStr.empty(), true ) << "Expected non-empty regex for int.";
//}

int main( int argc, char * *argv )
{
  testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
