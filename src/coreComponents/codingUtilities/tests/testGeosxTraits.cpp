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
#include "codingUtilities/traits.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geos;
using namespace geos::traits;



IS_VALID_EXPRESSION( HasOperatorPlusEquals, T, std::declval< T & >() += T() );
TEST( testGeosxTraits, HasOperatorPlusEquals )
{
  static_assert( HasOperatorPlusEquals< int >, "Should be true." );
  static_assert( HasOperatorPlusEquals< string >, "Should be true." );

  static_assert( !HasOperatorPlusEquals< int const >, "Should be false." );
  static_assert( !HasOperatorPlusEquals< std::vector< int > >, "Should be false." );
}

IS_VALID_EXPRESSION_2( HasOperatorPlusEquals2, T, U, std::declval< T & >() += std::declval< U >() );
TEST( testGeosxTraits, HasOperatorPlusEquals2 )
{
  static_assert( HasOperatorPlusEquals2< double, int >, "Should be true." );
  static_assert( HasOperatorPlusEquals2< double, int const >, "Should be true." );
  static_assert( HasOperatorPlusEquals2< string, string >, "Should be true." );

  static_assert( !HasOperatorPlusEquals2< double const, int >, "Should be false." );
  static_assert( !HasOperatorPlusEquals2< R1Tensor, string >, "Should be false." );
  static_assert( !HasOperatorPlusEquals2< std::vector< int >, int >, "Should be false." );
}

HAS_MEMBER_FUNCTION_NO_RTYPE( at, 55 );
TEST( testGeosxTraits, HasMemberFunction_at )
{
  static_assert( HasMemberFunction_at< string >, "Should be true." );
  static_assert( HasMemberFunction_at< std::vector< int > >, "Should be true." );
  static_assert( HasMemberFunction_at< std::vector< double > >, "Should be true." );
  static_assert( HasMemberFunction_at< std::map< int, string > >, "Should be true." );
  static_assert( HasMemberFunction_at< std::unordered_map< int, string > >, "Should be true." );

  static_assert( !HasMemberFunction_at< int >, "Should be false." );
  static_assert( !HasMemberFunction_at< std::map< string, string > >, "Should be false." );
  static_assert( !HasMemberFunction_at< array1d< localIndex > >, "Should be false." );
}

HAS_MEMBER_FUNCTION_NO_RTYPE( insert,
                              std::declval< typename CLASS::const_iterator >(),
                              std::declval< typename CLASS::value_type >() );
TEST( testGeosxTraits, HasMemberFunction_insert )
{
  static_assert( HasMemberFunction_insert< std::vector< int > >, "Should be true." );
  static_assert( HasMemberFunction_insert< std::map< string, int > >, "Should be true." );
  static_assert( HasMemberFunction_insert< std::list< std::vector< int > > >, "Should be true." );

  static_assert( !HasMemberFunction_insert< int >, "Should be false." );
  static_assert( !HasMemberFunction_insert< std::array< int, 5 > >, "Should be false." );
}


TEST( testGeosxTraits, HasAlias_value_type )
{
  static_assert( HasAlias_value_type< array1d< double > >, "Should be true." );
  static_assert( HasAlias_value_type< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( HasAlias_value_type< SortedArray< string > >, "Should be true." );
  static_assert( HasAlias_value_type< std::vector< int > >, "Should be true." );
  static_assert( HasAlias_value_type< std::map< string, string > >, "Should be true." );

  static_assert( !HasAlias_value_type< int >, "Should be false." );
  static_assert( !HasAlias_value_type< double >, "Should be false." );
}

TEST( testGeosxTraits, Pointer )
{
  static_assert( std::is_same< Pointer< int >, int * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< R1Tensor >, R1Tensor * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< std::vector< double > >, double * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< string >, char * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< array3d< string > >, string * >::value, "Should be true." );
  static_assert( std::is_same< Pointer< SortedArray< float > >, float const * >::value, "Should be true." );

  static_assert( std::is_same< ConstPointer< int >, int const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< R1Tensor >, R1Tensor const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< std::vector< double > >, double const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< string >, char const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< array3d< string > >, string const * >::value, "Should be true." );
  static_assert( std::is_same< ConstPointer< SortedArray< float > >, float const * >::value, "Should be true." );
}

TEST( testGeosxTraits, HasMemberFunction_data )
{
  static_assert( HasMemberFunction_data< array1d< double > >, "Should be true." );
  static_assert( HasMemberFunction_data< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( HasMemberFunction_data< std::vector< int > >, "Should be true." );
  static_assert( HasMemberFunction_data< SortedArray< string > >, "Should be true." );

  static_assert( !HasMemberFunction_data< std::map< string, string > >, "Should be false." );
  static_assert( !HasMemberFunction_data< int >, "Should be false." );
  static_assert( !HasMemberFunction_data< double >, "Should be false." );
}

TEST( testGeosxTraits, HasMemberFunction_size )
{
  static_assert( HasMemberFunction_size< array1d< double > >, "Should be true." );
  static_assert( HasMemberFunction_size< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( HasMemberFunction_size< SortedArray< string > >, "Should be true." );
  static_assert( HasMemberFunction_size< std::vector< int > >, "Should be true." );
  static_assert( HasMemberFunction_size< std::map< string, string > >, "Should be true." );

  static_assert( !HasMemberFunction_size< int >, "Should be false." );
  static_assert( !HasMemberFunction_size< double >, "Should be false." );
}

TEST( testGeosxTraits, HasMemberFunction_resize )
{
  static_assert( HasMemberFunction_resize< array1d< double > >, "Should be true." );
  static_assert( HasMemberFunction_resize< array5d< array1d< R1Tensor > > >, "Should be true." );
  static_assert( HasMemberFunction_resize< std::vector< int > >, "Should be true." );
  static_assert( HasMemberFunction_resize< ArrayOfArrays< int > >, "Should be true." );
  static_assert( HasMemberFunction_resize< ArrayOfSets< int > >, "Should be true." );

  static_assert( !HasMemberFunction_resize< SortedArray< string > >, "Should be false." );
  static_assert( !HasMemberFunction_resize< std::map< string, string > >, "Should be false." );
  static_assert( !HasMemberFunction_resize< int >, "Should be false." );
  static_assert( !HasMemberFunction_resize< double >, "Should be false." );

}

TEST( testGeosxTraits, CanStreamInto )
{
  static_assert( CanStreamInto< std::istringstream, int >, "Should be true." );
  static_assert( CanStreamInto< std::istringstream, double >, "Should be true." );
  static_assert( CanStreamInto< std::istringstream, string >, "Should be true." );

  static_assert( !CanStreamInto< std::istringstream, array1d< double > >, "Should be false." );
  static_assert( !CanStreamInto< std::istringstream, std::vector< int > >, "Should be false." );
}

TEST( testGeosxTraits, is_tensorT )
{
  static_assert( is_tensorT< R1Tensor >, "Should be true." );

  static_assert( !is_tensorT< int >, "Should be false." );
  static_assert( !is_tensorT< double >, "Should be false." );
  static_assert( !is_tensorT< void >, "Should be false." );
}


TEST( testGeosxTraits, is_string )
{
  static_assert( is_string< string >, "Should be true." );

  static_assert( !is_string< int >, "Should be false." );
  static_assert( !is_string< double >, "Should be false." );
  static_assert( !is_string< void >, "Should be false." );
}

TEST( testGeosxTraits, is_array )
{
  static_assert( is_array< array1d< int > >, "Should be true." );

  static_assert( !is_array< int >, "Should be false." );
  static_assert( !is_array< double >, "Should be false." );
  static_assert( !is_array< void >, "Should be false." );
}


struct Foo
{
  Foo & operator=( Foo const & ) = default;
  int m_value;
};

struct Bar
{
  Bar & operator=( Bar const & ) = delete;
  int m_value;
};

enum enumFoo
{
  blah,
  yada
};

enum class enumClassFoo
{
  blah,
  yada
};

TEST( testGeosxTraits, hasCopyAssignment )
{
  static_assert( hasCopyAssignmentOp< enumFoo >, "Should be true." );
  static_assert( hasCopyAssignmentOp< enumClassFoo >, "Should be true." );
  static_assert( hasCopyAssignmentOp< int >, "Should be true." );
  static_assert( hasCopyAssignmentOp< real64 >, "Should be true." );
  static_assert( hasCopyAssignmentOp< Foo >, "Should be true." );
  static_assert( !hasCopyAssignmentOp< Bar >, "Should be false." );
}

struct A
{};

struct B
{};

struct C
{};

TEST( testGeosxTraits, typeListIndex )
{
  static_assert( type_list_index< A, std::tuple< A, B > > == 0, "Should be true." );
  static_assert( type_list_index< B, std::tuple< A, B > > == 1, "Should be true." );
  static_assert( type_list_index< C, std::tuple< A, B > > == 2, "Should be true." );
  static_assert( type_list_index< A, std::tuple<> > == 0, "Should be true." );
}
