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
 * @file wrapperHelpers.hpp
 */

#ifndef GEOS_DATAREPOSITORY_WRAPPERHELPERS_HPP_
#define GEOS_DATAREPOSITORY_WRAPPERHELPERS_HPP_


/// Enables verbose logging of restart output
#define RESTART_TYPE_LOGGING 0

// Source includes
#include "BufferOps.hpp"
#include "BufferOpsDevice.hpp"
#include "DefaultValue.hpp"
#include "ConduitRestart.hpp"
#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"
#include "common/Span.hpp"
#include "codingUtilities/traits.hpp"

#if defined(GEOS_USE_PYGEOSX)
#include "LvArray/src/python/python.hpp"
#endif

// TPL includes
#include <conduit.hpp>

// System includes
#include <cstring>

#if RESTART_TYPE_LOGGING
#include <unordered_set>
#endif

namespace geos
{
namespace dataRepository
{
namespace wrapperHelpers
{
namespace internal
{

inline void logOutputType( string const & typeString, string const & msg )
{
#if RESTART_TYPE_LOGGING
  static std::unordered_set< string > m_types;

  if( !m_types.count( typeString ) )
  {
    m_types.insert( typeString );
    GEOS_LOG( msg << typeString );
  }
#else
  GEOS_DEBUG_VAR( typeString );
  GEOS_DEBUG_VAR( msg );
#endif
}

template< typename T, typename ... INDICES >
string getIndicesToComponent( T const &, int const component, INDICES const ... existingIndices )
{
  GEOS_ERROR_IF_NE( component, 0 );
  return LvArray::indexing::getIndexString( existingIndices ... );
}

template< typename ... INDICES >
string getIndicesToComponent( R1Tensor const &, int const component, INDICES const ... existingIndices )
{ return LvArray::indexing::getIndexString( existingIndices ..., component ); }

template< typename T >
T const * getPointerToComponent( T const & var, int const component )
{
  GEOS_ERROR_IF_NE( component, 0 );
  return &var;
}

inline
real64 const * getPointerToComponent( R1Tensor const & var, int const component )
{
  GEOS_ERROR_IF_GE( component, 3 );
  return &var[ component ];
}

} // namespace internal

template< typename T >
class ArrayDimLabels
{
public:

  void set( integer const, Span< string const > )
  {
    GEOS_ERROR( "Dimension labels are only available in Array wrappers" );
  }

  Span< string const > get( integer const ) const
  {
    GEOS_ERROR( "Dimension labels are only available in Array wrappers" );
    return {};
  }
};

template< typename T, int NDIM, typename PERM >
class ArrayDimLabels< Array< T, NDIM, PERM > >
{
public:

  void set( integer const dim, Span< string const > labels )
  {
    GEOS_ERROR_IF_LT( dim, 0 );
    GEOS_ERROR_IF_GE( dim, NDIM );
    m_values[dim].resize( labels.size() );
    std::copy( labels.begin(), labels.end(), m_values[dim].begin() );
  }

  Span< string const > get( integer const dim ) const
  {
    GEOS_ERROR_IF_LT( dim, 0 );
    GEOS_ERROR_IF_GE( dim, NDIM );
    return { m_values[dim].begin(), m_values[dim].end() };
  }

private:

  string_array m_values[NDIM]{};
};

template< typename T >
inline std::enable_if_t< traits::HasMemberFunction_size< T >, size_t >
size( T const & value )
{ return value.size(); }

template< typename T >
inline std::enable_if_t< !traits::HasMemberFunction_size< T >, size_t >
size( T const & GEOS_UNUSED_PARAM( value ) )
{ return 1; }


inline char *
dataPtr( string & var )
{ return const_cast< char * >( var.data() ); }

inline char *
dataPtr( Path & var )
{ return const_cast< char * >( var.data() ); }

template< typename T >
inline std::enable_if_t< traits::HasMemberFunction_data< T >, typename traits::Pointer< T > >
dataPtr( T & value )
{ return value.data(); }

template< typename T >
inline std::enable_if_t< !traits::HasMemberFunction_data< T >, typename traits::Pointer< T > >
dataPtr( T & value )
{ return &value; }

template< class T >
inline typename traits::ConstPointer< T >
dataPtr( T const & value )
{ return dataPtr( const_cast< T & >( value ) ); }


template< typename T >
inline std::enable_if_t< traits::HasMemberFunction_resize< T > >
resize( T & value, localIndex const newSize )
{ value.resize( newSize ); }

template< typename T >
inline std::enable_if_t< !traits::HasMemberFunction_resize< T > >
resize( T & GEOS_UNUSED_PARAM( value ),
        localIndex const GEOS_UNUSED_PARAM( newSize ) )
{}


template< typename T, int NDIM, typename PERMUTATION >
inline std::enable_if_t< DefaultValue< Array< T, NDIM, PERMUTATION > >::has_default_value >
resizeDefault( Array< T, NDIM, PERMUTATION > & value,
               localIndex const newSize,
               DefaultValue< Array< T, NDIM, PERMUTATION > > const & defaultValue )
{ value.resizeDefault( newSize, defaultValue.value ); }

template< typename T >
inline void
resizeDefault( T & value, localIndex const newSize, DefaultValue< T > const & GEOS_UNUSED_PARAM( defaultValue ) )
{ resize( value, newSize ); }


template< typename T, int NDIM, typename PERMUTATION >
inline void
resizeDimensions( Array< T, NDIM, PERMUTATION > & value, int num_dims, localIndex const * const dims )
{ value.resize( num_dims, dims ); }

template< typename T >
inline void
resizeDimensions( T & value, int num_dims, localIndex const * const dims )
{
  if( num_dims != 1 )
  {
    GEOS_ERROR( "Data is not multidimensional" );
    return;
  }
  resize( value, dims[ 0 ] );
}


template< typename T >
inline localIndex
byteSizeOfElement()
{ return sizeof( *dataPtr( std::declval< T >() ) ); }


template< typename T >
inline size_t
byteSize( T const & value )
{ return wrapperHelpers::size( value ) * byteSizeOfElement< T >(); }


template< typename T >
inline localIndex
numElementsFromByteSize( localIndex const byteSize )
{
  GEOS_ERROR_IF_NE( byteSize % byteSizeOfElement< T >(), 0 );
  return byteSize / byteSizeOfElement< T >();
}


template< typename T >
std::enable_if_t< traits::HasMemberFunction_reserve< T > >
reserve( T & value, localIndex const newCapacity )
{ value.reserve( newCapacity ); }

template< typename T >
std::enable_if_t< !traits::HasMemberFunction_reserve< T > >
reserve( T & GEOS_UNUSED_PARAM( value ), localIndex const GEOS_UNUSED_PARAM( newCapacity ) )
{}


template< typename T >
std::enable_if_t< traits::HasMemberFunction_capacity< T const >, localIndex >
capacity( T const & value )
{ return value.capacity(); }

template< typename T >
std::enable_if_t< !traits::HasMemberFunction_capacity< T const >, localIndex >
capacity( T const & value )
{ return wrapperHelpers::size( value ); }



template< typename T >
std::enable_if_t< traits::HasMemberFunction_setName< T > >
setName( T & value, string const & name )
{ value.setName( name ); }

template< typename T >
std::enable_if_t< !traits::HasMemberFunction_setName< T > >
setName( T & GEOS_UNUSED_PARAM( value ), string const & GEOS_UNUSED_PARAM( name ) )
{}

template< typename T >
std::enable_if_t< traits::HasMemberFunction_move< T > >
move( T & value, LvArray::MemorySpace const space, bool const touch )
{ value.move( space, touch ); }

template< typename T >
std::enable_if_t< !traits::HasMemberFunction_move< T > >
move( T & GEOS_UNUSED_PARAM( value ),
      LvArray::MemorySpace const GEOS_UNUSED_PARAM( space ),
      bool const GEOS_UNUSED_PARAM( touch ) )
{}

// This is for an object that needs to be packed.
template< typename T >
std::enable_if_t< !bufferOps::can_memcpy< typename traits::Pointer< T > > >
pushDataToConduitNode( T const & var, conduit::Node & node )
{
  internal::logOutputType( LvArray::system::demangleType( var ), "Packing for output: " );

  // Get the number of bytes in the packed object.
  localIndex const byteSize = bufferOps::PackSize( var );

  // Create a conduit data type that describes the array.
  conduit::DataType const dtype( conduitTypeInfo< buffer_unit_type >::id, byteSize );

  // Allocate the array in the "__values__" child.
  conduit::Node & valuesNode = node[ "__values__" ];
  valuesNode.set( dtype );

  // Get the pointer to the array and pack the object into it.
  buffer_unit_type * buffer = valuesNode.value();
  bufferOps::Pack< true >( buffer, var );
}

// This is for an object that needs to be packed.
template< typename T >
std::enable_if_t< !bufferOps::can_memcpy< typename traits::Pointer< T > > >
pullDataFromConduitNode( T & var, conduit::Node const & node )
{
  conduit::Node const & valuesNode = node.fetch_existing( "__values__" );

  // Get the number of bytes in the array and a pointer to the array.
  localIndex const byteSize = valuesNode.dtype().number_of_elements();
  buffer_unit_type const * buffer = valuesNode.value();

  // Unpack the object from the array.
  localIndex const bytesRead = bufferOps::Unpack( buffer, var );
  GEOS_ERROR_IF_NE( bytesRead, byteSize );
}

// This is for an string since the type of char is different on different platforms :(.
inline
void
pushDataToConduitNode( string const & var, conduit::Node & node )
{
  internal::logOutputType( LvArray::system::demangleType( var ), "Output via external pointer: " );

  constexpr int conduitTypeID = conduitTypeInfo< signed char >::id;
  conduit::DataType const dtype( conduitTypeID, var.size() );

  signed char * const ptr = const_cast< signed char * >( reinterpret_cast< signed char const * >( var.data() ) );
  node[ "__values__" ].set_external( dtype, ptr );
}

// This is for Path since it derives from string. See overload for string.
inline
void
pushDataToConduitNode( Path const & var, conduit::Node & node )
{
  pushDataToConduitNode( static_cast< string const & >(var), node );
}

// This is for an object that doesn't need to be packed but isn't an LvArray.
template< typename T >
std::enable_if_t< bufferOps::can_memcpy< typename traits::Pointer< T > > >
pushDataToConduitNode( T const & var, conduit::Node & node )
{
  internal::logOutputType( LvArray::system::demangleType( var ), "Output via external pointer: " );

  constexpr int conduitTypeID = conduitTypeInfo< typename traits::Pointer< T > >::id;
  constexpr int sizeofConduitType = conduitTypeInfo< typename traits::Pointer< T > >::sizeOfConduitType;
  localIndex const numBytes = byteSize( var );
  conduit::DataType const dtype( conduitTypeID, numBytes / sizeofConduitType );

  void * const ptr = const_cast< void * >( static_cast< void const * >( dataPtr( var ) ) );
  node[ "__values__" ].set_external( dtype, ptr );
}

// This is for an object that doesn't need to be packed but isn't an LvArray or a SortedArray.
template< typename T >
std::enable_if_t< bufferOps::can_memcpy< typename traits::Pointer< T > > >
pullDataFromConduitNode( T & var, conduit::Node const & node )
{
  conduit::Node const & valuesNode = node.fetch_existing( "__values__" );

  localIndex const byteSize = LvArray::integerConversion< localIndex >( valuesNode.dtype().strided_bytes() );
  localIndex const numElements = numElementsFromByteSize< T >( byteSize );

  resize( var, numElements );

  std::memcpy( dataPtr( var ), valuesNode.data_ptr(), byteSize );
}

// This is for a SortedArray that doesn't need to be packed.
template< typename T >
std::enable_if_t< bufferOps::can_memcpy< T > >
pullDataFromConduitNode( SortedArray< T > & var, conduit::Node const & node )
{
  conduit::Node const & valuesNode = node.fetch_existing( "__values__" );

  localIndex const byteSize = LvArray::integerConversion< localIndex >( valuesNode.dtype().strided_bytes() );
  localIndex const numElements = numElementsFromByteSize< T >( byteSize );

  T const * const values = reinterpret_cast< T const * >( valuesNode.data_ptr() );
  var.insert( values, values + numElements );
}


// This is an LvArray that doesn't need to be packed.
template< typename T, int NDIM, typename PERMUTATION >
std::enable_if_t< bufferOps::can_memcpy< T > >
pushDataToConduitNode( Array< T, NDIM, PERMUTATION > const & var,
                       conduit::Node & node )
{
  internal::logOutputType( LvArray::system::demangleType( var ), "Output array via external pointer: " );

  // Push the data into conduit
  constexpr int conduitTypeID = conduitTypeInfo< T >::id;
  constexpr int sizeofConduitType = conduitTypeInfo< T >::sizeOfConduitType;
  conduit::DataType const dtype( conduitTypeID, var.size() * sizeof( T ) / sizeofConduitType );
  void * const ptr = const_cast< void * >( static_cast< void const * >( var.data() ) );
  node[ "__values__" ].set_external( dtype, ptr );

  // Create a copy of the dimensions
  camp::idx_t temp[ NDIM + 1 ];
  for( int i = 0; i < NDIM; ++i )
  {
    temp[ i ] = var.size( i );
  }

  // If T is something like a Tensor than there is an extra implicit dimension.
  constexpr int const implicitDimensionLength = conduitTypeInfo< T >::numConduitValues;
  constexpr bool const hasImplicitDimension = implicitDimensionLength != 1;
  constexpr int totalNumDimensions = NDIM + hasImplicitDimension;
  if( hasImplicitDimension )
  {
    temp[ NDIM ] = implicitDimensionLength;
  }

  // push the dimensions into the node
  conduit::DataType const dimensionType( conduitTypeInfo< camp::idx_t >::id, totalNumDimensions );
  node[ "__dimensions__" ].set( dimensionType, temp );

  // Create a copy of the permutation
  constexpr std::array< camp::idx_t, NDIM > const perm = RAJA::as_array< PERMUTATION >::get();
  for( int i = 0; i < NDIM; ++i )
  {
    temp[ i ] = perm[ i ];
  }

  if( hasImplicitDimension )
  {
    temp[ NDIM ] = NDIM;
  }

  node[ "__permutation__" ].set( dimensionType, temp );
}

// This is an LvArray that doesn't need to be packed.
template< typename T, int NDIM, typename PERMUTATION >
std::enable_if_t< bufferOps::can_memcpy< T > >
pullDataFromConduitNode( Array< T, NDIM, PERMUTATION > & var,
                         conduit::Node const & node )
{
  // Get the number of dimensions written out, accounting for an implicit dimension and the permutation.
  constexpr int const implicitDimensionLength = conduitTypeInfo< T >::numConduitValues;
  constexpr bool const hasImplicitDimension = implicitDimensionLength != 1;
  constexpr int totalNumDimensions = NDIM + hasImplicitDimension;

  // Check that the permutations match.
  conduit::Node const & permutationNode = node.fetch_existing( "__permutation__" );
  GEOS_ERROR_IF_NE( permutationNode.dtype().number_of_elements(), totalNumDimensions );

  constexpr std::array< camp::idx_t, NDIM > const perm = RAJA::as_array< PERMUTATION >::get();
  camp::idx_t const * const permFromConduit = permutationNode.value();
  for( int i = 0; i < NDIM; ++i )
  {
    GEOS_ERROR_IF_NE_MSG( permFromConduit[ i ], perm[ i ],
                          "The permutation of the data in conduit and the provided Array don't match." );
  }

  if( hasImplicitDimension )
  {
    GEOS_ERROR_IF_NE_MSG( permFromConduit[ NDIM ], NDIM,
                          "The permutation of the data in conduit and the provided Array don't match." );
  }

  // Now pull out the dimensions and resize the array.
  conduit::Node const & dimensionNode = node.fetch_existing( "__dimensions__" );
  GEOS_ERROR_IF_NE( dimensionNode.dtype().number_of_elements(), totalNumDimensions );
  camp::idx_t const * const dims = dimensionNode.value();

  if( hasImplicitDimension )
  {
    GEOS_ERROR_IF_NE( dims[ NDIM ], implicitDimensionLength );
  }

  var.resize( NDIM, dims );

  // Finally memcpy
  conduit::Node const & valuesNode = node.fetch_existing( "__values__" );
  localIndex numBytesFromArray =  var.size() * sizeof( T );
  GEOS_ERROR_IF_NE( numBytesFromArray, valuesNode.dtype().strided_bytes() );
  std::memcpy( var.data(), valuesNode.data_ptr(), numBytesFromArray );
}



template< typename T, typename INDEX_TYPE >
std::enable_if_t< bufferOps::can_memcpy< T > >
pushDataToConduitNode( ArrayOfArrays< T, INDEX_TYPE > const & var2,
                       conduit::Node & node )
{
  ArrayOfArraysView< T const, INDEX_TYPE > const & var = var2.toViewConst();
  internal::logOutputType( LvArray::system::demangleType( var ), "Output array via external pointer: " );

  // ArrayOfArray::m_numArrays
  INDEX_TYPE const numArrays = var.size();
  conduit::DataType const numArraysType( conduitTypeInfo< INDEX_TYPE >::id, 1 );
  node[ "__numberOfArrays__" ].set( numArraysType, const_cast< void * >( static_cast< void const * >(&numArrays) ) );

  // ArrayOfArray::m_offsets
  INDEX_TYPE const * const offsets = var.getOffsets();
  conduit::DataType const offsetsType( conduitTypeInfo< INDEX_TYPE >::id, numArrays+1 );
  node[ "__offsets__" ].set_external( offsetsType, const_cast< void * >( static_cast< void const * >( offsets ) ) );

  // ArrayOfArray::m_sizes
  INDEX_TYPE const * const sizes = var.getSizes();
  conduit::DataType const sizesType( conduitTypeInfo< INDEX_TYPE >::id, numArrays );
  node[ "__sizes__" ].set_external( sizesType, const_cast< void * >( static_cast< void const * >( sizes ) ) );

  // Push the data into conduit
  constexpr int conduitTypeID = conduitTypeInfo< T >::id;
  constexpr int sizeofConduitType = conduitTypeInfo< T >::sizeOfConduitType;
  conduit::DataType const dtype( conduitTypeID, offsets[numArrays] * sizeof( T ) / sizeofConduitType );
  void * const ptr = const_cast< void * >( static_cast< void const * >( var.getValues() ) );
  node[ "__values__" ].set_external( dtype, ptr );
}

// This is an LvArray that doesn't need to be packed.
template< typename T, typename INDEX_TYPE >
std::enable_if_t< bufferOps::can_memcpy< T > >
pullDataFromConduitNode( ArrayOfArrays< T, INDEX_TYPE > & var,
                         conduit::Node const & node )
{


  conduit::Node const & numArraysNode = node.fetch_existing( "__numberOfArrays__" );
  INDEX_TYPE const * const numArrays = numArraysNode.value();

  conduit::Node const & offsetsNode = node.fetch_existing( "__offsets__" );
  conduit::DataType const & offsetsDataType = offsetsNode.dtype();
  INDEX_TYPE const * const offsets = offsetsNode.value();
  INDEX_TYPE const sizeOffsets = offsetsDataType.number_of_elements();

  conduit::Node const & sizesNode = node.fetch_existing( "__sizes__" );
  conduit::DataType const & sizesDataType = sizesNode.dtype();
  INDEX_TYPE const * const sizes = sizesNode.value();
  INDEX_TYPE const sizeSizes = sizesDataType.number_of_elements();

  GEOS_ERROR_IF_NE( *numArrays, sizeSizes );
  GEOS_ERROR_IF_NE( *numArrays+1, sizeOffsets );

  // values node
  conduit::Node const & valuesNode = node.fetch_existing( "__values__" );
  conduit::DataType const & valuesDataType = valuesNode.dtype();
  const INDEX_TYPE valuesSize = valuesDataType.number_of_elements();

  var.resize( *numArrays, valuesSize/(*numArrays) );
  localIndex allocatedSize = 0;
  for( INDEX_TYPE i = 0; i < *numArrays; ++i )
  {
    INDEX_TYPE const arrayAllocation = offsets[i+1] - offsets[i];
    var.setCapacityOfArray( i, arrayAllocation );
    var.resizeArray( i, sizes[ i ] );
    allocatedSize += arrayAllocation;
  }
  ArrayOfArraysView< T const, INDEX_TYPE > const & varView = var.toViewConst();

  GEOS_ERROR_IF_NE( valuesSize, allocatedSize );
  GEOS_ERROR_IF_NE( allocatedSize, offsets[sizeOffsets] );

  INDEX_TYPE const * const varOffsets = varView.getOffsets();
  INDEX_TYPE const * const varSizes = varView.getSizes();

  GEOS_ERROR_IF_NE( varOffsets[0], offsets[0] );
  for( INDEX_TYPE i = 0; i<*numArrays; ++i )
  {
    GEOS_ERROR_IF_NE( varOffsets[i+1], offsets[i+1] );
    GEOS_ERROR_IF_NE( varSizes[i], sizes[i] );
  }


  localIndex numBytesFromArray =  allocatedSize * sizeof( T );
  // GEOS_ERROR_IF_NE( numBytesFromArray, valuesNode.dtype().strided_bytes() );
  std::memcpy( &var( 0, 0 ), valuesNode.data_ptr(), numBytesFromArray );
}



template< typename T >
void pushDataToConduitNode( InterObjectRelation< T > const & var,
                            conduit::Node & node )
{return pushDataToConduitNode( var.base(), node ); }

template< typename T >
void pullDataFromConduitNode( InterObjectRelation< T > & var,
                              conduit::Node const & node )
{ return pullDataFromConduitNode( var.base(), node ); }


/// TODO: Remove this function once https://github.com/visit-dav/visit/issues/4637 is fixed and released.
template< typename T, int NDIM, int USD >
std::enable_if_t< std::is_arithmetic< T >::value || traits::is_tensorT< T > >
addBlueprintField( ArrayView< T const, NDIM, USD > const & var,
                   conduit::Node & fields,
                   string const & fieldName,
                   string const & topology,
                   std::vector< string > const & componentNames )
{
  GEOS_ERROR_IF_LE( var.size(), 0 );

  using ConduitType = typename conduitTypeInfo< T >::type;
  constexpr int conduitTypeID = conduitTypeInfo< T >::id;
  constexpr int numComponentsPerValue = conduitTypeInfo< T >::numConduitValues;

  localIndex const totalNumberOfComponents = numComponentsPerValue * var.size() / var.size( 0 );
  if( !componentNames.empty() )
  {
    GEOS_ERROR_IF_NE( localIndex( componentNames.size() ), totalNumberOfComponents );
  }

  var.move( hostMemorySpace, false );

  conduit::DataType dtype( conduitTypeID, var.size( 0 ) );
  dtype.set_stride( sizeof( ConduitType ) * numComponentsPerValue * var.strides()[ 0 ] );

  localIndex curComponent = 0;
  LvArray::forValuesInSliceWithIndices( var[ 0 ], [&fields, &fieldName, &topology, &componentNames, totalNumberOfComponents, &dtype, &curComponent]
                                          ( T const & val, auto const ... indices )
  {
    for( int i = 0; i < numComponentsPerValue; ++i )
    {
      string name;
      if( totalNumberOfComponents == 1 )
      {
        name = fieldName;
      }
      else if( componentNames.empty() )
      {
        string indexString = internal::getIndicesToComponent( val, i, indices ... );
        indexString.erase( indexString.begin() );
        indexString.pop_back();
        indexString.pop_back();
        name = fieldName + indexString;
      }
      else
      {
        name = componentNames[ curComponent++ ];
      }

      conduit::Node & field = fields[ name ];
      field[ "association" ] = "element";
      field[ "volume_dependent" ] = "false";
      field[ "topology" ] = topology;

      void const * pointer = internal::getPointerToComponent( val, i );
      field[ "values" ].set_external( dtype, const_cast< void * >( pointer ) );
    }
  } );
}

template< typename T >
void addBlueprintField( T const &,
                        conduit::Node & fields,
                        string const &,
                        string const &,
                        std::vector< string > const & )
{
  GEOS_ERROR( "Cannot create a mcarray out of " << LvArray::system::demangleType< T >() <<
              "\nWas trying to write it to " << fields.path() );
  GEOS_UNUSED_VAR( fields );
}

template< typename T, int NDIM, int USD >
std::enable_if_t< std::is_arithmetic< T >::value || traits::is_tensorT< T > >
populateMCArray( ArrayView< T const, NDIM, USD > const & var,
                 conduit::Node & node,
                 std::vector< string > const & componentNames )
{
  GEOS_ERROR_IF_LE( var.size(), 0 );

  using ConduitType = typename conduitTypeInfo< T >::type;
  constexpr int conduitTypeID = conduitTypeInfo< T >::id;
  constexpr int numComponentsPerValue = conduitTypeInfo< T >::numConduitValues;

  if( !componentNames.empty() )
  {
    GEOS_ERROR_IF_NE( localIndex( componentNames.size() ), numComponentsPerValue * var.size() / var.size( 0 ) );
  }

  var.move( hostMemorySpace, false );

  conduit::DataType dtype( conduitTypeID, var.size( 0 ) );
  dtype.set_stride( sizeof( ConduitType ) * numComponentsPerValue * var.strides()[ 0 ] );

  localIndex curComponent = 0;
  LvArray::forValuesInSliceWithIndices( var[ 0 ], [&componentNames, &node, &dtype, &curComponent]
                                          ( T const & val, auto const ... indices )
  {
    for( int i = 0; i < numComponentsPerValue; ++i )
    {
      string const name = componentNames.empty() ? internal::getIndicesToComponent( val, i, indices ... ) :
                          componentNames[ curComponent++ ];

      void const * pointer = internal::getPointerToComponent( val, i );
      node[ name ].set_external( dtype, const_cast< void * >( pointer ) );
    }
  } );
}

template< typename T >
void populateMCArray( T const &,
                      conduit::Node & node,
                      std::vector< string > const & )
{
  GEOS_ERROR( "Cannot create a mcarray out of " << LvArray::system::demangleType< T >() <<
              "\nWas trying to write it to " << node.path() );
  GEOS_UNUSED_VAR( node );
}

template< typename T >
std::enable_if_t< std::is_arithmetic< T >::value, std::unique_ptr< Array< T, 1 > > >
averageOverSecondDim( ArrayView< T const, 1, 0 > const & var )
{
  std::unique_ptr< Array< T, 1 > > ret = std::make_unique< Array< T, 1 > >();

  ret->resize( var.size() );
  ret->template setValues< serialPolicy >( var );

  return ret;
}

template< typename T, int NDIM, int USD >
std::enable_if_t< std::is_arithmetic< T >::value, std::unique_ptr< Array< T, NDIM - 1 > > >
averageOverSecondDim( ArrayView< T const, NDIM, USD > const & var )
{
  std::unique_ptr< Array< T, NDIM - 1 > > ret = std::make_unique< Array< T, NDIM - 1 > >();

  localIndex newDims[ NDIM - 1 ];
  newDims[ 0 ] = var.size( 0 );
  for( int i = 2; i < NDIM; ++i )
  {
    newDims[ i - 1 ] = var.size( i );
  }

  ret->resize( NDIM - 1, newDims );

  ArrayView< T, NDIM - 1 > const & output = *ret;

  localIndex const numSamples = var.size( 1 );
  forAll< serialPolicy >( var.size( 0 ), [var, numSamples, &output] ( localIndex const i )
  {
    LvArray::sumOverFirstDimension( var[ i ], output[ i ] );

    LvArray::forValuesInSlice( output[ i ], [numSamples] ( T & val )
    {
      val /= numSamples;
    } );
  } );

  return ret;
}

template< typename T >
std::unique_ptr< int > averageOverSecondDim( T const & )
{
  GEOS_ERROR( "Cannot average over the second dimension of " << LvArray::system::demangleType< T >() );
  return std::unique_ptr< int >( nullptr );
}

template< typename T, int NDIM, int USD >
int numArrayDims( ArrayView< T const, NDIM, USD > const & GEOS_UNUSED_PARAM( var ) )
{
  return NDIM;
}

template< typename T >
int numArrayDims( T const & GEOS_UNUSED_PARAM( var ) )
{
  return 0;
}

template< typename T, int NDIM, int USD >
localIndex numArrayComp( ArrayView< T const, NDIM, USD > const & var )
{
  return LvArray::indexing::multiplyAll< NDIM - 1 >( var.dims() + 1 );
}

template< typename T >
localIndex numArrayComp( ArrayView< T const, 1, 0 > const & GEOS_UNUSED_PARAM( var ) )
{
  return 1;
}

template< typename T >
localIndex numArrayComp( T const & GEOS_UNUSED_PARAM( var ) )
{
  return 0;
}

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_packable_by_index< T >, localIndex >
PackByIndex( buffer_unit_type * & buffer, T & var, IDX & idx )
{ return bufferOps::PackByIndex< DO_PACKING >( buffer, var, idx ); }

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_packable_by_index< T >, localIndex >
PackByIndex( buffer_unit_type * &, T &, IDX & )
{
  GEOS_ERROR( "Trying to pack data type (" << LvArray::system::demangleType< T >() << ") by index. Operation not supported." );
  return 0;
}

template< typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_packable_by_index< T >, localIndex >
UnpackByIndex( buffer_unit_type const * & buffer, T & var, IDX & idx )
{ return bufferOps::UnpackByIndex( buffer, var, idx ); }

template< typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_packable_by_index< T >, localIndex >
UnpackByIndex( buffer_unit_type const * &, T &, IDX & )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") by index. Operation not supported." );
  return 0;
}


template< bool DO_PACKING, typename T >
inline std::enable_if_t< bufferOps::is_container< T > || bufferOps::can_memcpy< T >, localIndex >
PackDevice( buffer_unit_type * & buffer, T const & var, parallelDeviceEvents & events )
{ return bufferOps::PackDevice< DO_PACKING >( buffer, var, events ); }


template< bool DO_PACKING, typename T >
inline std::enable_if_t< !bufferOps::is_container< T > && !bufferOps::can_memcpy< T >, localIndex >
PackDevice( buffer_unit_type * &, T const &, parallelDeviceEvents & )
{
  GEOS_ERROR( "Trying to pack data type (" << LvArray::system::demangleType< T >() << ") on device. Operation not supported." );
  return 0;
}

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
PackByIndexDevice( buffer_unit_type * & buffer, T const & var, IDX & idx, parallelDeviceEvents & events )
{ return bufferOps::PackByIndexDevice< DO_PACKING >( buffer, var, idx, events ); }

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
PackByIndexDevice( buffer_unit_type * &, T const &, IDX &, parallelDeviceEvents & )
{
  GEOS_ERROR( "Trying to pack data type (" << LvArray::system::demangleType< T >() << ") by index on device. Operation not supported." );
  return 0;
}

template< typename T >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
UnpackDevice( buffer_unit_type const * & buffer, T const & var, parallelDeviceEvents & events )
{ return bufferOps::UnpackDevice( buffer, var, events ); }

template< typename T >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
UnpackDevice( buffer_unit_type const * &, T const &, parallelDeviceEvents & )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") on device. Operation not supported." );
  return 0;
}

template< typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
UnpackByIndexDevice( buffer_unit_type const * & buffer, T const & var, IDX & idx, parallelDeviceEvents & events, MPI_Op op=MPI_REPLACE )
{ return bufferOps::UnpackByIndexDevice( buffer, var, idx, events, op ); }

template< typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
UnpackByIndexDevice( buffer_unit_type const * &, T &, IDX &, parallelDeviceEvents &, MPI_Op )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") by index on device. Operation not supported." );
  return 0;
}


template< bool DO_PACKING, typename T >
localIndex
PackDataDevice( buffer_unit_type * & buffer, T const & var, parallelDeviceEvents & events )
{ return bufferOps::PackDataDevice< DO_PACKING >( buffer, var, events ); }

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
PackDataByIndexDevice( buffer_unit_type * & buffer, T const & var, IDX & idx, parallelDeviceEvents & events )
{ return bufferOps::PackDataByIndexDevice< DO_PACKING >( buffer, var, idx, events ); }

template< bool DO_PACKING, typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
PackDataByIndexDevice( buffer_unit_type * &, T const &, IDX &, parallelDeviceEvents & )
{
  GEOS_ERROR( "Trying to pack data type (" << LvArray::system::demangleType< T >() << ") by index on device. Operation not supported." );
  return 0;
}

template< typename T >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
UnpackDataDevice( buffer_unit_type const * & buffer, T const & var, parallelDeviceEvents & events )
{ return bufferOps::UnpackDataDevice( buffer, var, events ); }

template< typename T >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
UnpackDataDevice( buffer_unit_type const * &, T const &, parallelDeviceEvents & )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") on device. Operation not supported." );
  return 0;
}

template< typename T, typename IDX >
inline std::enable_if_t< bufferOps::is_container< T >, localIndex >
UnpackDataByIndexDevice( buffer_unit_type const * & buffer, T const & var, IDX & idx, parallelDeviceEvents & events, MPI_Op op )
{ return bufferOps::UnpackDataByIndexDevice( buffer, var, idx, events, op ); }

template< typename T, typename IDX >
inline std::enable_if_t< !bufferOps::is_container< T >, localIndex >
UnpackDataByIndexDevice( buffer_unit_type const * &, T const &, IDX &, parallelDeviceEvents &, MPI_Op )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") by index on device. Operation not supported." );
  return 0;
}

#if defined(GEOS_USE_PYGEOSX)

template< typename T >
inline std::enable_if_t< LvArray::python::CanCreate< T >, PyObject * >
createPythonObject( T & object )
{ return LvArray::python::create( object ); }

template< typename T >
inline std::enable_if_t< !LvArray::python::CanCreate< T >, PyObject * >
createPythonObject( T & )
{ return nullptr; }

#endif

} // namespace wrapperHelpers
} // namespace dataRepository
} // namespace geos

#undef RESTART_TYPE_LOGGING

#endif // GEOS_DATAREPOSITORY_WRAPPERHELPERS_HPP_
