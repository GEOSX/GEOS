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

#ifndef GEOS_DATAREPOSITORY_BUFFEROPSDEVICE_H_
#define GEOS_DATAREPOSITORY_BUFFEROPSDEVICE_H_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/traits.hpp"
#include "LvArray/src/limits.hpp"
#include "BufferOps.hpp"

#include "common/GEOS_RAJA_Interface.hpp"

#include <type_traits>

namespace geos
{

namespace bufferOps
{

namespace internal
{
// Forward decl so we can use this for contained types
template< typename T >
struct is_device_packable_helper;

/// Whether an object of type T is itself packable on device
template< typename T >
constexpr bool is_device_packable_object = std::is_arithmetic< T >::value ||
                                           std::is_enum< T >::value ||
                                           traits::is_tensorT< T >;

/// Whether an object is an lvarray arrayview which contains device-packable values when fully indexed
template< typename >
constexpr bool is_device_packable_array = false;

template< typename T, int NDIM, int USD >
constexpr bool is_device_packable_array< ArrayView< T, NDIM, USD > > = is_device_packable_helper< typename std::remove_const< T >::type >::value;

template< typename T >
struct is_device_packable_helper
{
  static constexpr bool value = is_device_packable_object< T > || is_device_packable_array< T >;
};

}

/// Whether an object is device packable
template< typename T >
constexpr bool is_device_packable_v = internal::is_device_packable_helper< std::remove_const_t< std::remove_pointer_t< T > > >::value;

/// Whether an object can be indexed to pack a subset of the contained values on device
template< typename T >
constexpr bool is_device_packable_by_index_v = internal::is_device_packable_array< typename std::decay< T >::type >;

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
GEOS_HOST_DEVICE
localIndex
PackPointerDevice( buffer_unit_type * & buffer,
                   T const * const GEOS_RESTRICT var,
                   localIndex const length );

//------------------------------------------------------------------------------
template< typename T >
GEOS_HOST_DEVICE
localIndex
UnpackPointerDevice( buffer_unit_type const * & buffer,
                     T * const GEOS_RESTRICT var,
                     localIndex const expectedLength,
                     MPI_Op op );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< is_device_packable_v< T >, localIndex >::type
PackDevice( buffer_unit_type * & buffer,
            ArrayView< T const, NDIM, USD > const & var,
            parallelDeviceEvents & events );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
GEOS_HOST_DEVICE
localIndex
PackDevice( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ),
            T const & GEOS_UNUSED_PARAM( var ),
            parallelDeviceEvents & GEOS_UNUSED_PARAM( events ) )
{
  GEOS_ERROR( "Trying to pack data type (" << LvArray::system::demangleType< T >() << ") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< is_device_packable_by_index_v< ArrayView< T const, NDIM, USD > const & >, localIndex >::type
PackByIndexDevice( buffer_unit_type * & buffer,
                   ArrayView< T const, NDIM, USD > const & var,
                   T_INDICES const & indices,
                   parallelDeviceEvents & events );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_INDICES >
localIndex
PackByIndexDevice( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ),
                   T const & GEOS_UNUSED_PARAM( var ),
                   T_INDICES const & GEOS_UNUSED_PARAM( indices ),
                   parallelDeviceEvents & GEOS_UNUSED_PARAM( events ) )
{
  GEOS_ERROR( "Trying to pack data type (" << LvArray::system::demangleType< T >() << ") on device but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD >
typename std::enable_if< is_device_packable_v< T >, localIndex >::type
UnpackDevice( buffer_unit_type const * & buffer,
              ArrayView< T, NDIM, USD > const & var,
              parallelDeviceEvents & events,
              MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
localIndex
UnpackDevice( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ),
              T & GEOS_UNUSED_PARAM( var ),
              parallelDeviceEvents & GEOS_UNUSED_PARAM( events ),
              MPI_Op GEOS_UNUSED_PARAM( op ) )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< is_device_packable_by_index_v< ArrayView< T, NDIM, USD > const & >, localIndex >::type
UnpackByIndexDevice ( buffer_unit_type const * & buffer,
                      ArrayView< T, NDIM, USD > const & var,
                      T_INDICES const & indices,
                      parallelDeviceEvents & events,
                      MPI_Op op );

//------------------------------------------------------------------------------
template< typename T, typename T_INDICES >
localIndex
UnpackByIndexDevice( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ),
                     T & GEOS_UNUSED_PARAM( var ),
                     T_INDICES const & GEOS_UNUSED_PARAM( indices ),
                     parallelDeviceEvents & GEOS_UNUSED_PARAM( events ),
                     MPI_Op GEOS_UNUSED_PARAM( op ) )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") but type is not packable by index." );
  return 0;
}


//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
GEOS_HOST_DEVICE
localIndex
PackDataPointerDevice( buffer_unit_type * & buffer,
                       T const * const GEOS_RESTRICT var,
                       localIndex const length,
                       parallelDeviceEvents & events );

//------------------------------------------------------------------------------
template< typename T >
GEOS_HOST_DEVICE
localIndex
UnpackDataPointerDevice( buffer_unit_type const * & buffer,
                         T * const GEOS_RESTRICT var,
                         localIndex const expectedLength,
                         parallelDeviceEvents & events,
                         MPI_Op op );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< is_device_packable_v< T >, localIndex >::type
PackDataDevice( buffer_unit_type * & buffer,
                ArrayView< T const, NDIM, USD > const & var,
                parallelDeviceEvents & events );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
GEOS_HOST_DEVICE
localIndex
PackDataDevice( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ),
                T const & GEOS_UNUSED_PARAM( var ),
                parallelDeviceEvents & GEOS_UNUSED_PARAM( events ) )
{
  GEOS_ERROR( "Trying to pack data type (" << LvArray::system::demangleType< T >() << ") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< is_device_packable_by_index_v< ArrayView< T const, NDIM, USD > const & >, localIndex >::type
PackDataByIndexDevice ( buffer_unit_type * & buffer,
                        ArrayView< T const, NDIM, USD > const & var,
                        T_INDICES const & indices,
                        parallelDeviceEvents & events );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_INDICES >
localIndex
PackDataByIndexDevice( buffer_unit_type * & GEOS_UNUSED_PARAM( buffer ),
                       T const & GEOS_UNUSED_PARAM( var ),
                       T_INDICES const & GEOS_UNUSED_PARAM( indices ),
                       parallelDeviceEvents & GEOS_UNUSED_PARAM( events ) )
{
  GEOS_ERROR( "Trying to pack data type (" << LvArray::system::demangleType< T >() << ") on device but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD >
typename std::enable_if< is_device_packable_v< T >, localIndex >::type
UnpackDataDevice( buffer_unit_type const * & buffer,
                  ArrayView< T, NDIM, USD > const & var,
                  parallelDeviceEvents & events,
                  MPI_Op op );

//------------------------------------------------------------------------------
template< typename T >
localIndex
UnpackDataDevice( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ),
                  T & GEOS_UNUSED_PARAM( var ),
                  parallelDeviceEvents & GEOS_UNUSED_PARAM( events ),
                  MPI_Op GEOS_UNUSED_PARAM( op ) )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< is_device_packable_by_index_v< ArrayView< T, NDIM, USD > const & >, localIndex >::type
UnpackDataByIndexDevice ( buffer_unit_type const * & buffer,
                          ArrayView< T, NDIM, USD > const & var,
                          T_INDICES const & indices,
                          parallelDeviceEvents & events,
                          MPI_Op op );

//------------------------------------------------------------------------------
template< typename T, typename T_INDICES >
localIndex
UnpackDataByIndexDevice( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ),
                         T & GEOS_UNUSED_PARAM( var ),
                         T_INDICES const & GEOS_UNUSED_PARAM( indices ),
                         parallelDeviceEvents & GEOS_UNUSED_PARAM( events ),
                         MPI_Op GEOS_UNUSED_PARAM( op ) )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") but type is not packable by index." );
  return 0;
}

} // namespace bufferOps
} // namespace geos

#endif // GEOS_DATAREPOSITORY_BUFFEROPSDEVICE_H_
