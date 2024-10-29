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
                     localIndex const expectedLength );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< can_memcpy< T >, localIndex >::type
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
typename std::enable_if< can_memcpy< T >, localIndex >::type
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
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDevice( buffer_unit_type const * & buffer,
              ArrayView< T, NDIM, USD > const & var,
              parallelDeviceEvents & events );

//------------------------------------------------------------------------------
template< typename T >
localIndex
UnpackDevice( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ),
              T & GEOS_UNUSED_PARAM( var ),
              parallelDeviceEvents & GEOS_UNUSED_PARAM( events ) )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackByIndexDevice ( buffer_unit_type const * & buffer,
                      ArrayView< T, NDIM, USD > const & var,
                      T_INDICES const & indices,
                      parallelDeviceEvents & events,
                      MPI_Op op=MPI_REPLACE );

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
                         parallelDeviceEvents & events );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< can_memcpy< T >, localIndex >::type
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
typename std::enable_if< can_memcpy< T >, localIndex >::type
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
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDataDevice( buffer_unit_type const * & buffer,
                  ArrayView< T, NDIM, USD > const & var,
                  parallelDeviceEvents & events );

//------------------------------------------------------------------------------
template< typename T >
localIndex
UnpackDataDevice( buffer_unit_type const * & GEOS_UNUSED_PARAM( buffer ),
                  T & GEOS_UNUSED_PARAM( var ),
                  parallelDeviceEvents & GEOS_UNUSED_PARAM( events ) )
{
  GEOS_ERROR( "Trying to unpack data type (" << LvArray::system::demangleType< T >() << ") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDataByIndexDevice ( buffer_unit_type const * & buffer,
                          ArrayView< T, NDIM, USD > const & var,
                          T_INDICES const & indices,
                          parallelDeviceEvents & events,
                          MPI_Op op=MPI_REPLACE );

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
