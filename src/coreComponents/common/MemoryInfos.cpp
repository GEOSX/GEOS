/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-LiCense-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "MemoryInfos.hpp"

namespace geos
{
MemoryInfos::MemoryInfos( umpire::MemoryResourceTraits::resource_type resourceType ):
  m_totalMemory( 0 ),
  m_availableMemory( 0 ),
  m_physicalMemoryHandled( 1 )
{
  switch( resourceType )
  {
    case umpire::MemoryResourceTraits::resource_type::host:
    case umpire::MemoryResourceTraits::resource_type::pinned:
      #if defined( _SC_PHYS_PAGES ) && defined( _SC_PAGESIZE )
      m_totalMemory = sysconf( _SC_PHYS_PAGES ) * sysconf( _SC_PAGESIZE );
      #if defined(_SC_AVPHYS_PAGES)
      m_availableMemory = sysconf( _SC_AVPHYS_PAGES ) * sysconf( _SC_PAGESIZE );
      #else
      GEOS_WARNING( "Unknown device avaialable memory size getter for this system." );
      m_availableMemory = 0;
      #endif
      #else
      GEOS_WARNING( "Unknown device physical memory size getter for this compiler." );
      m_physicalMemoryHandled = 0;
      #endif
      break;
    case umpire::MemoryResourceTraits::resource_type::device:
    case umpire::MemoryResourceTraits::resource_type::device_const:
    case umpire::MemoryResourceTraits::resource_type::um:
      #if defined( GEOS_USE_CUDA )
      cudaMemGetInfo( &m_availableMemory, &m_totalMemory );
      #else
      GEOS_WARNING( "Unknown device physical memory size getter for this compiler." );
      m_physicalMemoryHandled = 0;
      #endif
      break;
    default:
      GEOS_WARNING( "Physical memory lookup not implemented" );
      m_physicalMemoryHandled = 0;
      break;
  }
}

size_t MemoryInfos::getTotalMemory() const
{
  return m_totalMemory;
}

size_t MemoryInfos::getAvailableMemory() const
{
  return m_availableMemory;
}

bool MemoryInfos::isPhysicalMemoryHandled() const
{
  return m_physicalMemoryHandled;
}

}
