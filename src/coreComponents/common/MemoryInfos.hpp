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

#ifndef GEOS_COMMON_MemoryInfos_HPP_
#define GEOS_COMMON_MemoryInfos_HPP_

#include "umpire/util/MemoryResourceTraits.hpp"
#include "common/Logger.hpp"
#include <unistd.h>
#include <iostream>
#if defined( GEOS_USE_CUDA )
#include <cuda.h>
#endif

namespace geos
{

/**
 * @class MemoryInfos
 * @brief Class to fetch and store memory information for different resource types.
 */
class MemoryInfos
{
public:

  /**
   * @brief Constructor for MemoryInfos.
   * @param resourceType The type of memory resource.
   */
  MemoryInfos( umpire::MemoryResourceTraits::resource_type resourceType );

  /**
   * @brief Get the total memory available for the resource type.
   * @return Total memory in bytes.
   */
  size_t getTotalMemory() const;

  /**
   * @brief Get the available memory for the resource type.
   * @return Available memory in bytes.
   */
  size_t getAvailableMemory() const;

  /**
   * @brief Check if physical memory is handled.
   * @return True if physical memory is handled, false otherwise.
   */
  bool isPhysicalMemoryHandled() const;
private:

  ///total memory available.
  size_t m_totalMemory;
  ///Available memory.
  size_t m_availableMemory;
  ///Flag indicating if physical memory is handled.
  bool m_physicalMemoryHandled;
};

}

#endif
