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
 * @file Timer.hpp
 */

#ifndef GEOS_COMMON_TIMER_HPP
#define GEOS_COMMON_TIMER_HPP

#include <chrono>

namespace geos
{

/**
 * @class Timer
 * @brief Object that times the duration of its existence.
 */
class Timer
{
public:

  /**
   * @brief Constructor. The time the object is alive is added to @p duration.
   * @param duration A reference to the duration to add to.
   */
  Timer( std::chrono::system_clock::duration & duration ):
    m_start( std::chrono::system_clock::now() ),
    m_duration( duration )
  {}

  /// Destructor. Adds to the referenced duration.
  ~Timer()
  { m_duration += std::chrono::system_clock::now() - m_start; }

private:
  /// The time at which this object was constructed.
  std::chrono::system_clock::time_point const m_start;
  /// A reference to the duration to add to.
  std::chrono::system_clock::duration & m_duration;
};

}

#endif // GEOS_COMMON_TIMER_HPP
