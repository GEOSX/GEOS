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



/* UNCRUSTIFY-OFF */

/**
 * @file TimingMacros.hpp
 *
 * A collection of timing-related macros that wrap Caliper.
 */

#ifndef GEOS_COMMON_TIMINGMACROS_HPP_
#define GEOS_COMMON_TIMINGMACROS_HPP_

#include "common/GeosxConfig.hpp"
#include "GeosxMacros.hpp"

namespace timingHelpers
{
  inline std::string stripPF( char const * prettyFunction )
  {
    std::string input(prettyFunction);
    std::string::size_type const end = input.find_first_of( '(' );
    std::string::size_type const beg = input.find_last_of( ' ', end)+1;
    return input.substr( beg, end-beg );
  }
}


#if defined( GEOS_USE_CUDA ) && defined( GEOS_USE_CUDA_NVTOOLSEXT )

#include "nvToolsExt.h"

namespace timingHelpers
{
  enum NVTXColors {
    GREEN = 0xff00ff00,
    BLUE = 0xff0000ff,
    YELLOW = 0xffffff00,
    PURPLE = 0xffff00ff,
    AQUA = 0xff00ffff,
    RED = 0xffff0000,
    WHITE = 0xffffffff
  };

  class NVTXScopeTracer {
  public:
    NVTXScopeTracer( const char* name, NVTXColors color ) {
      nvtxEventAttributes_t eventAttrib = { 0 };
      eventAttrib.version = NVTX_VERSION;
      eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
      eventAttrib.colorType = NVTX_COLOR_ARGB;
      eventAttrib.color = color;
      eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
      eventAttrib.message.ascii = name;
      nvtxRangePushEx(&eventAttrib);
    }
    ~NVTXScopeTracer() {
      nvtxRangePop();
    }
  };

}
/// Mark a scope with NVTX with name and color (only a local helper should not be used elsewhere)
#  define GEOS_NVTX_MARK_SCOPE_COLORED(name, color) timingHelpers::NVTXScopeTracer __FILE__ ## _ ## __LINE__ ## _ ## scopeTracer = timingHelpers::NVTXScopeTracer(name, color)
/// Mark a scope with NVTX with a given name and color purple
#  define GEOS_NVTX_MARK_SCOPE(name) GEOS_NVTX_MARK_SCOPE_COLORED(STRINGIZE_NX(name), timingHelpers::PURPLE)
/// Mark a function with NVTX using function name and color blue
#  define GEOS_NVTX_MARK_FUNCTION GEOS_NVTX_MARK_SCOPE_COLORED(timingHelpers::stripPF(__PRETTY_FUNCTION__).c_str(), timingHelpers::BLUE)
#else
/// @cond DO_NOT_DOCUMENT
#  define GEOS_NVTX_MARK_SCOPE(name)
#  define GEOS_NVTX_MARK_FUNCTION
/// @endcond
#endif /* USE_CUDA */


#ifdef GEOS_USE_CALIPER
#include <caliper/cali.h>
#include <sys/time.h>
#include <string>
#include <iostream>

/// Mark a function or scope for timing with a given name
#define GEOS_CALIPER_MARK_SCOPE(name) cali::Function __cali_ann##__LINE__(STRINGIZE_NX(name))

/// Mark a function for timing using a compiler-provided name
#define GEOS_CALIPER_MARK_FUNCTION cali::Function __cali_ann##__func__(timingHelpers::stripPF(__PRETTY_FUNCTION__).c_str())

/// Mark the beginning of timed statement group
#define GEOS_CALIPER_MARK_BEGIN(name) CALI_MARK_BEGIN(STRINGIZE(name))

/// Mark the end of timed statements group
#define GEOS_CALIPER_MARK_END(name) CALI_MARK_END(STRINGIZE(name))

/// Mark the beginning of function, only useful when you don't want to or can't mark the whole function.
#define GEOS_CALIPER_MARK_FUNCTION_BEGIN CALI_MARK_FUNCTION_BEGIN

/// Mark the end of function, only useful when you don't want to or can't mark the whole function.
#define GEOS_CALIPER_MARK_FUNCTION_END CALI_MARK_FUNCTION_END

#else // GEOS_USE_CALIPER

/// @cond DO_NOT_DOCUMENT
#define GEOS_CALIPER_MARK_SCOPE(name)
#define GEOS_CALIPER_MARK_FUNCTION

#define GEOS_CALIPER_MARK_BEGIN(name)
#define GEOS_CALIPER_MARK_END(name)

#define GEOS_CALIPER_MARK_FUNCTION_BEGIN
#define GEOS_CALIPER_MARK_FUNCTION_END
/// @endcond

#endif // GEOS_USE_CALIPER

/// Mark scope with both Caliper and NVTX if enabled
#define GEOS_MARK_SCOPE(name) GEOS_CALIPER_MARK_SCOPE(name); GEOS_NVTX_MARK_SCOPE(name)
/// Mark function with both Caliper and NVTX if enabled
#define GEOS_MARK_FUNCTION GEOS_CALIPER_MARK_FUNCTION; GEOS_NVTX_MARK_FUNCTION
/// Mark the beginning of function, only useful when you don't want to or can't mark the whole function.
#define GEOS_MARK_FUNCTION_BEGIN(name) GEOS_CALIPER_MARK_FUNCTION_BEGIN(name)
/// Mark the end of function, only useful when you don't want to or can't mark the whole function.
#define GEOS_MARK_FUNCTION_END(name) GEOS_CALIPER_MARK_FUNCTION_END(name)
/// Mark the beginning of timed statement group
#define GEOS_MARK_BEGIN(name) GEOS_CALIPER_MARK_BEGIN(name)
/// Mark the end of timed statements group
#define GEOS_MARK_END(name) GEOS_CALIPER_MARK_END(name)

/// Get current time of day as a floating point number of seconds in a variable @p time.
#ifdef GEOS_USE_TIMERS
#define GEOS_GET_TIME( time )                                                 \
  real64 time;                                                                 \
  do                                                                           \
  {                                                                            \
    timeval tim;                                                               \
    gettimeofday(&tim, nullptr);                                               \
    time = tim.tv_sec + (tim.tv_usec / 1000000.0);                             \
  } while (false)
#else
#define GEOS_GET_TIME( time )
#endif

#endif // GEOS_COMMON_TIMINGMACROS_HPP_
