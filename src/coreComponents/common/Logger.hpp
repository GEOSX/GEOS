/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * Logger.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: settgast1
 */

#ifndef SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_

#include <string>
#include "GeosxConfig.hpp"
#include <sstream>
#ifdef GEOSX_USE_ATK
#include "slic/slic.hpp"
#endif


namespace geosx
{

void geos_abort( std::string message );

#ifdef GEOSX_USE_ATK

#define GEOS_ERROR_IF(EXP, msg) SLIC_ERROR_IF(EXP, msg)

#define GEOS_WARNING_IF(EXP, msg) SLIC_WARNING_IF(EXP, msg)

#define GEOS_ASSERT(EXP, msg) SLIC_ASSERT_MSG(EXP, msg)

#define GEOS_CHECK(EXP, msg) SLIC_CHECK_MSG(EXP, msg)

#define GEOS_INFO_IF(EXP, msg) SLIC_INFO_IF(EXP, msg)

#else /* GEOSX_USE_ATK */

#define GEOS_ERROR_IF(EXP, msg)
  do {                                                        \
    if (EXP)                                                  \
    {                                                         \
      std::cerr<< "***** GEOS_ERROR "<<std::endl;             \
      std::cerr<< "***** FILE: " << __FILE__ << std::endl;    \
      std::cerr<< "***** LINE: " << __LINE__ << std::endl;    \
      std::ostringstream oss;                                 \
      oss << msg;                                             \
      geosx::geos_abort(oss.str());                           \
    }                                                         \
  } while ( false )

#define GEOS_WARNING_IF(EXP, msg)
  do {                                                        \
    if (EXP)                                                  \
    {                                                         \
      std::cerr<< "***** GEOS_WARNING "<<std::endl;           \
      std::cerr<< "***** FILE: " << __FILE__ << std::endl;    \
      std::cerr<< "***** LINE: " << __LINE__ << std::endl;    \
      std::ostringstream oss;                                 \
      oss << msg;                                             \
    }                                                         \
  } while ( false )

#define GEOS_INFO_IF(EXP, msg)
  do {                                                        \
    if (EXP)                                                  \
    {                                                         \
      std::cerr<< "***** GEOS_INFO "<<std::endl;              \
      std::cerr<< "***** FILE: " << __FILE__ << std::endl;    \
      std::cerr<< "***** LINE: " << __LINE__ << std::endl;    \
      std::ostringstream oss;                                 \
      oss << msg;                                             \
    }                                                         \
  } while ( false )

#ifdef GEOSX_DEBUG

#define GEOS_ASSERT(EXP, msg) GEOS_ERROR_IF(EXP, msg)

#define GEOS_CHECK(EXP, msg) GEOS_WARNING_IF(EXP, msg)

#else /* #ifdef GEOSX_DEBUG */

#define GEOS_ASSERT(EXP, msg) ((void) 0)

#define GEOS_CHECK(EXP, msg) ((void) 0)

#endif /* #ifdef GEOSX_DEBUG */

#endif /* #ifdef GEOSX_USE_ATK */

#define GEOS_ERROR(msg) GEOS_ERROR_IF(true, msg)

#define GEOS_WARNING(msg) GEOS_WARNING_IF(true, msg)

#define GEOS_INFO(msg) GEOS_INFO_IF(true, msg)

}

#endif /* SRC_COMPONENTS_CORE_SRC_COMMON_LOGGER_HPP_ */
