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

/**
 * @file KrylovUtils.hpp
 */
#ifndef GEOS_LINEARALGEBRA_SOLVERS_KRYLOVUTILS_HPP_
#define GEOS_LINEARALGEBRA_SOLVERS_KRYLOVUTILS_HPP_

#include "codingUtilities/Utilities.hpp"

/**
 * @brief Exit solver iteration and report a breakdown if value too close to zero.
 * @param VAR the variable or expression
 */
#define GEOS_KRYLOV_BREAKDOWN_IF_ZERO( VAR ) \
  if( isZero( VAR, 0.0 ) )                  \
  {                                         \
    if( m_params.logLevel >= 1 )            \
    {                                       \
      GEOS_LOG_RANK_0( "Breakdown in " << methodName() << ": " << #VAR << " = " << VAR ); \
    }                                       \
    m_result.status = LinearSolverResult::Status::Breakdown; \
    break;                                  \
  }                                         \

#endif //GEOS_LINEARALGEBRA_SOLVERS_KRYLOVUTILS_HPP_
