/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_

#include "common/DataTypes.hpp"
#include "HYPRE_utilities.h"

namespace geosx
{

inline HYPRE_BigInt * toHYPRE_BigInt( globalIndex * const index )
{
  return reinterpret_cast<HYPRE_BigInt*>(index);
}

inline HYPRE_BigInt const * toHYPRE_BigInt( globalIndex const * const index )
{
  return reinterpret_cast<HYPRE_BigInt const*>(index);
}

inline HYPRE_Real * toHYPRE_Real( real64 * const value )
{
  return reinterpret_cast<HYPRE_Real*>(value);
}

inline HYPRE_Real const * toHYPRE_Real( real64 const * const value )
{
  return reinterpret_cast<HYPRE_Real const*>(value);
}

}

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_*/
