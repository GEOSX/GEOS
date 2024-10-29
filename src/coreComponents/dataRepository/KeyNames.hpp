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
 * @file KeyNames.hpp
 */

#ifndef GEOS_DATAREPOSITORY__KEYNAMES_HPP_
#define GEOS_DATAREPOSITORY__KEYNAMES_HPP_

#include <string>

namespace geos
{
namespace dataRepository
{
namespace keys
{

/// @cond DO_NOT_DOCUMENT

static constexpr auto ProblemManager = "Problem";
static constexpr auto cellManager = "cellManager";
static constexpr auto particleManager = "particleManager";

/// @endcond

}
}
}
#endif /* GEOS_DATAREPOSITORY__KEYNAMES_HPP_ */
