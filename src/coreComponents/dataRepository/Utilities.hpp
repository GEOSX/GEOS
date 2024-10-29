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
 * @file Utilities.hpp
 */

#include "common/DataTypes.hpp"

namespace geos
{
namespace dataRepository
{
class Group;

/**
 */
/**
 * @brief Prints the memory allocations for a group in the data repository recursively
 *
 * @param group The group to print
 * @param[in] indent The level of indentation to add to this level of output.
 * @param[in] threshold The allocation size required to output a allocation size.
 */
void printMemoryAllocation( Group const & group, integer const indent, real64 const threshold );

}
}
