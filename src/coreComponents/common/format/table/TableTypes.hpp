

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
 * @file TableData.hpp
 */


#ifndef GEOS_COMMON_FORMAT_TABLETYPES_HPP
#define GEOS_COMMON_FORMAT_TABLETYPES_HPP

#include "common/DataTypes.hpp"
#include "common/format/Format.hpp"

namespace geos
{

    

enum class CellType : char  
{
  MERGE = '\x01',
  SEPARATOR = '\x02',
  Header = '\x03',
  Value = '\x04',
};

struct DataType
{
  CellType type;
  string value;
};

}

#endif /* GEOS_COMMON_FORMAT_TABLETYPES_HPP */
