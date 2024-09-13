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
 * @file layouts.hpp
 */

#ifndef GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_LAYOUTS_HPP
#define GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_LAYOUTS_HPP

#include "common/GeosxConfig.hpp"

#include "LvArray/src/typeManipulation.hpp"
#include "RAJA/RAJA.hpp"

namespace geos
{
namespace constitutive
{
namespace relperm
{

#if defined( GEOS_USE_DEVICE )
/// Constitutive model phase array layout
using LAYOUT_PHASE = RAJA::PERM_JKI;
/// Constitutive model phase relative permeability array layout
using LAYOUT_RELPERM = RAJA::PERM_JKLI;
/// Constitutive model phase relative permeability saturation derivative array layout
using LAYOUT_RELPERM_DS = RAJA::PERM_JKLIM;//FIX ME don't know really

#else

/// Constitutive model phase array layout
using LAYOUT_PHASE = RAJA::PERM_IJK;
using LAYOUT_MOB = RAJA::PERM_IJK;
/// Constitutive model phase relative permeability array layout
using LAYOUT_RELPERM = RAJA::PERM_IJKL;
using LAYOUT_MOB_DC = RAJA::PERM_IJKL;
/// Constitutive model phase relative permeability saturation derivative array layout
using LAYOUT_RELPERM_DS = RAJA::PERM_IJKLM;


#endif

/// Constitutive model phase relative permeability unit stride dimension
static constexpr int USD_PHASE = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE {} );
static constexpr int USD_MOB = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_MOB {} );
/// Constitutive model phase relative permeability unit stride dimension
static constexpr int USD_RELPERM = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_RELPERM{} );
static constexpr int USD_MOB_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_MOB_DC{} );
/// Constitutive model phase relative permeability saturation derivative unit stride dimension
static constexpr int USD_RELPERM_DS = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_RELPERM_DS{} );

} // namespace relperm
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_LAYOUTS_HPP
