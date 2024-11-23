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
 * @file Layouts.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_LAYOUTS_HPP
#define GEOS_CONSTITUTIVE_FLUID_SINGLEFLUID_LAYOUTS_HPP

#include "common/DataTypes.hpp"
#include "common/GeosxConfig.hpp"

#include "LvArray/src/typeManipulation.hpp"
#include "RAJA/RAJA.hpp"

namespace geos
{
namespace constitutive
{

namespace singlefluid
{
struct DerivativeOffset
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// index of derivative wrt temperature
  static integer constexpr dT = 1;

};


/// indices of pressure, temperature
template< integer IS_THERMAL >
struct DerivativeOffsetC {};

template<>
struct DerivativeOffsetC< 1 >
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// index of derivative wrt temperature
  static integer constexpr dT = dP + 1;
  /// number of derivatives
  static integer constexpr nDer =  2;
};
template<>
struct DerivativeOffsetC< 0 >
{
  /// index of derivative wrt pressure
  static integer constexpr dP = 0;
  /// number of derivatives
  static integer constexpr nDer = 1;
};

#if defined( GEOS_USE_DEVICE )

/// Constitutive model phase property array layout
using LAYOUT_PHASE = RAJA::PERM_JKI;
/// Constitutive model phase property compositional derivative array layout
using LAYOUT_PHASE_DC = RAJA::PERM_JKLI;

/// Constitutive model phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_JKLI;
/// Constitutive model phase composition compositional derivative array layout
using LAYOUT_PHASE_COMP_DC = RAJA::PERM_JKLMI;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID = RAJA::PERM_JI;
/// Constitutive model fluid property compositional derivative array layout
using LAYOUT_FLUID_DC = RAJA::PERM_JKI;

///Constitutive model singe fluid property derivative array layout
using LAYOUT_SINGLEFLUID_DC = RAJA::PERMJI;

#else

/// Constitutive model single phase  property array layout with derivatives
/// Constitutive model phase property array layout
using LAYOUT_PHASE = RAJA::PERM_IJK;
/// Constitutive model phase property compositional derivative array layout
using LAYOUT_PHASE_DC = RAJA::PERM_IJKL;

/// Constitutive model phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_IJKL;
/// Constitutive model phase composition compositional derivative array layout
using LAYOUT_PHASE_COMP_DC = RAJA::PERM_IJKLM;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID = RAJA::PERM_IJ;
/// Constitutive model fluid property compositional derivative array layout
using LAYOUT_FLUID_DC = RAJA::PERM_IJK;

/// Constitutive model singe fluid property derivative array layout
using LAYOUT_SINGLEFLUID_DC = RAJA::PERM_IJ;

#endif

/// Constitutive model phase property unit stride dimension
static constexpr int USD_PHASE = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE{} );
/// Constitutive model phase property compositional derivative unit stride dimension
static constexpr int USD_PHASE_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_DC{} );

/// Constitutive model phase composition unit stride dimension
static constexpr int USD_PHASE_COMP = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_COMP{} );
/// Constitutive model phase composition compositional derivative unit stride dimension
static constexpr int USD_PHASE_COMP_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_COMP_DC{} );

/// Constitutive model fluid property unit stride dimension
static constexpr int USD_FLUID = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID{} );
/// Constitutive model fluid property compositional derivative unit stride dimension
static constexpr int USD_FLUID_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID_DC{} );

} // namespace singefluid
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_SINGLEFLUID_LAYOUTS_HPP
