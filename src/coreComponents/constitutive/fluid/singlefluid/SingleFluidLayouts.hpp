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

/// Constitutive model phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_JKLI;
/// Constitutive model phase composition compositional derivative array layout
using LAYOUT_PHASE_COMP_DC = RAJA::PERM_JKLMI;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID = RAJA::PERM_JI;
/// Constitutive model fluid property  derivative array layout
using LAYOUT_FLUID_DER = RAJA::PERM_JKI;

#else

/// Constitutive model fluid property array layout
using LAYOUT_FLUID = RAJA::PERM_IJ;
/// Constitutive model fluid property  derivative array layout
using LAYOUT_FLUID_DER = RAJA::PERM_IJK;

#endif

/// Constitutive model fluid property unit stride dimension
static constexpr int USD_FLUID = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID{} );
/// Constitutive model fluid property compositional derivative unit stride dimension
static constexpr int USD_FLUID_DER = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID_DER{} );

} // namespace singefluid
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_SINGLEFLUID_LAYOUTS_HPP
