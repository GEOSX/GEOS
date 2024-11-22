/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file WellTags.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLTAGS_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLTAGS_HPP


namespace geos
{

namespace wellTags
{

static constexpr real64 minDensForDivision = 1e-10;

// tag to access well and reservoir elements in perforation rates computation
struct SubRegionTag
{
  static constexpr integer RES  = 0;
  static constexpr integer WELL = 1;
};

// tag to access the next and current well elements of a connection
struct ElemTag
{
  static constexpr integer CURRENT = 0;
  static constexpr integer NEXT    = 1;
};

// define the column offset of the derivatives
struct ColOffset
{
  static constexpr integer DPRES = 0;
  static constexpr integer DCOMP = 1;
};

template< integer NC, integer IS_THERMAL >
struct ColOffset_WellJac;

template< integer NC >
struct ColOffset_WellJac< NC, 0 >
{
  static constexpr integer dP = 0;
  static constexpr integer dC = 1;
  static constexpr integer dQ = dC + NC;
  static integer constexpr nDer =  dQ + 1;

};

template< integer NC >
struct ColOffset_WellJac< NC, 1 >
{
  static constexpr integer dP = 0;
  static constexpr integer dC = 1;
  static constexpr integer dQ = dC + NC;
  static constexpr integer dT = dQ+1;
  /// number of derivatives
  static integer constexpr nDer =  dT + 1;
};

// define the row offset of the residual equations
struct RowOffset
{
  static constexpr integer CONTROL = 0;
  static constexpr integer MASSBAL = 1;
};

template< integer NC, integer IS_THERMAL >
struct RowOffset_WellJac;

template< integer NC >
struct RowOffset_WellJac< NC, 0 >
{
  static constexpr integer CONTROL   = 0;
  static constexpr integer MASSBAL   = 1;
  static constexpr integer VOLBAL    = MASSBAL + NC;
  static constexpr integer nEqn      = VOLBAL+1;
};

template< integer NC >
struct RowOffset_WellJac< NC, 1 >
{
  static constexpr integer CONTROL   = 0;
  static constexpr integer MASSBAL   = 1;
  static constexpr integer VOLBAL    = MASSBAL + NC;
  static constexpr integer ENERGYBAL = VOLBAL+1;
  static constexpr integer nEqn      = ENERGYBAL+1;

};

} // end namespace wellTags

} // end namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLTAGS_HPP
