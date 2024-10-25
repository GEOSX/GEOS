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
 * @file CO2EOSSolver.hpp
 */

#include "common/DataTypes.hpp"

#ifndef GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2EOSSOLVER_HPP
#define GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2EOSSOLVER_HPP

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

/**
 * @class CO2EOSSolver
 *
 * A class to solve a nonlinear CO2 EOS equation for a given pair (pres, temp)
 * This class is used to:
 * 1- Solve the CO2 equation of state (to get the CO2 solubility in brine) in CO2Solubility.cpp
 * 2- Solve the Helmholtz energy equation (to get the CO2 density) in SpanWagnerCO2Density.cpp
 */
class CO2EOSSolver
{
public:

  /**
   * @brief Find the solution to f(temp,pres,var) = 0, and return var
   * @param[in] name the catalog name of the calling class
   * @param[in] maxNumNewtonIter the maximum number of Newton iterations
   * @param[in] maxNumBacktrackIter the maximum number of backtracking iterations
   * @param[in] tolerance the tolerance used in the convergence check
   * @param[in] minAbsDeriv the minimum value of the derivative that we allow (since we use FD, the deriv may be zero...)
   * @param[in] maxAbsUpdate the maximum value of the Newton update that we allow
   * @param[in] allowedMinValue the minimum value that we allow (otherwise, we chop back)
   * @param[in] initialGuess the initial guess for Newton's methode
   * @param[in] temp temperature
   * @param[in] pres pressure
   * @param[in] presMultiplierForReporting pressure multiplier used to report problems (since the CO2Solubility solve does use Pascal)
   * @param[in] f the objective function
   *
   * @detail var is reduced volume in the CO2 equation of state, and density in the Helmholtz energy equation
   */
  static real64
  solve( string const & name,
         integer const maxNumNewtonIter,
         integer const maxNumBacktrackIter,
         real64 const & tolerance,
         real64 const & minAbsDeriv,
         real64 const & maxAbsUpdate,
         real64 const & allowedMinValue,
         real64 const & initialGuess,
         real64 const & temp,
         real64 const & pres,
         real64 const & presMultiplierForReporting,
         real64 (* f)( real64 const & t, real64 const & p, real64 const & var ));

};


} // namespace PVTProps

} // namespace constitutive

} // namespace geos


#endif
