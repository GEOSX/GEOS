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
 * @file MultiFluidConstants.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_MULTIFLUIDCONSTANTS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_MULTIFLUIDCONSTANTS_HPP_

#include "LvArray/src/Macros.hpp"
#include "common/PhysicsConstants.hpp"

namespace geos
{
namespace constitutive
{

struct MultiFluidConstants
{
  /**
   * @brief Maximum supported number of fluid components (species)
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr integer MAX_NUM_COMPONENTS = 9;

  /**
   * @brief Maximum supported number of fluid phases
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr integer MAX_NUM_PHASES = 4;

  /**
   * @brief Epsilon used in the calculations to check against zero
   */
  static constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;

  /**
   * @brief Max number of SSI iterations
   */
  static constexpr integer maxSSIIterations = 1000;

  /**
   * @brief Max number of Newton iterations
   */
  static constexpr integer maxNewtonIterations = 30;

  /**
   * @brief Tolerance for checking fugacity ratio convergence
   */
  static constexpr real64 fugacityTolerance = 1.0e-8;

  /**
   * @brief Tolerance for checking SSI iterations for convergence
   */
  static constexpr real64 SSITolerance = 1.0e-3;

  /**
   * @brief Tolerance for checking Newton iterations for convergence
   */
  static constexpr real64 newtonTolerance = 1.0e-12;

  /**
   * @brief Tolerance for baseline comparisons
   * @note Used by PVTDriver
   */
  static constexpr real64 baselineTolerance = 1.0e-3;

  /**
   * @brief Minimum saturation or mole fraction for phase or component presence
   */
  static constexpr real64 minForSpeciesPresence = 1.0e-10;

};

} //namespace constitutive

} //namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_MULTIFLUIDCONSTANTS_HPP_
