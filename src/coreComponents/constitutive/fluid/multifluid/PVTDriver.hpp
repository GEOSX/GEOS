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
 * @file PVTDriver.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_PVTDRIVER_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_PVTDRIVER_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geos
{

namespace constitutive
{
class MultiFluidBase;
}

/**
 * @class PVTDriver
 *
 * Class to allow for testing PVT behavior without the
 * complexity of setting up a full simulation.
 *
 */
class PVTDriver : public TaskBase
{
public:
  PVTDriver( const string & name,
             Group * const parent );

  static string catalogName() { return "PVTDriver"; }

  void postInputInitialization() override;

  virtual bool execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                        real64 const GEOS_UNUSED_PARAM( dt ),
                        integer const GEOS_UNUSED_PARAM( cycleNumber ),
                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                        DomainPartition & GEOS_UNUSED_PARAM( domain ) ) override;

  /**
   * @brief Run test using loading protocol in table
   * @param i Fluid constitutive model
   * @param table Table with input/output time history
   */
  template< typename FLUID_TYPE >
  void runTest( FLUID_TYPE & fluid, arrayView2d< real64 > const & table );

  /**
   * @brief Ouput table to file for easy plotting
   */
  void outputResults();

  /**
   * @brief Read in a baseline table from file and compare with computed one (for unit testing purposes)
   */
  void compareWithBaseline();

private:

  /**
   * @brief Get the fluid model from the catalog
   */
  constitutive::MultiFluidBase & getFluid();

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * fluidNameString() { return "fluid"; }
    constexpr static char const * pressureFunctionString() { return "pressureControl"; }
    constexpr static char const * temperatureFunctionString() { return "temperatureControl"; }
    constexpr static char const * numStepsString() { return "steps"; }
    constexpr static char const * outputString() { return "output"; }
    constexpr static char const * baselineString() { return "baseline"; }
    constexpr static char const * feedString() { return "feedComposition"; }
    constexpr static char const * outputMassDensityString() { return "outputMassDensity"; }
    constexpr static char const * outputCompressibilityString() { return "outputCompressibility"; }
    constexpr static char const * outputPhaseCompositionString() { return "outputPhaseComposition"; }
  };

  integer m_numSteps;      ///< Number of load steps
  integer m_numPhases;     ///< Number of fluid phases
  integer m_numComponents; ///< Number of fluid components

  string m_fluidName;                   ///< Fluid identifier
  string m_pressureFunctionName;        ///< Time-dependent function controlling pressure
  string m_temperatureFunctionName;     ///< Time-dependent function controlling temperature
  string m_outputFile;                  ///< Output file (optional, no output if not specified)
  integer m_outputMassDensity{0};       ///< Flag to indicate that the mass density of each phase should be output
  integer m_outputCompressibility{0};   ///< Flag to indicate that the total compressibility should be output
  integer m_outputPhaseComposition{0};  ///< Flag to indicate that phase compositions should be output

  array1d< real64 > m_feed;  ///< User specified feed composition
  array2d< real64 > m_table; ///< Table storing time-history of input/output

  Path m_baselineFile; ///< Baseline file (optional, for unit testing of solid models)

  enum columnKeys { TIME, PRES, TEMP }; ///< Enumeration of "input" column keys for readability

};

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_FLUID_PVTDRIVER_HPP_ */
