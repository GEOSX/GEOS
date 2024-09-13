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

#ifndef GEOS_RELPERMDRIVER_HPP_
#define GEOS_RELPERMDRIVER_HPP_

#include "events/tasks/TaskBase.hpp"
#include "constitutive/relativePermeability/TableRelativePermeabilityHysteresis.hpp"

namespace geos
{

class RelpermDriver : public TaskBase
{

public:
  RelpermDriver( const string & name,
                 Group * const parent );

  static string catalogName()
  { return "RelpermDriver"; }

  void postInputInitialization() override;

  virtual bool execute( real64 const GEOS_UNUSED_PARAM( time_n ),
                        real64 const GEOS_UNUSED_PARAM( dt ),
                        integer const GEOS_UNUSED_PARAM( cycleNumber ),
                        integer const GEOS_UNUSED_PARAM( eventCounter ),
                        real64 const GEOS_UNUSED_PARAM( eventProgress ),
                        DomainPartition &
                        GEOS_UNUSED_PARAM( domain ) ) override;

  /**
   * @brief Run test using loading protocol in table
   * @param i Relperm constitutive model
   * @param table Table with input/output time history
   */
  template< typename RELPERM_TYPE >
  std::enable_if_t< std::is_same< constitutive::TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
  runTest( RELPERM_TYPE & relperm,
           const arrayView2d< real64, 1 > & table );

  template< typename RELPERM_TYPE >
  std::enable_if_t< !std::is_same< constitutive::TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
  runTest( RELPERM_TYPE & relperm,
           const arrayView2d< real64, 1 > & table );

  /**
   * @brief Output table to file for easy plotting
   */
  void outputResults();

  /**
   * @brief Read in a baseline table from file and compare with computed one (for unit testing purposes)
   */
  void compareWithBaseline();

private:

  template< typename RELPERM_TYPE >
  void resizeTables();

  template< typename RELPERM_TYPE >
  std::enable_if_t< std::is_same< constitutive::TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
  resizeTable();

  template< typename RELPERM_TYPE >
  std::enable_if_t< !std::is_same< constitutive::TableRelativePermeabilityHysteresis, RELPERM_TYPE >::value, void >
  resizeTable();

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    constexpr static char const * relpermNameString()
    { return "relperm"; }

    constexpr static char const * numStepsString()
    { return "steps"; }

    constexpr static char const * outputString()
    { return "output"; }

    constexpr static char const * baselineString()
    { return "baseline"; }
  };

  integer m_numSteps;      ///< Number of load steps
  integer m_numColumns;    ///< Number of columns in data table (depends on number of fluid phases)
  integer m_numPhases;     ///< Number of fluid phases

  string m_relpermName;               ///< relPermType identifier
  string m_outputFile;              ///< Output file (optional, no output if not specified)

  array2d< real64 > m_table; ///< Table storing time-history of input/output

  Path m_baselineFile; ///< Baseline file (optional, for unit testing of solid models)

  enum columnKeys
  {
    TIME
  }; ///< Enumeration of "input" column keys for readability

  static constexpr real64 m_baselineTol = 1e-3; ///< Comparison tolerance for baseline results
};


}

#endif //GEOS_RELPERMDRIVER_HPP_
