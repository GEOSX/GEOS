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
 * @file SinglePhaseStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"

namespace geos
{

class SinglePhaseBase;

/**
 * @class SinglePhaseStatistics
 *
 * Task class allowing for the computation of aggregate statistics in single-phase simulations
 */
class SinglePhaseStatistics : public FieldStatisticsBase< SinglePhaseBase >
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  SinglePhaseStatistics( const string & name,
                         Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "SinglePhaseStatistics"; }

  /**
   * @defgroup Tasks Interface Functions
   *
   * This function implements the interface defined by the abstract TaskBase class
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**@}*/


  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the region statistics
    constexpr static char const * regionStatisticsString() { return "regionStatistics"; }
  };

  struct RegionStatistics
  {
    /// average region pressure
    real64 averagePressure;
    /// minimum region pressure
    real64 minPressure;
    /// maximum region pressure
    real64 maxPressure;

    /// minimum region delta pressure
    real64 minDeltaPressure;
    /// maximum region delta pressure
    real64 maxDeltaPressure;

    // fluid mass
    real64 totalMass;

    /// average region temperature
    real64 averageTemperature;
    /// minimum region temperature
    real64 minTemperature;
    /// maximum region temperature
    real64 maxTemperature;

    /// total region pore volume
    real64 totalPoreVolume;
    /// total region uncompacted pore volume
    real64 totalUncompactedPoreVolume;
  };

private:

  using Base = FieldStatisticsBase< SinglePhaseBase >;

  /**
   * @brief Compute some statistics on the reservoir (average field pressure, etc)
   * @param[in] mesh the mesh level object
   * @param[in] regionNames the array of target region names
   */
  void computeRegionStatistics( real64 const time,
                                MeshLevel & mesh,
                                arrayView1d< string const > const & regionNames ) const;


  void registerDataOnMesh( Group & meshBodies ) override;

};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASESTATISTICS_HPP_ */
