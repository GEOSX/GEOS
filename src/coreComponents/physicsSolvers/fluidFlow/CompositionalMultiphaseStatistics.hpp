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
 * @file CompositionalMultiphaseStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICS_HPP_

#include "physicsSolvers/FieldStatisticsBase.hpp"

namespace geos
{

class CompositionalMultiphaseBase;

/**
 * @class CompositionalMultiphaseStatistics
 *
 * Task class allowing for the computation of aggregate statistics in compositional multiphase simulations
 */
class CompositionalMultiphaseStatistics : public FieldStatisticsBase< CompositionalMultiphaseBase >
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  CompositionalMultiphaseStatistics( const string & name,
                                     Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "CompositionalMultiphaseStatistics"; }

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

private:

  using Base = FieldStatisticsBase< CompositionalMultiphaseBase >;

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the flag deciding the computation of the CFL numbers
    constexpr static char const * computeCFLNumbersString() { return "computeCFLNumbers"; }
    /// String for the flag deciding the computation of the region statistics
    constexpr static char const * computeRegionStatisticsString() { return "computeRegionStatistics"; }
    /// String for the region statistics
    constexpr static char const * regionStatisticsString() { return "regionStatistics"; }
    /// String for the relperm threshold
    constexpr static char const * relpermThresholdString() { return "relpermThreshold"; }
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
    /// phase region phase pore volume
    array1d< real64 > phasePoreVolume;

    /// region phase mass (trapped and non-trapped, immobile and mobile)
    array1d< real64 > phaseMass;
    /// trapped region phase mass
    array1d< real64 > trappedPhaseMass;
    /// immobile region phase mass
    array1d< real64 > immobilePhaseMass;
    /// region component mass
    array2d< real64 > componentMass;


  };

  /**
   * @brief Compute some statistics on the reservoir (average field pressure, etc)
   * @param[in] time current time
   * @param[in] mesh the mesh level object
   * @param[in] regionNames the array of target region names
   */
  void computeRegionStatistics( real64 const time,
                                MeshLevel & mesh,
                                arrayView1d< string const > const & regionNames ) const;

  /**
   * @brief Compute CFL numbers
   * @param[in] time current time
   * @param[in] dt the time step size
   * @param[in] domain the domain partition
   */
  void computeCFLNumbers( real64 const time,
                          real64 const dt,
                          DomainPartition & domain ) const;

  void postInputInitialization() override;

  void registerDataOnMesh( Group & meshBodies ) override;

  /// Flag to decide whether CFL numbers are computed or not
  integer m_computeCFLNumbers;

  /// Flag to decide whether region statistics are computed or not
  integer m_computeRegionStatistics;

  /// Threshold to decide whether a phase is considered "mobile" or not
  real64 m_relpermThreshold;

};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICS_HPP_ */
