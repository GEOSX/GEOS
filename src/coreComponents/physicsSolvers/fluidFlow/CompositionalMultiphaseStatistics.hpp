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
  struct RegionStatistics : public dataRepository::Group
  {

    /**
     * @brief Constructor for the RegionStatistics struct
     * @param[in] name the name of the task coming from the xml
     * @param[in] parent the parent group of the task
     */
    RegionStatistics( const string & name,
                      Group * const parent );

    struct viewKeyStruct
    {
      constexpr static char const * averagePressureString() { return "averagePressure"; }
      constexpr static char const * minPressureString() { return "minPressure"; }
      constexpr static char const * maxPressureString() { return "maxPressure"; }

      constexpr static char const * minDeltaPressureString() { return "minDeltaPressure"; }
      constexpr static char const * maxDeltaPressureString() { return "maxDeltaPressure"; }

      constexpr static char const * averageTemperatureString() { return "averageTemperature"; }
      constexpr static char const * minTemperatureString() { return "minTemperature"; }
      constexpr static char const * maxTemperatureString() { return "maxTemperature"; }

      constexpr static char const * totalPoreVolumeString() { return "totalPoreVolume"; }
      constexpr static char const * totalUncompactedPoreVolumeString() { return "totalUncompactedPoreVolume"; }

      constexpr static char const * phasePoreVolumeString() { return "phasePoreVolume"; }
      constexpr static char const * phaseMassString() { return "phaseMass"; }
      constexpr static char const * trappedPhaseMassString() { return "trappedPhaseMass"; }
      constexpr static char const * immobilePhaseMassString() { return "immobilePhaseMass"; }
      constexpr static char const * dissolvedComponentMassString() { return "dissolvedComponentMass"; }
    };

    /**
     * @brief Intialiaze region statistics
     * @param numPhases The number of fluid phases
     * @param numComps The number of components
     */
    void init( integer const numPhases, integer const numComps );

    RegionStatistics() = delete;
    /// average region pressure
    real64 m_averagePressure;
    /// minimum region pressure
    real64 m_minPressure;
    /// maximum region pressure
    real64 m_maxPressure;

    /// minimum region delta pressure
    real64 m_minDeltaPressure;
    /// maximum region delta pressure
    real64 m_maxDeltaPressure;

    /// average region temperature
    real64 m_averageTemperature;
    /// minimum region temperature
    real64 m_minTemperature;
    /// maximum region temperature
    real64 m_maxTemperature;

    /// total region pore volume
    real64 m_totalPoreVolume;
    /// total region uncompacted pore volume
    real64 m_totalUncompactedPoreVolume;
    /// phase region phase pore volume
    array1d< real64 > m_phasePoreVolume;

    /// region phase mass (trapped and non-trapped, immobile and mobile)
    array1d< real64 > m_phaseMass;
    /// trapped region phase mass
    array1d< real64 > m_trappedPhaseMass;
    /// immobile region phase mass
    array1d< real64 > m_immobilePhaseMass;
    /// dissolved region component mass
    array2d< real64 > m_dissolvedComponentMass;
  };

  using statsVKS = RegionStatistics::viewKeyStruct;

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
    constexpr static char const * relpermThresholdString() { return "relpermThreshold";}
  };


  /**
   * @brief Compute some statistics on the reservoir (average field pressure, etc)
   * @param[in] time current time
   * @param[in] mesh the mesh level object
   * @param[in] regionNames the array of target region names
   */
  void computeRegionStatistics( real64 const time,
                                MeshLevel & mesh,
                                arrayView1d< string const > const & regionNames );

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

using RegionCompStatsClass = CompositionalMultiphaseStatistics::RegionStatistics;

} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASESTATISTICS_HPP_ */
