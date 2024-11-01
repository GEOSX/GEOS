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
 * @file SourceFluxStatistics.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SOURCEFLUXSTATISTICS_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SOURCEFLUXSTATISTICS_HPP_

#include "../FieldStatisticsBase.hpp"
#include "FlowSolverBase.hpp"
#include "mesh/DomainPartition.hpp"


namespace geos
{

/**
 * @class SourceFluxStatsAggregator
 *
 * Task class allowing for the computation of aggregate statistics of SourceFluxBoundaryCondition
 */
class SourceFluxStatsAggregator final : public FieldStatisticsBase< FlowSolverBase >
{
public:

  /**
   * @brief Aggregated flux statistics data.
   */
  struct StatData
  {
    /// Amount of fluid produced by the flux(es). Negative if injecting. One value for each fluid phase.
    /// In kg by default, or in mol if useMass = 0 on the solver.
    array1d< real64 > m_producedMass;
    /// Flux(es) production rate. Negative if injecting. One value for each fluid phase.
    /// In kg/s by default, or in mol/s if useMass = 0 on the solver.
    array1d< real64 > m_productionRate;
    /// Number of elements in which we are producing / injecting
    integer m_elementCount = 0;

    /**
     * @brief resize the phase data arrays if needed
     * @param phaseCount the count of phases
     */
    void allocate( integer phaseCount );
    /**
     * @brief reset the stats to 0.0
     */
    void reset();
    /**
     * @return the phase count for phasic stats
     */
    integer getPhaseCount() const
    { return m_producedMass.size(); }
    /**
     * @brief Aggregate the statistics of the instance with those of another one.
     * @param other the other stats structure.
     */
    void combine( StatData const & other );
    /**
     * @brief Aggregate the statistics of the instance with those from all instances from other MPI ranks.
     * Must be called only once per timestep.
     */
    void mpiReduce();
  };
  /**
   * @brief Class that aggregate statistics of a flux over multiple time-steps for a given
   * SourceFluxStatsAggregator and a for a given mesh part (i.e. a subregion, a region...).
   */
  class WrappedStats
  {
public:
    /**
     * @brief Set the subjects targeted by the stats.
     * @param aggregatorName the name of the targeted SourceFluxStatsAggregator
     * @param fluxName the name of the targeted SourceFluxBoundaryCondition
     */
    void setTarget( string_view aggregatorName, string_view fluxName );

    /**
     * @brief Set the current time step stats. Accumulate the statistics only if the time is strictly
     * after the previous time + dt (override the statistics if the current timestep is cut).
     * @param currentTime  time of the timestep start since simulation starting
     * @param dt           time delta of the current timestep
     * @param producedMass Amount of produced fluid (see StatData::m_producedMass).
     * @param elementCount number of cell elements concerned by this instance
     */
    void gatherTimeStepStats( real64 currentTime, real64 dt,
                              arrayView1d< real64 const > const & producedMass,
                              integer elementCount );

    /**
     * @brief Finalize the period statistics of each timestep gathering and render data over all mpi ranks.
     * The result can be read by the data() assessor.
     * @note This method must be synchronously called by all MPI ranks.
     */
    void finalizePeriod();

    /**
     * @return the reference to the wrapped stats data collected over the last period (one timestep or more), computed by finalizePeriod()
     */
    StatData & stats()
    { return m_stats; }

    /**
     * @return the reference to the wrapped stats data collected over the last period (one timestep or more), computed by finalizePeriod()
     */
    StatData const & stats() const
    { return m_stats; }

    /**
     * @return the start time of the wrapped stats period (in s)
     */
    real64 getStatsPeriodStart() const
    { return m_statsPeriodStart; }

    /**
     * @return the duration of the wrapped stats period (in s)
     */
    real64 getStatsPeriodDuration() const
    { return m_statsPeriodDT; }

    /**
     * @return the name of the SourceFluxStatsAggregator that want to collect data on this instance.
     */
    string_view getAggregatorName() const
    { return m_aggregatorName; }

    /**
     * @return the name of the SourceFluxBoundaryCondition from which we are collecting data on this instance.
     */
    string_view getFluxName() const
    { return m_fluxName; }
private:
    /// stats data collected over the last period (one timestep or more), computed by finalizePeriod()
    StatData m_stats;
    /// the start time of the wrapped stats period (in s)
    real64 m_statsPeriodStart;
    /// the duration of the wrapped stats period (in s)
    real64 m_statsPeriodDT;

    /**
     * @brief This struct is used to accumulate statistics along one or more timesteps.
     */
    struct PeriodStats
    {
      /// Fluid production during current time-step. Same unit as StatData::productedMass.
      array1d< real64 > m_timeStepMass;
      /// Fluid production during all previous time-step of the current period. Same unit as StatData::productedMass.
      array1d< real64 > m_periodPendingMass;
      /// time of when the timestep starts (since the simulation start).
      real64 m_timeStepStart = 0.0;
      /// time that the current timestep is simulating.
      real64 m_timeStepDeltaTime = 0.0;
      /// start time of the current period.
      real64 m_periodStart = 0.0;
      /// delta time from all previous time-step of the current period.
      real64 m_periodPendingDeltaTime = 0.0;
      /// number of cell elements targeted by this instance
      integer m_elementCount = 0;
      /// Did the period statistics gathering started ?
      bool m_isGathering = false;

      /**
       * @brief resize the phase data arrays if needed
       * @param phaseCount the count of phases
       */
      void allocate( integer phaseCount );
      /**
       * @brief reset the stats to 0.0
       */
      void reset();
      /**
       * @return the phase count for phasic stats
       */
      integer getPhaseCount() const
      { return m_timeStepMass.size(); }
    } m_periodStats;

    /// Name of the SourceFluxStatsAggregator that want to collect data on this instance.
    string m_aggregatorName;
    /// Name of the SourceFluxBoundaryCondition from which we are collecting data on this instance.
    string m_fluxName;
  };


  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  SourceFluxStatsAggregator( const string & name,
                             Group * const parent );

  /**
   * @return The catalog name
   */
  static string catalogName() { return "SourceFluxStatistics"; }

  /**
   * @copydoc ExecutableGroup::execute()
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) override;

  /**
   * @brief Apply a functor to WrappedStats that combines all stats for each target solver
   *        discretization mesh levels.
   * @param domain   the domain for which we want the statistics
   * @param lambda   the functor that will be called for each WrappedStats. Takes in parameter the MeshLevel
   *                 reference and the reference to the WrappedStats that combines all stats for the instance.
   * @tparam LAMBDA  the type of lambda function to call in the function
   * @note To be retrieved, the WrappedStats structs must be registered on the container during the
   * registerDataOnMesh() call.
   */
  template< typename LAMBDA >
  void forMeshLevelStatsWrapper( DomainPartition & domain, LAMBDA && lambda );
  /**
   * @brief Apply a functor to each WrappedStats that combines the stats over all region for a flux.
   * @param meshLevel the mesh level for which we want the statistics.
   * @param lambda    the functor that will be called for each WrappedStats. Takes in parameter the MeshLevel
   *                  reference and the reference to one of the flux WrappedStats.
   * @tparam LAMBDA   the type of lambda function to call in the function
   * @note To be retrieved, the WrappedStats structs must be registered on the container during the
   * registerDataOnMesh() call.
   */
  template< typename LAMBDA >
  void forAllFluxStatsWrappers( MeshLevel & meshLevel, LAMBDA && lambda );
  /**
   * @brief Apply a functor to all simulated region WrappedStats (of the given MeshLevel) that target a
   *        given flux.
   * @param meshLevel the mesh level we want to loop over all its regions.
   * @param fluxName  the name of the flux from which we want the statistics.
   * @param lambda    the functor that will be called for each WrappedStats. Takes in parameter the
   *                  ElementRegionBase reference and the reference of its WrappedStats.
   * @tparam LAMBDA   the type of lambda function to call in the function
   * @note To be retrieved, the WrappedStats structs must be registered on the container during the
   * registerDataOnMesh() call.
   */
  template< typename LAMBDA >
  void forAllRegionStatsWrappers( MeshLevel & meshLevel, string_view fluxName, LAMBDA && lambda );
  /**
   * @brief Apply a functor to all subregion WrappedStats (of the given region) that target a given flux.
   * @param region    the region from which we want to execute the lambda for each of its sub-region.
   * @param fluxName  the name of the flux from which we want the statistics.
   * @param lambda    the functor that will be called for each WrappedStats. Takes in parameter the
   *                  ElementSubRegionBase reference and the reference of its WrappedStats.
   * @tparam LAMBDA   the type of lambda function to call in the function
   * @note To be retrieved, the WrappedStats structs must be registered on the container during the
   * registerDataOnMesh() call.
   */
  template< typename LAMBDA >
  void forAllSubRegionStatsWrappers( ElementRegionBase & region, string_view fluxName, LAMBDA && lambda );

  /**
   * @brief Apply a functor to all WrappedStats of the given group that target a given flux (and
   * potentially multiple SourceFluxStatsAggregator).
   * The functor takes in parameter the reference to the currently processed WrappedStats.
   * @note To be retrieved, the WrappedStats structs must be registered on the container during the
   * registerDataOnMesh() call.
   * @tparam LAMBDA   the type of lambda function to call in the function
   * @param container the container from which we want the statistics.
   * @param fluxName  the name of the flux from which we want the statistics.
   * @param lambda    the functor that will be called for each WrappedStats. Takes in parameter the
   *                  reference to the currently processed WrappedStats.
   */
  template< typename LAMBDA >
  static void forAllFluxStatWrappers( Group & container, string_view fluxName, LAMBDA && lambda );

  /**
   * @return a string used to name the wrapper that is added to each region that is simulated by the solver.
   * The string is unique within the region for the SourceFluxBoundaryCondition and the SourceFluxStatsAggregator.
   * @param fluxName The name of the flux. For the mesh-level global stats, fluxSetWrapperString() can be used.
   */
  inline string getStatWrapperName( string_view fluxName ) const;


  /**
   * @brief View keys
   */
  struct viewKeyStruct
  {
    /// @return The key for setName
    constexpr inline static string_view fluxNamesString() { return "fluxNames"; }
    /// @return The key for statistics wrapper name that targets all region
    constexpr inline static string_view allRegionWrapperString() { return "all_regions"; }
    /// @return The key for statistics wrapper name that targets all fluxes of the set
    constexpr inline static string_view fluxSetWrapperString() { return "flux_set"; }
  };


private:
  using Base = FieldStatisticsBase< FlowSolverBase >;

  /// the names of the SourceFlux(s) for which we want the statistics
  string_array m_fluxNames;

  /**
   * @copydoc Group::registerDataOnMesh(Group &)
   */
  void registerDataOnMesh( Group & meshBodies ) override;

  /**
   * @copydoc Group::postInputInitialization()
   */
  void postInputInitialization() override;

  dataRepository::Wrapper< WrappedStats > & registerWrappedStats( Group & group,
                                                                  string_view fluxName,
                                                                  string_view elementSetName );

  /**
   * @brief If requested, output in the log and CSV the given statistics.
   * @param minLogLevel    the min log level to output any line.
   * @param elementSetName the region / sub-subregion name concerned by the statistics.
   * @param stats          the statistics that must be output in the log.
   */
  void writeStatsToLog( integer minLogLevel, string_view elementSetName, WrappedStats const & stats );
  /**
   * @brief If CSV is enabled, create or append a new CSV file.
   * @param elementSetName the region / sub-subregion name concerned by the statistics.
   * @param stats          the statistics that must be output in the log.
   * @param writeHeader    If true, create the CSV with the header. If false, append it with the statistics.
   */
  void writeStatsToCSV( string_view elementSetName, WrappedStats const & stats, bool writeHeader );

};


template< typename LAMBDA >
void SourceFluxStatsAggregator::forAllFluxStatWrappers( Group & container,
                                                        string_view fluxName,
                                                        LAMBDA && lambda )
{
  container.forWrappers< WrappedStats >( [&]( dataRepository::Wrapper< WrappedStats > & statsWrapper )
  {
    if( statsWrapper.referenceAsView().getFluxName() == fluxName )
    {
      lambda( statsWrapper.reference() );
    }
  } );
}

template< typename LAMBDA >
void SourceFluxStatsAggregator::forMeshLevelStatsWrapper( DomainPartition & domain,
                                                          LAMBDA && lambda )
{
  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                            [&] ( string const &,
                                                  MeshLevel & meshLevel,
                                                  arrayView1d< string const > const & )
  {
    string const wrapperName = getStatWrapperName( viewKeyStruct::fluxSetWrapperString() );
    WrappedStats & stats = meshLevel.getWrapper< WrappedStats >( wrapperName ).reference();

    lambda( meshLevel, stats );
  } );
}
template< typename LAMBDA >
void SourceFluxStatsAggregator::forAllFluxStatsWrappers( MeshLevel & meshLevel,
                                                         LAMBDA && lambda )
{
  for( string const & fluxName : m_fluxNames )
  {
    string const wrapperName = getStatWrapperName( fluxName );
    WrappedStats & stats = meshLevel.getWrapper< WrappedStats >( wrapperName ).reference();

    lambda( meshLevel, stats );
  }
}
template< typename LAMBDA >
void SourceFluxStatsAggregator::forAllRegionStatsWrappers( MeshLevel & meshLevel,
                                                           string_view fluxName,
                                                           LAMBDA && lambda )
{
  string const wrapperName = getStatWrapperName( fluxName );
  meshLevel.getElemManager().forElementRegions( [&]( ElementRegionBase & region )
  {
    WrappedStats & stats = region.getWrapper< WrappedStats >( wrapperName ).reference();

    lambda( region, stats );
  } );
}
template< typename LAMBDA >
void SourceFluxStatsAggregator::forAllSubRegionStatsWrappers( ElementRegionBase & region,
                                                              string_view fluxName,
                                                              LAMBDA && lambda )
{
  string const wrapperName = getStatWrapperName( fluxName );
  region.forElementSubRegions( [&]( ElementSubRegionBase & subRegion )
  {
    WrappedStats & stats = subRegion.getWrapper< WrappedStats >( wrapperName ).reference();

    lambda( subRegion, stats );
  } );
}

inline string SourceFluxStatsAggregator::getStatWrapperName( string_view fluxName ) const
{ return GEOS_FMT( "{}_region_stats_for_{}", fluxName, getName() ); }


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_SOURCEFLUXSTATISTICS_HPP_ */
