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
 * @file StencilDataCollection.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_CELLTOCELLDATACOLLECTOR_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_CELLTOCELLDATACOLLECTOR_HPP_

#include "events/tasks/TaskBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geos
{

/**
 * @class StencilDataCollection
 * Task class allowing the ouput of communicating elements data.
 * Only current usage is to output the cell-cell transmissibility (= "CellToCellDataCollection"), but we could add more
 * implementations to output the transmissibility of boundary or surface connections (or everything that is managed
 * by the stencils).
 */
class StencilDataCollection : public TaskBase
{
public:

  /**
   * @brief Constructor for the statistics class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  StencilDataCollection( const string & name,
                         dataRepository::Group * const parent );

  /// Accessor for the catalog name
  static string catalogName() { return "CellToCellDataCollection"; }

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

  struct viewKeyStruct
  {
    static constexpr char const * solverNameString() { return "flowSolverName"; }
    static constexpr char const * meshNameString() { return "meshBody"; }
    static constexpr char const * cellAGlobalIdString() { return "cellAGlobalId"; }
    static constexpr char const * cellBGlobalIdString() { return "cellBGlobalId"; }
    static constexpr char const * transmissibilityABString() { return "transmissibilityAB"; }
    static constexpr char const * transmissibilityBAString() { return "transmissibilityBA"; }
  };

  class Kernel;

  /**
   * @brief Element-element connection data extracted from the kernel.
   */
  struct KernelConnectionData
  {
    real64 m_transmissibility[2];
    localIndex m_regionId[2];
    localIndex m_subRegionId[2];
    localIndex m_elementId[2];
  };

private:

  using Base = TaskBase;
  using LocalToGlobalMap = ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >;

  /**
   * @brief Element-element connection data to be output by the class.
   */
  struct ConnectionData
  {
    real64 m_transmissibility[2];
    globalIndex m_globalId[2];


    /**
     * @return Construct a ConnectionData from the given kernel data.
     * @param kernelData the kernel extracted data to be converted.
     * @param localToGlobalMap the local to global id map, useful to convert the kernel extracted ids of the elements.
     */
    static ConnectionData fromKernel( KernelConnectionData const & kernelData,
                                      LocalToGlobalMap const & localToGlobalMap );

    /**
     * @brief minus operator to be able to sort instances before outputing them.
     * The sorting is done by the first global id, then by the second global id.
     * @param other compared instance
     * @return true if this instance is to sort before the other instance.
     */
    bool operator<( ConnectionData const & other ) const
    {
      return m_globalId[0] != other.m_globalId[0] ?
             m_globalId[0] < other.m_globalId[0] :
             m_globalId[1] < other.m_globalId[1];
    }
  };

  array1d< globalIndex > m_cellAGlobalId;
  array1d< globalIndex > m_cellBGlobalId;
  array1d< real64 > m_transmissibilityAB;
  array1d< real64 > m_transmissibilityBA;

  /// Name of the target mesh body
  string m_meshName;
  /// Name of the target discretization
  string m_solverName;

  /// Pointer to the target flow solver
  FlowSolverBase const * m_solver = nullptr;
  /// Pointer to the target mesh body
  MeshLevel const * m_meshLevel = nullptr;
  /// Pointer to the target discretization
  FluxApproximationBase const * m_discretization = nullptr;



  void postInputInitialization() override;

  /**
   * @brief Initialization of the internal buffers, must happen:
   *  - after FluxApproximationBase::computeCellStencil()
   *  - before TimeHistoryOutput::initializePostInitialConditionsPostSubGroups()
   */
  void initializePostInitialConditionsPostSubGroups() override;

  /**
   * @brief gather the element-element connection data of the current timestep using a given StencilWrapper.
   * @tparam STENCILWRAPPER_T the type of the StencilWrapper
   * @param mesh the mesh for which we want the data
   * @param stencilWrapper the StencilWrapper to use to compute the element-element connection data
   * @return Return the gathered data in an LvArray
   */
  template< typename STENCILWRAPPER_T >
  array1d< KernelConnectionData > gatherConnectionData( STENCILWRAPPER_T const & stencilWrapper ) const;

  /**
   * @brief Output the element-element connection data of the current timestep.
   * @param mesh the specific mesh for which we output the data. We will also need it to convert the ids to global ids.
   * @param stencil the specific mesh for which we output the data.
   * @param kernelData the connection data, gathered by a kernel.
   */
  void storeConnectionData( string_view stencilName,
                            arrayView1d< KernelConnectionData > const & kernelData );

  void logStoredConnections( string_view stencilName );
};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_FLUIDFLOW_CELLTOCELLDATACOLLECTOR_HPP_ */
