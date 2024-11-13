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
 * @file PoromechanicsInitialization.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSINITIALIZATION_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSINITIALIZATION_HPP_

#include "events/tasks/TaskBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsStateReset.hpp"

namespace geos
{
class SolidMechanicsStatistics;

/**
 * @class PoromechanicsInitialization
 *
 * Task class used to perform stress initialization in the poromechanical simulations
 * A common workflow might be:
 * (1) performStressInitialization = true
 * (2) solve an equilibration problem for initial conditions,
 * (3) performStressInitialization = false
 * (4) run a normal simulation
 *
 */
template< typename POROMECHANICS_SOLVER >
class PoromechanicsInitialization : public TaskBase
{
public:

  /**
   * @brief Constructor for the initialization class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  PoromechanicsInitialization( const string & name,
                               Group * const parent );

  /// Destructor for the class
  ~PoromechanicsInitialization() override;

  /// Accessor for the catalog name
  static string catalogName()
  {
    return POROMECHANICS_SOLVER::catalogName() + "Initialization";
  }

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

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct
  {
    /// String for the poromechanics solver name
    constexpr static char const * poromechanicsSolverNameString() { return "poromechanicsSolverName"; }
    /// String for the solid mechanics statistics name
    constexpr static char const * solidMechanicsStatisticsNameString() { return "solidMechanicsStatisticsName"; }
  };

  void postInputInitialization() override;

//  void registerDataOnMesh( Group & meshBodies ) override;

  /// Name of the poromechanics solver
  string m_poromechanicsSolverName;

  /// Name of the solid mechanics statistics
  string m_solidMechanicsStatisticsName;

  /// Pointer to the poromechanics solver
  POROMECHANICS_SOLVER * m_poromechanicsSolver;

  /// Pointer to the solid mechanics statistics
  SolidMechanicsStatistics * m_solidMechanicsStatistics;

  SolidMechanicsStateReset m_solidMechanicsStateResetTask;

};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSINITIALIZATION_HPP_ */
