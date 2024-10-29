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
 * @file SolidMechanicsStateReset.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSTATERESET_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSTATERESET_HPP_

#include "events/tasks/TaskBase.hpp"

namespace geos
{

class SolidMechanicsLagrangianFEM;

/**
 * @class SolidMechanicsStateReset
 *
 * Task class allowing for various "reset" actions on state variables, useful for pre-
 * and post-equilibration actions. Currently, this class gives the user the option
 * to reset displacements to zero and enable/disable inelastic behavior. A common
 * workflow might be:
 * (1) resetDisplacments = false, disableInelasticity = true,
 * (2) solve an equilibration problem for initial conditions,
 * (3) resetDisplacements = true, disablePlasticity = false,
 * (4) run a normal simulation
 *
 */
class SolidMechanicsStateReset : public TaskBase
{
  template< typename > friend class PoromechanicsInitialization;

public:

  /**
   * @brief Constructor for the state reset class
   * @param[in] name the name of the task coming from the xml
   * @param[in] parent the parent group of the task
   */
  SolidMechanicsStateReset( const string & name,
                            Group * const parent );

  /// Destructor for the class
  ~SolidMechanicsStateReset() override;

  /// Accessor for the catalog name
  static string catalogName() { return "SolidMechanicsStateReset"; }

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
    /// String for the solid solver name
    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }
    /// String for the displacement reset flag
    constexpr static char const * resetDisplacementsString() { return "resetDisplacements"; }
    /// String for the inelasticity flag
    constexpr static char const * disableInelasticityString() { return "disableInelasticity"; }
  };


  void postInputInitialization() override;

  /// Name of the solid mechanics solver
  string m_solidSolverName;

  /// Pointer to the solid mechanics solver
  SolidMechanicsLagrangianFEM * m_solidSolver;

  /// Flag for resetting displacements (and velocities)
  int m_resetDisplacements;

  /// Flag for disabling / enabling inelasticity
  int m_disableInelasticity;
};


} /* namespace geos */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSTATERESET_HPP_ */
