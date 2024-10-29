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
 * @file ExecutableGroup.hpp
 */

#ifndef GEOS_DATAREPOSITORY_EXECUTABLEGROUP_HPP_
#define GEOS_DATAREPOSITORY_EXECUTABLEGROUP_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "common/DataTypes.hpp"
#include "Group.hpp"
#include "mesh/DomainPartition.hpp"

namespace geos
{

// class DomainPartition;

/**
 * @class ExecutableGroup
 *
 * Objects that are executable and/or able to request changes dt
 * should be derived from this type.
 */
class ExecutableGroup : public dataRepository::Group
{
public:

  using dataRepository::Group::Group;

  /**
   * @brief Main extension point of executable targets.
   * @param[in] time_n        current time level
   * @param[in] dt            time step to be taken
   * @param[in] cycleNumber   global cycle number
   * @param[in] eventCounter  index of event that triggered execution
   * @param[in] eventProgress fractional progress in current cycle
   * @param[in,out] domain    the physical domain up-casted to a Group.
   * @return True iff this event requires exiting the event loop.
   *
   * @details If the start criteria are satisfied, then the event manager
   * will call this method.
   */
  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain ) = 0;

  /**
   * @brief Inform the object that it expects to execute during the next timestep.
   * @param[in] time_n        current time level
   * @param[in] dt            time step to be taken
   * @param[in] cycle         global cycle number
   * @param[in,out] domain    the physical domain
   */
  virtual void signalToPrepareForExecution( real64 const time_n,
                                            real64 const dt,
                                            integer const cycle,
                                            DomainPartition & domain );
  /**
   * @brief Called as the code exits the main run loop.
   * @param[in] time_n        current time level
   * @param[in] cycleNumber   global cycle number
   * @param[in] eventCounter  index of event that triggered execution
   * @param[in] eventProgress fractional progress in current cycle
   * @param[in,out] domain    the physical domain
   */
  virtual void cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain );

  /**
   * @brief Supplies the timestep request for this target to the event manager.
   * @param[in] time current time level
   * @return         desired time step size
   */
  virtual real64 getTimestepRequest( real64 const time )
  {
    GEOS_UNUSED_VAR( time );
    return 1e99;
  }

  /**
   * @brief Timestepping type.
   */
  enum class TimesteppingBehavior : integer
  {
    DeterminesTimeStepSize, ///< The group (say, the solver) does the timestepping
    DoesNotDetermineTimeStepSize ///< The event targetting this group does the timestepping
  };

  /**
   * @brief Set the timestep behavior for a target.
   * @param[in] timesteppingBehavior the timestepping behavior
   */
  void setTimesteppingBehavior( TimesteppingBehavior const timesteppingBehavior ) { m_timesteppingBehavior = timesteppingBehavior; }

  /**
   * @brief Get the target's time step behavior.
   * @return The time stepping type
   */
  TimesteppingBehavior getTimesteppingBehavior() const { return m_timesteppingBehavior; }

private:

  TimesteppingBehavior m_timesteppingBehavior = TimesteppingBehavior::DoesNotDetermineTimeStepSize;
};

/** @cond DO_NOT_DOCUMENT */
ENUM_STRINGS( ExecutableGroup::TimesteppingBehavior,
              "DeterminesTimeStepSize",
              "DoesNotDetermineTimeStepSize" );
/** @endcond */

}


#endif /* GEOS_DATAREPOSITORY_EXECUTABLEGROUP_HPP_ */
