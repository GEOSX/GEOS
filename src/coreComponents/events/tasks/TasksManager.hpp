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
 * @file TasksManager.hpp
 */

#ifndef GEOS_EVENTS_TASKS_TASKSMANAGER_HPP_
#define GEOS_EVENTS_TASKS_TASKSMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "TaskBase.hpp"

namespace geos
{

/**
 * @class TasksManager
 * @brief A class to manage and execute tasks.
 */
class TasksManager : public dataRepository::Group
{
public:
  /// @copydoc geos::dataRepository::Group::Group
  TasksManager( string const & name, Group * const parent );
  /// Destructor
  virtual ~TasksManager() override;

  /// @copydoc geos::dataRepository::Group::createChild
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void expandObjectCatalogs() override;

private:
  TasksManager() = delete;

};

} /* namespace geos */

#endif
