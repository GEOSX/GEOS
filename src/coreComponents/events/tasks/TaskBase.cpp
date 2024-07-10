/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TasksBase.hpp
 */

#include "TaskBase.hpp"

namespace geos
{

using namespace dataRepository;

TaskBase::TaskBase( string const & name,
                    Group * const parent ):
  ExecutableGroup( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

TaskBase::~TaskBase()
{ }

TaskBase::CatalogInterface::CatalogType & TaskBase::getCatalog()
{
  static TaskBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void TaskBase::postInputInitialization()
{ }

} /* namespace */
