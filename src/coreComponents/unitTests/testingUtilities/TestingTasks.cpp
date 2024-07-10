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
 * @file TestingTasks.cpp
 */

#include "TestingTasks.hpp"

namespace geos
{
namespace testing
{


TimeStepChecker::TimeStepChecker( string const & name, Group * const parent ):
  TaskBase( name, parent )
{}

bool TimeStepChecker::execute( real64 const time_n,
                               real64 const GEOS_UNUSED_PARAM( dt ),
                               integer const GEOS_UNUSED_PARAM( cycleNumber ),
                               integer const GEOS_UNUSED_PARAM( eventCounter ),
                               real64 const GEOS_UNUSED_PARAM( eventProgress ),
                               DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  EXPECT_TRUE( m_checkTimeStepFunction );
  m_checkTimeStepFunction( time_n );

  ++m_timestepId;
  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase, TimeStepChecker, string const &, geos::dataRepository::Group * const )


} // namespace testing

} // namespace geos
