/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ExecutableGroup.cpp
 */

#include "ExecutableGroup.hpp"

namespace geos
{

void ExecutableGroup::signalToPrepareForExecution( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                   real64 const GEOS_UNUSED_PARAM( dt ),
                                                   integer const GEOS_UNUSED_PARAM( cycle ),
                                                   DomainPartition &
                                                   GEOS_UNUSED_PARAM( domain ) )
{}

void ExecutableGroup::cleanup( real64 const GEOS_UNUSED_PARAM( time_n ),
                               integer const GEOS_UNUSED_PARAM( cycleNumber ),
                               integer const GEOS_UNUSED_PARAM( eventCounter ),
                               real64 const GEOS_UNUSED_PARAM( eventProgress ),
                               DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{}

}
