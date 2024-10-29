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
 * @file InsertPeriodicContactSurfacesMPMEvent.cpp
 */

#include "InsertPeriodicContactSurfacesMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  InsertPeriodicContactSurfacesMPMEvent::InsertPeriodicContactSurfacesMPMEvent( const string & name,
                                                                                Group * const parent ) :
                                                                                MPMEventBase(  name, parent )
  {}

  InsertPeriodicContactSurfacesMPMEvent::~InsertPeriodicContactSurfacesMPMEvent() 
  {}

  REGISTER_CATALOG_ENTRY( MPMEventBase, InsertPeriodicContactSurfacesMPMEvent, string const &, Group * const )

} /* namespace geos */
