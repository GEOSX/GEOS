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
 * @file InitializeStressMPMEvent.cpp
 */

#include "InitializeStressMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  InitializeStressMPMEvent::InitializeStressMPMEvent( const string & name,
                                  Group * const parent ) :
                                  MPMEventBase(  name, parent ),
                                  m_pressure( 0.0 ),
                                  m_targetRegion( "mat1" )
  {  
    registerWrapper( viewKeyStruct::pressureString(), &m_pressure ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Starting temperature to ramp from" );
        
    registerWrapper( viewKeyStruct::targetRegionString(), &m_targetRegion ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Particle region to perform anneal on" );

  }

  InitializeStressMPMEvent::~InitializeStressMPMEvent() 
  {}

  void InitializeStressMPMEvent::postInputInitialization()
  {
    GEOS_LOG_RANK_0( "InitializeStressEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "pressure=" << m_pressure << ", " << 
                     "targetRegion=" << m_targetRegion );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, InitializeStressMPMEvent, string const &, Group * const )

} /* namespace geos */
