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
 * @file HealMPMEvent.cpp
 */

#include "HealMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  HealMPMEvent::HealMPMEvent( const string & name,
                              Group * const parent ) :
                              MPMEventBase(  name, parent ),
                              m_targetRegion( "mat1" )
  {
    registerWrapper( viewKeyStruct::targetRegionString(), &m_targetRegion ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Particle region to perform heal on" );
  }

  HealMPMEvent::~HealMPMEvent() 
  {}

  void HealMPMEvent::postInputInitialization()
  {
    GEOS_LOG_RANK_0( "HealEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "targetRegion=" << m_targetRegion );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, HealMPMEvent, string const &, Group * const )

} /* namespace geos */
