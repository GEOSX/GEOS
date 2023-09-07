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
 * @file AnnealMPMEvent.cpp
 */

#include "AnnealMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  AnnealMPMEvent::AnnealMPMEvent( const string & name,
                                  Group * const parent ) :
                                  MPMEventBase(  name, parent ),
                                  m_source( "mat1" )
  {  
    registerWrapper( viewKeyStruct::sourceString(), &m_source ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Particle region to perform anneal on" );
  }

  AnnealMPMEvent::~AnnealMPMEvent() 
  {}

  void AnnealMPMEvent::postProcessInput()
  {
    GEOS_LOG_RANK_0( "AnnealEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "Source=" << m_source );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, AnnealMPMEvent, string const &, Group * const )

} /* namespace geos */
