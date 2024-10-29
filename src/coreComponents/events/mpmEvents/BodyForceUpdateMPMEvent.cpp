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
 * @file BodyForceUpdateMPMEvent.cpp
 */

#include "BodyForceUpdateMPMEvent.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

  using namespace dataRepository;
  
  BodyForceUpdateMPMEvent::BodyForceUpdateMPMEvent( const string & name,
                                  Group * const parent ) :
                                  MPMEventBase(  name, parent ),
                                  m_bodyForce()
  {  
    registerWrapper( viewKeyStruct::bodyForceString(), &m_bodyForce ).
        setInputFlag( InputFlags::OPTIONAL ).
        setDescription( "Body force vector" );
  }

  BodyForceUpdateMPMEvent::~BodyForceUpdateMPMEvent() 
  {}

  void BodyForceUpdateMPMEvent::postInputInitialization()
  {
    GEOS_ERROR_IF( m_bodyForce.size() != 3 && m_bodyForce.size() > 0,
                   "bodyForce must be of length 3. ");
    
    //Initialize body force if they're not specified by the user
    if( m_bodyForce.size() == 0)
    {
      m_bodyForce.resize(3);
      LvArray::tensorOps::fill< 3 >( m_bodyForce, 0.0 );
    }

    GEOS_LOG_RANK_0( "BodyForceUpdateEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "bodyForce=" << m_bodyForce );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, BodyForceUpdateMPMEvent, string const &, Group * const )

} /* namespace geos */
