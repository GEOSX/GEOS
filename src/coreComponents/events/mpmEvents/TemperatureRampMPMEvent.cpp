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
 * @file TemperatureRampMPMEvent.cpp
 */

#include "TemperatureRampMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  TemperatureRampMPMEvent::TemperatureRampMPMEvent( const string & name,
                                                    Group * const parent ) :
                                                    MPMEventBase(  name, parent ),
                                                    m_startTemperature( 0.0 ),
                                                    m_endTemperature( 0.0 ),
                                                    m_interpType( 1 )
  {  
    registerWrapper( viewKeyStruct::startTemperatureString(), &m_startTemperature ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Starting temperature to ramp from" );

    registerWrapper( viewKeyStruct::endTemperatureString(), &m_endTemperature ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "End temperature to ramp to" );

    registerWrapper( viewKeyStruct::interpTypeString(), &m_interpType ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue( m_interpType ).
        setDescription( "Interpolation scheme: 0 (Linear), 1 (Cosine), 2 (Sigmoid)" ); // CC: double check this!
  }

  TemperatureRampMPMEvent::~TemperatureRampMPMEvent() 
  {}

  void TemperatureRampMPMEvent::postInputInitialization()
  {
    GEOS_ERROR_IF( m_startTemperature < 0.0 || m_endTemperature < 0.0  , "Temperatures must be positive!");

    GEOS_LOG_RANK_0( "TemperatureRampEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "startTemperature=" << m_startTemperature << ", " << 
                     "endTemperature=" << m_endTemperature << ", " << 
                     "interpType=" << m_interpType );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, TemperatureRampMPMEvent, string const &, Group * const )

} /* namespace geos */
