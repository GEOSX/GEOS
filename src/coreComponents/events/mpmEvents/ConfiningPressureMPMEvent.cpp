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
 * @file ConfiningPressureMPMEvent.cpp
 */

#include "ConfiningPressureMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  ConfiningPressureMPMEvent::ConfiningPressureMPMEvent( const string & name,
                                                    Group * const parent ) :
                                                    MPMEventBase(  name, parent ),
                                                    m_confiningPressureBoxMin(  ),
                                                    m_confiningPressureBoxMax(  ),
                                                    m_startPressure( 0.0 ),
                                                    m_endPressure( 0.0 ),
                                                    m_interpType( 1 )
  {  
    registerWrapper( viewKeyStruct::confiningPressureBoxMinString(), &m_confiningPressureBoxMin ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Min corner of box defining confining pressure region" );

    registerWrapper( viewKeyStruct::confiningPressureBoxMaxString(), &m_confiningPressureBoxMax ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Max corner of box defining confining pressure region" );

    registerWrapper( viewKeyStruct::startPressureString(), &m_startPressure ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Starting pressure to ramp from" );

    registerWrapper( viewKeyStruct::endPressureString(), &m_endPressure ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "End pressure to ramp to" );

    registerWrapper( viewKeyStruct::interpTypeString(), &m_interpType ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue( m_interpType ).
        setDescription( "Interpolation scheme: 0 (Linear), 1 (Cosine), 2 (Smooth-step)" ); // CC: double check this!
  }

  ConfiningPressureMPMEvent::~ConfiningPressureMPMEvent() 
  {}

  void ConfiningPressureMPMEvent::postInputInitialization()
  {
    GEOS_ERROR_IF( m_interval < 0.0   , "Interval must be positive!");
    GEOS_ERROR_IF( m_confiningPressureBoxMin.size() != 3 , "confiningPressureBoxMin must be of length 3. ");
    GEOS_ERROR_IF( m_confiningPressureBoxMax.size() != 3 , "confiningPressureBoxMax must be of length 3. ");

    GEOS_LOG_RANK_0( "ConfiningPressureEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "confiningPressureBoxMin=" << m_confiningPressureBoxMin << ", " << 
                     "confiningPressureBoxMax=" << m_confiningPressureBoxMax << ", " << 
                     "startPressure=" << m_startPressure << ", " << 
                     "endPressure=" << m_endPressure << ", " << 
                     "interpType=" << m_interpType );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, ConfiningPressureMPMEvent, string const &, Group * const )

} /* namespace geos */
