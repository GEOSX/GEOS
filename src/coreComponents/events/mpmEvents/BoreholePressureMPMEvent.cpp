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
 * @file BoreholePressureMPMEvent.cpp
 */

#include "BoreholePressureMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  BoreholePressureMPMEvent::BoreholePressureMPMEvent( const string & name,
                                                    Group * const parent ) :
                                                    MPMEventBase(  name, parent ),
                                                    m_boreholeRadius( 0.0 ),
                                                    m_startPressure( 0.0 ),
                                                    m_endPressure( 0.0 ),
                                                    m_interpType( 1 )
  {  
    registerWrapper( viewKeyStruct::boreholeRadiusString(), &m_boreholeRadius ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Starting temperature to ramp from" );

    registerWrapper( viewKeyStruct::startPressureString(), &m_startPressure ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Starting temperature to ramp from" );

    registerWrapper( viewKeyStruct::endPressureString(), &m_endPressure ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "End temperature to ramp to" );

    registerWrapper( viewKeyStruct::interpTypeString(), &m_interpType ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue( m_interpType ).
        setDescription( "Interpolation scheme: 0 (Linear), 1 (Cosine), 2 (Smooth-step)" ); // CC: double check this!
  }

  BoreholePressureMPMEvent::~BoreholePressureMPMEvent() 
  {}

  void BoreholePressureMPMEvent::postInputInitialization()
  {
    GEOS_ERROR_IF( m_boreholeRadius < 0.0   , "Borehole radius must be positive!");

    GEOS_LOG_RANK_0( "BoreholePressureEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "boreholeRadius=" << m_boreholeRadius << ", " << 
                     "startPressure=" << m_startPressure << ", " << 
                     "endPressure=" << m_endPressure << ", " << 
                     "interpType=" << m_interpType );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, BoreholePressureMPMEvent, string const &, Group * const )

} /* namespace geos */
