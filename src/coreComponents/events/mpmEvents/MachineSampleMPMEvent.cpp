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
 * @file MachineSampleMPMEvent.cpp
 */

#include "MachineSampleMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  MachineSampleMPMEvent::MachineSampleMPMEvent( const string & name,
                                              Group * const parent ) :
                                              MPMEventBase(  name, parent ),
                                              m_sampleType( "dogbone" ),
                                              m_filletRadius( -1 ),
                                              m_gaugeLength( -1 ),
                                              m_gaugeRadius( -1 ),
                                              m_diskRadius( -1 )

  {
    registerWrapper( viewKeyStruct::sampleTypeString(), &m_sampleType ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Type of sample to machine" );

    registerWrapper( viewKeyStruct::filletRadiusString(), &m_filletRadius ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue(-1).
        setDescription( "Fillet radius for machining dogbone sample" );

    registerWrapper( viewKeyStruct::gaugeLengthString(), &m_gaugeLength ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue(-1).
        setDescription( "Gauge length for machining dogbone sample" );

    registerWrapper( viewKeyStruct::gaugeRadiusString(), &m_gaugeRadius ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue(-1).
        setDescription( "Gauge radius for machining dogbone sample" );

    registerWrapper( viewKeyStruct::diskRadiusString(), &m_diskRadius ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue(-1).
        setDescription( "Disk radius for machining brazil disk sample" );
  }

  MachineSampleMPMEvent::~MachineSampleMPMEvent() 
  {}

  void MachineSampleMPMEvent::postInputInitialization()
  {
    if( m_sampleType == "dogbone" )
    {
        GEOS_ERROR_IF( m_filletRadius < 0.0, "Fillet radius must be specified for dogbone sample type and positive" );
        GEOS_ERROR_IF( m_gaugeLength < 0.0, "Gauge length must be specified for dogbone sample type and positive" );
        GEOS_ERROR_IF( m_gaugeRadius < 0.0, "Gauge radius must be specified for dogbone sample type and positive" );
    }

    if( m_sampleType == "brazilDisk")
    {
        GEOS_ERROR_IF( m_diskRadius < 0.0, "Disk radius must be specified for brazil disk sample type and positive" ); 
    }

    if( m_sampleType == "cylinder" )
    {
        GEOS_ERROR_IF( m_gaugeRadius < 0.0, "Gauge radius must be specified for cylinder sample type and positive" );
    }
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, MachineSampleMPMEvent, string const &, Group * const )

} /* namespace geos */
