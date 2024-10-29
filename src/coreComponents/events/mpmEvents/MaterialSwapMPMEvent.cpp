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
 * @file MaterialSwapMPMEvent.cpp
 */

#include "MaterialSwapMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  MaterialSwapMPMEvent::MaterialSwapMPMEvent( const string & name,
                                              Group * const parent ) :
                                              MPMEventBase(  name, parent ),
                                              m_sourceRegion( "mat1" ),
                                              m_destinationRegion( "mat2" )
  {
    registerWrapper( viewKeyStruct::sourceRegionString(), &m_sourceRegion ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Particle region to transfer particles from" );

    registerWrapper( viewKeyStruct::destinationRegionString(), &m_destinationRegion ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Particle region to transfer particles to" );
  }

  MaterialSwapMPMEvent::~MaterialSwapMPMEvent() 
  {}

  void MaterialSwapMPMEvent::postInputInitialization()
  {
    GEOS_LOG_RANK_0( "MaterialSwapEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "Source region=" << m_sourceRegion << ", " << 
                     "Destination region=" << m_destinationRegion );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, MaterialSwapMPMEvent, string const &, Group * const )

} /* namespace geos */
