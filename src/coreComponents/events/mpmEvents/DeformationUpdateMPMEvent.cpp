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
 * @file DeformationUpdateMPMEvent.cpp
 */

#include "DeformationUpdateMPMEvent.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

  using namespace dataRepository;
  
  DeformationUpdateMPMEvent::DeformationUpdateMPMEvent( const string & name,
                                  Group * const parent ) :
                                  MPMEventBase(  name, parent ),
                                  m_prescribedFTable( 0 ),
                                  m_prescribedBoundaryFTable( 0 ),
                                  m_stressControl()
  {      
    registerWrapper( viewKeyStruct::prescribedFTableString(), &m_prescribedFTable ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue( 0 ).
        setDescription( "Flag to enable prescribed F table" );

    registerWrapper( viewKeyStruct::prescribedBoundaryFTableString(), &m_prescribedBoundaryFTable ).
        setInputFlag( InputFlags::OPTIONAL ).
        setApplyDefaultValue( 0 ).
        setDescription( "Flag to enable prescribed boundary F table" );

    registerWrapper( viewKeyStruct::stressControlString(), &m_stressControl ).
        setInputFlag( InputFlags::OPTIONAL ).
        setDescription( "Flags to enable stress control in x, y and z directions" );
  }

  DeformationUpdateMPMEvent::~DeformationUpdateMPMEvent() 
  {}

  void DeformationUpdateMPMEvent::postInputInitialization()
  {
    GEOS_ERROR_IF( m_stressControl.size() != 3 && m_stressControl.size() > 0,
                   "stressControl must be of length 3. ");

    //Initialize body force if they're not specified by the user
    if( m_stressControl.size() == 0)
    {
      m_stressControl.resize(3);
      LvArray::tensorOps::fill< 3 >( m_stressControl, 0 );
    }

    GEOS_LOG_RANK_0( "DeformationUpdateEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval << ", " << 
                     "prescribedBoundaryFTable=" << m_prescribedBoundaryFTable << ", "
                     "prescribedFTable=" << m_prescribedFTable << ", " <<
                     "stressControl=" << m_stressControl);
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, DeformationUpdateMPMEvent, string const &, Group * const )

} /* namespace geos */
