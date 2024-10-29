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
 * @file FrictionCoefficientSwapMPMEvent.cpp
 */

#include "FrictionCoefficientSwapMPMEvent.hpp"

namespace geos
{

  using namespace dataRepository;
  
  FrictionCoefficientSwapMPMEvent::FrictionCoefficientSwapMPMEvent( const string & name,
                              Group * const parent ) :
                              MPMEventBase(  name, parent ),
                              m_frictionCoefficient( -1.0 ),
                              m_frictionCoefficientTable()
  {
    registerWrapper( "frictionCoefficient", &m_frictionCoefficient ).
        setApplyDefaultValue( -1 ).
        setInputFlag( InputFlags::OPTIONAL ).
        setDescription( "Coefficient of friction, currently assumed to be the same everywhere" );

    registerWrapper( "frictionCoefficientTable", &m_frictionCoefficientTable ).
        setInputFlag( InputFlags::OPTIONAL ).
        setDescription( "Friction coefficient table for different groups" );
  }

  FrictionCoefficientSwapMPMEvent::~FrictionCoefficientSwapMPMEvent() 
  {}

  void FrictionCoefficientSwapMPMEvent::postInputInitialization()
  {
    GEOS_LOG_RANK_0( "FrictionCoefficientSwapEvent: " << 
                     "Time=" << m_time << ", " << 
                     "Interval=" << m_interval );
  }

  REGISTER_CATALOG_ENTRY( MPMEventBase, FrictionCoefficientSwapMPMEvent, string const &, Group * const )

} /* namespace geos */
