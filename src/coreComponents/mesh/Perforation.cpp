/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file Perforation.cpp
 */

#include "Perforation.hpp"
#include "dataRepository/InputFlags.hpp"

namespace geos
{

using namespace dataRepository;

Perforation::Perforation( string const & name, Group * const parent )
  : Group( name, parent ),
  m_distanceFromHead( 0 ),
  m_wellTransmissibility( 0 ),
  m_wellSkinFactor( 0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::distanceFromHeadString(), &m_distanceFromHead ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Linear distance from well head to the perforation" );

  registerWrapper( viewKeyStruct::wellTransmissibilityString(), &m_wellTransmissibility ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Perforation transmissibility" );

  registerWrapper( viewKeyStruct::wellSkinFactorString(), &m_wellSkinFactor ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Perforation skin factor" );
}


void Perforation::postInputInitialization()
{
  GEOS_ERROR_IF( m_distanceFromHead <= 0,
                 getWrapperDataContext( viewKeyStruct::distanceFromHeadString() ) <<
                 ": distance from well head to perforation cannot be negative." );
}


Perforation::~Perforation()
{}

} //namespace geos
