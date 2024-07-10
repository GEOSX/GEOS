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
 *  @file CoulombContact.cpp
 */

#include "CoulombContact.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CoulombContact::CoulombContact( string const & name, Group * const parent ):
  ContactBase( name, parent ),
  m_cohesion(),
  m_frictionCoefficient(),
  m_elasticSlip()
{
  registerWrapper( viewKeyStruct::cohesionString(), &m_cohesion ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Cohesion" );

  registerWrapper( viewKeyStruct::frictionCoefficientString(), &m_frictionCoefficient ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Friction coefficient" );

  registerWrapper( viewKeyStruct::elasticSlipString(), &m_elasticSlip ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Elastic Slip" );

}

CoulombContact::~CoulombContact()
{}

void CoulombContact::postInputInitialization()
{
  GEOS_THROW_IF( m_frictionCoefficient < 0.0,
                 getFullName() << ": The provided friction coefficient is less than zero. Value: " << m_frictionCoefficient,
                 InputError );

}

void CoulombContact::allocateConstitutiveData( Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  m_elasticSlip.resize( 0, 2 );

  ContactBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


CoulombContactUpdates CoulombContact::createKernelWrapper() const
{
  return CoulombContactUpdates( m_penaltyStiffness,
                                m_shearStiffness,
                                m_displacementJumpThreshold,
                                *m_apertureTable,
                                m_cohesion,
                                m_frictionCoefficient,
                                m_elasticSlip );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CoulombContact, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
