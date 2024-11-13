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

/**
 *  @file CoulombFriction.cpp
 */

#include "CoulombFriction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CoulombFriction::CoulombFriction( string const & name, Group * const parent ):
  FrictionBase( name, parent ),
  m_cohesion(),
  m_frictionCoefficient(),
  m_elasticSlip()
{
  registerWrapper( viewKeyStruct::shearStiffnessString(), &m_shearStiffness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the shear elastic stiffness. Units of Pressure/length" );

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

CoulombFriction::~CoulombFriction()
{}

void CoulombFriction::postInputInitialization()
{
  GEOS_THROW_IF( m_frictionCoefficient < 0.0,
                 getFullName() << ": The provided friction coefficient is less than zero. Value: " << m_frictionCoefficient,
                 InputError );

}

void CoulombFriction::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  m_elasticSlip.resize( 0, 2 );

  FrictionBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


CoulombFrictionUpdates CoulombFriction::createKernelUpdates() const
{
  return CoulombFrictionUpdates( m_displacementJumpThreshold,
                                 m_shearStiffness,
                                 m_cohesion,
                                 m_frictionCoefficient,
                                 m_elasticSlip );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CoulombFriction, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
