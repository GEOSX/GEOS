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
 *  @file RateAndStateFriction.cpp
 */

#include "RateAndStateFriction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

RateAndStateFriction::RateAndStateFriction( string const & name, Group * const parent ):
  FrictionBase( name, parent )
{
  registerWrapper( viewKeyStruct::aCoefficientString(), &m_a ).
    setDescription( "Rate- and State-dependent friction coefficient a." );

  registerWrapper( viewKeyStruct::bCoefficientString(), &m_b ).
    setDescription( "Rate- and State-dependent friction coefficient b." );

  registerWrapper( viewKeyStruct::DcCoefficientString(), &m_Dc ).
    setDescription( "Rate- and State-dependent friction characteristic length." );

  registerWrapper( viewKeyStruct::referenceVelocityString(), &m_V0 ).
    setDescription( "Rate- and State-dependent friction reference slip rate." );

  registerWrapper( viewKeyStruct::referenceFrictionCoefficientString(), &m_mu0 ).
    setDescription( "Rate- and State-dependent friction reference friction coefficient." );

  registerWrapper( viewKeyStruct::frictionCoefficientString(), &m_frictionCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Friction coefficient" );

  /// Default values
  registerWrapper( viewKeyStruct::defaultACoefficientString(), &m_defaultA ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the Rate- and State-dependent friction coefficient a." );

  registerWrapper( viewKeyStruct::defaultBCoefficientString(), &m_defaultB ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the Rate- and State-dependent friction coefficient b." );

  registerWrapper( viewKeyStruct::defaultDcCoefficientString(), &m_defaultDc ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the Rate- and State-dependent friction characteristic length." );

  registerWrapper( viewKeyStruct::defaultReferenceVelocityString(), &m_defaultV0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the Rate- and State-dependent friction reference slip rate." );

  registerWrapper( viewKeyStruct::defaultReferenceFrictionCoefficientString(), &m_defaultMu0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the Rate- and State-dependent friction reference friction coefficient." );
}

RateAndStateFriction::~RateAndStateFriction()
{}

void RateAndStateFriction::postInputInitialization()
{
  this->getWrapper< array1d< real64 > >( viewKeyStruct::aCoefficientString() ).
    setApplyDefaultValue( m_defaultA );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::bCoefficientString() ).
    setApplyDefaultValue( m_defaultB );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::DcCoefficientString() ).
    setApplyDefaultValue( m_defaultDc );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::referenceVelocityString() ).
    setApplyDefaultValue( m_defaultV0 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::referenceFrictionCoefficientString() ).
    setApplyDefaultValue( m_defaultMu0 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::frictionCoefficientString() ).
    setApplyDefaultValue( m_defaultMu0 );
}

void RateAndStateFriction::allocateConstitutiveData( Group & parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  FrictionBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

using RateAndStateFrictionUpdates = RateAndStateFriction::KernelWrapper;
RateAndStateFrictionUpdates RateAndStateFriction::createKernelUpdates() const
{
  return RateAndStateFrictionUpdates( m_displacementJumpThreshold,
                                      m_frictionCoefficient, m_a, m_b,
                                      m_Dc, m_V0, m_mu0 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, RateAndStateFriction, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
