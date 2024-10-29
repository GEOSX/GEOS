/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FrictionlessContact.cpp
 */

#include "FrictionlessContact.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

FrictionlessContact::FrictionlessContact( string const & name,
                                          Group * const parent ):
  ContactBase( name, parent )
{}

FrictionlessContact::~FrictionlessContact()
{}

void FrictionlessContact::allocateConstitutiveData( Group & parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  ContactBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

FrictionlessContactUpdates FrictionlessContact::createKernelWrapper() const
{
  return FrictionlessContactUpdates( m_penaltyStiffness,
                                     m_shearStiffness,
                                     m_displacementJumpThreshold,
                                     *m_apertureTable );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, FrictionlessContact, string const &, Group * const )

} /* end namespace constitutive */

} /* end namespace geos */
