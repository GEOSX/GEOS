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
  FrictionBase( name, parent )
{}

FrictionlessContact::~FrictionlessContact()
{}

FrictionlessContactUpdates FrictionlessContact::createKernelWrapper() const
{
  return FrictionlessContactUpdates( m_displacementJumpThreshold );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, FrictionlessContact, string const &, Group * const )

} /* end namespace constitutive */

} /* end namespace geos */
