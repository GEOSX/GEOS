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
 * @file DispersionBase.cpp
 */

#include "constitutive/dispersion/DispersionBase.hpp"
#include "constitutive/dispersion/DispersionFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

DispersionBase::DispersionBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerField( fields::dispersion::dispersivity{}, &m_dispersivity );
}

void DispersionBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  m_dispersivity.resize( 0, 0, 3 );
}

void DispersionBase::allocateConstitutiveData( dataRepository::Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_dispersivity.resize( 0, 1, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

} // namespace constitutive

} // namespace geos
