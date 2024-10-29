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
 * @file ConstantDiffusion.cpp
 */

#include "ConstantDiffusion.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ConstantDiffusion::ConstantDiffusion( string const & name, Group * const parent ):
  DiffusionBase( name, parent )
{
  registerWrapper( viewKeyStruct::diffusivityComponentsString(), &m_diffusivityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy, and zz components of a diffusivity tensor [m^2/s]" );
}

std::unique_ptr< ConstitutiveBase >
ConstantDiffusion::deliverClone( string const & name,
                                 Group * const parent ) const
{
  return DiffusionBase::deliverClone( name, parent );
}

void ConstantDiffusion::allocateConstitutiveData( dataRepository::Group & parent,
                                                  localIndex const numConstitutivePointsPerParentIndex )
{
  DiffusionBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    // NOTE: enforcing 1 quadrature point
    for( localIndex q = 0; q < 1; ++q )
    {
      m_diffusivity[ei][q][0] = m_diffusivityComponents[0];
      m_diffusivity[ei][q][1] = m_diffusivityComponents[1];
      m_diffusivity[ei][q][2] = m_diffusivityComponents[2];
    }
  }
}

void ConstantDiffusion::postInputInitialization()
{
  GEOS_THROW_IF( m_diffusivityComponents.size() != 3,
                 GEOS_FMT( "{}: the size of the diffusivity must be equal to 3",
                           getFullName() ),
                 InputError );

  GEOS_THROW_IF( m_diffusivityComponents[0] < 0 ||
                 m_diffusivityComponents[1] < 0 ||
                 m_diffusivityComponents[2] < 0,
                 GEOS_FMT( "{}: the components of the diffusivity tensor must be non-negative",
                           getFullName() ),
                 InputError );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ConstantDiffusion, string const &, Group * const )

} // namespace constitutive

} // namespace geos
