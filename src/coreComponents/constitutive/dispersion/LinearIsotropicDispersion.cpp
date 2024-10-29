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
 * @file LinearIsotropicDispersion.cpp
 */

#include "LinearIsotropicDispersion.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

LinearIsotropicDispersion::LinearIsotropicDispersion( string const & name, Group * const parent ):
  DispersionBase( name, parent )
{
  registerWrapper( viewKeyStruct::longitudinalDispersivityString(), &m_longitudinalDispersivity ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Longitudinal dispersivity [m]" );
}

std::unique_ptr< ConstitutiveBase >
LinearIsotropicDispersion::deliverClone( string const & name,
                                         Group * const parent ) const
{
  return DispersionBase::deliverClone( name, parent );
}

void LinearIsotropicDispersion::postInputInitialization()
{
  GEOS_THROW_IF( m_longitudinalDispersivity < 0,
                 GEOS_FMT( "{}: longitudinal dispersivity must be positive",
                           getFullName() ),
                 InputError );
}

void LinearIsotropicDispersion::initializeVelocityState( arrayView2d< real64 const > const & initialVelocity ) const
{
  saveConvergedVelocityState( initialVelocity );
}

void LinearIsotropicDispersion::saveConvergedVelocityState( arrayView2d< real64 const > const & convergedVelocity ) const
{
  // note that the update function is called here, and not in the solver, because total velocity is treated explicitly
  KernelWrapper dispersionWrapper = createKernelWrapper();

  forAll< parallelDevicePolicy<> >( dispersionWrapper.numElems(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < dispersionWrapper.numGauss(); ++q )
    {
      dispersionWrapper.update( k, q, convergedVelocity[k] );
    }
  } );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearIsotropicDispersion, string const &, Group * const )

} // namespace constitutive

} // namespace geos
