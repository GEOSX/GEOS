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
 * @file DiffusionBase.cpp
 */

#include "constitutive/diffusion/DiffusionBase.hpp"
#include "constitutive/diffusion/DiffusionFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

DiffusionBase::DiffusionBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::defaultPhaseDiffusivityMultiplierString(), &m_defaultPhaseDiffusivityMultiplier ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "List of phase diffusivity multipliers" );

  registerField( fields::diffusion::diffusivity{}, &m_diffusivity );
  registerField( fields::diffusion::dDiffusivity_dTemperature{}, &m_dDiffusivity_dTemperature );
  registerField( fields::diffusion::phaseDiffusivityMultiplier{}, &m_phaseDiffusivityMultiplier );
}

void DiffusionBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  integer const numPhases = numFluidPhases();
  GEOS_THROW_IF_LT_MSG( numPhases, 2,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF_GT_MSG( numPhases, MAX_NUM_PHASES,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );

  GEOS_THROW_IF( numPhases != m_defaultPhaseDiffusivityMultiplier.size(),
                 GEOS_FMT( "{}: the arrays in `{}` and `{}` must have the same size",
                           getFullName(), viewKeyStruct::phaseNamesString(), viewKeyStruct::defaultPhaseDiffusivityMultiplierString() ),
                 InputError );

  m_diffusivity.resize( 0, 0, 3 );
  m_dDiffusivity_dTemperature.resize( 0, 0, 3 );
  m_phaseDiffusivityMultiplier.resize( 0, 0, 3 );
}

void DiffusionBase::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  // NOTE: enforcing 1 quadrature point
  m_diffusivity.resize( 0, 1, 3 );
  m_dDiffusivity_dTemperature.resize( 0, 1, 3 );
  m_phaseDiffusivityMultiplier.resize( 0, 1, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  for( localIndex ei = 0; ei < parent.size(); ++ei )
  {
    // NOTE: enforcing 1 quadrature point
    for( localIndex q = 0; q < 1; ++q )
    {
      for( integer ip = 0; ip < numFluidPhases(); ++ip )
      {
        m_phaseDiffusivityMultiplier[ei][q][ip] = m_defaultPhaseDiffusivityMultiplier[ip];
      }
    }
  }

}

} // namespace constitutive

} // namespace geos
