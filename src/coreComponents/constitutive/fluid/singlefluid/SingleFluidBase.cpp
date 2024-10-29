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
 * @file SingleFluidBase.cpp
 */

#include "SingleFluidBase.hpp"

#include "SingleFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SingleFluidBase::SingleFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerField( fields::singlefluid::density{}, &m_density );
  registerField( fields::singlefluid::dDensity_dPressure{}, &m_dDensity_dPressure );
  registerField( fields::singlefluid::dDensity_dTemperature{}, &m_dDensity_dTemperature );
  registerField( fields::singlefluid::density_n{}, &m_density_n );

  registerField( fields::singlefluid::viscosity{}, &m_viscosity );
  registerField( fields::singlefluid::dViscosity_dPressure{}, &m_dViscosity_dPressure );
  registerField( fields::singlefluid::dViscosity_dTemperature{}, &m_dViscosity_dTemperature );

  registerField( fields::singlefluid::internalEnergy{}, &m_internalEnergy );
  registerField( fields::singlefluid::internalEnergy_n{}, &m_internalEnergy_n );
  registerField( fields::singlefluid::dInternalEnergy_dPressure{}, &m_dInternalEnergy_dPressure );
  registerField( fields::singlefluid::dInternalEnergy_dTemperature{}, &m_dInternalEnergy_dTemperature );

  registerField( fields::singlefluid::enthalpy{}, &m_enthalpy );
  registerField( fields::singlefluid::dEnthalpy_dPressure{}, &m_dEnthalpy_dPressure );
  registerField( fields::singlefluid::dEnthalpy_dTemperature{}, &m_dEnthalpy_dTemperature );

}

void SingleFluidBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  // for fracture elements, set the default value
  getField< fields::singlefluid::density_n >().
    setDefaultValue( defaultDensity() );
}

void SingleFluidBase::initializeState() const
{
  saveConvergedState();
}

void SingleFluidBase::saveConvergedState() const
{
  m_density_n.setValues< parallelDevicePolicy<> >( m_density.toViewConst() );
  m_internalEnergy_n.setValues< parallelDevicePolicy<> >( m_internalEnergy.toViewConst() );
}

//START_SPHINX_INCLUDE_00
void SingleFluidBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  resize( parent.size() );

  m_density.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dTemperature.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_density_n.resize( parent.size(), numConstitutivePointsPerParentIndex );

  m_viscosity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dTemperature.resize( parent.size(), numConstitutivePointsPerParentIndex );

  m_internalEnergy.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_internalEnergy_n.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dInternalEnergy_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dInternalEnergy_dTemperature.resize( parent.size(), numConstitutivePointsPerParentIndex );

  m_enthalpy.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dEnthalpy_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dEnthalpy_dTemperature.resize( parent.size(), numConstitutivePointsPerParentIndex );
}
//END_SPHINX_INCLUDE_00

} //namespace constitutive

} //namespace geos
