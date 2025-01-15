/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
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
  : ConstitutiveBase( name, parent ),
  m_numDOF( 1 )
{
  registerField( fields::singlefluid::density{}, &m_density.value );
  registerField( fields::singlefluid::dDensity{}, &m_density.derivs );
  registerField( fields::singlefluid::density_n{}, &m_density_n );

  registerField( fields::singlefluid::viscosity{}, &m_viscosity.value );
  registerField( fields::singlefluid::dViscosity{}, &m_viscosity.derivs );

  registerField( fields::singlefluid::internalEnergy{}, &m_internalEnergy.value );
  registerField( fields::singlefluid::dInternalEnergy{}, &m_internalEnergy.derivs );
  registerField( fields::singlefluid::internalEnergy_n{}, &m_internalEnergy_n );

  registerField( fields::singlefluid::enthalpy{}, &m_enthalpy.value );
  registerField( fields::singlefluid::dEnthalpy{}, &m_enthalpy.derivs );
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
  localIndex const numElem = m_density.value.size( 0 );
  localIndex const numGauss = m_density.value.size( 1 );

  SingleFluidProp::ViewTypeConst const density = m_density.toViewConst();
  SingleFluidProp::ViewTypeConst const internalEnergy = m_internalEnergy.toViewConst();

  arrayView2d< real64, singlefluid::USD_FLUID > const density_n = m_density_n.toView();
  arrayView2d< real64, singlefluid::USD_FLUID > const internalEnergy_n = m_internalEnergy_n.toView();

  forAll< parallelDevicePolicy<> >( numElem, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numGauss; ++q )
    {
      density_n[k][q] = density.value[k][q];
      internalEnergy_n[k][q] = internalEnergy.value[k][q];
    }
  } );
}

//START_SPHINX_INCLUDE_00
void SingleFluidBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  resize( parent.size() );

  // density
  m_density.value.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_density.derivs.resize( parent.size(), numConstitutivePointsPerParentIndex, m_numDOF );
  m_density_n.resize( parent.size(), numConstitutivePointsPerParentIndex );

  // viscosity
  m_viscosity.value.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_viscosity.derivs.resize( parent.size(), numConstitutivePointsPerParentIndex, m_numDOF );

  // internal energy
  m_internalEnergy.value.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_internalEnergy.derivs.resize( parent.size(), numConstitutivePointsPerParentIndex, m_numDOF );
  m_internalEnergy_n.resize( parent.size(), numConstitutivePointsPerParentIndex );

  // enthalpy
  m_enthalpy.value.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_enthalpy.derivs.resize( parent.size(), numConstitutivePointsPerParentIndex, m_numDOF );
}
//END_SPHINX_INCLUDE_00

} //namespace constitutive

} //namespace geos
