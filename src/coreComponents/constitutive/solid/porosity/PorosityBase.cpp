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
 * @file PorosityBase.cpp
 */

#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{


PorosityBase::PorosityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newPorosity(),
  m_porosity_n(),
  m_dPorosity_dPressure(),
  m_dPorosity_dTemperature(),
  m_initialPorosity(),
  m_referencePorosity(),
  m_defaultReferencePorosity()
{
  registerWrapper( viewKeyStruct::defaultReferencePorosityString(), &m_defaultReferencePorosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the reference porosity" );

  registerField( fields::porosity::porosity{}, &m_newPorosity );

  registerField( fields::porosity::porosity_n{}, &m_porosity_n );

  registerField( fields::porosity::dPorosity_dPressure{}, &m_dPorosity_dPressure );

  registerField( fields::porosity::dPorosity_dTemperature{}, &m_dPorosity_dTemperature );

  registerField( fields::porosity::initialPorosity{}, &m_initialPorosity );

  registerField( fields::porosity::referencePorosity{}, &m_referencePorosity );
}

void PorosityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  m_newPorosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_porosity_n.resize( 0, numConstitutivePointsPerParentIndex );
  m_dPorosity_dPressure.resize( 0, numConstitutivePointsPerParentIndex );
  m_dPorosity_dTemperature.resize( 0, numConstitutivePointsPerParentIndex );
  m_initialPorosity.resize( 0, numConstitutivePointsPerParentIndex );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void PorosityBase::postInputInitialization()
{
  getField< fields::porosity::referencePorosity >().
    setApplyDefaultValue( m_defaultReferencePorosity );
}

void PorosityBase::scaleReferencePorosity( arrayView1d< real64 const > scalingFactors ) const
{
  localIndex const numE = numElem();

  arrayView1d< real64 > referencePorosity = m_referencePorosity;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    referencePorosity[k] *= scalingFactors[k];
  } );
}

void PorosityBase::saveConvergedState() const
{
  m_porosity_n.setValues< parallelDevicePolicy<> >( m_newPorosity.toViewConst() );
}

void PorosityBase::ignoreConvergedState() const
{
  m_newPorosity.setValues< parallelDevicePolicy<> >( m_porosity_n.toViewConst() );
}

void PorosityBase::initializeState() const
{
  m_porosity_n.setValues< parallelDevicePolicy<> >( m_newPorosity.toViewConst() );
  m_initialPorosity.setValues< parallelDevicePolicy<> >( m_newPorosity.toViewConst() );
}

}
} /* namespace geos */
