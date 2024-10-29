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
 * @file SolidBase.cpp
 */

#include "SolidBase.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

SolidBase::SolidBase( string const & name, Group * const parent ):
  ContinuumBase( name, parent ),
  m_thermalExpansionCoefficient()
{
  string const voightLabels[6] = { "XX", "YY", "ZZ", "YZ", "XZ", "XY" };

  registerWrapper( viewKeyStruct::defaultThermalExpansionCoefficientString(), &m_defaultThermalExpansionCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Linear Thermal Expansion Coefficient of the Solid Rock Frame" );

  registerWrapper( viewKeyStruct::thermalExpansionCoefficientString(), &m_thermalExpansionCoefficient ).
    setApplyDefaultValue( -1.0 ). // will be overwritten
    setDescription( "Linear Thermal Expansion Coefficient Field" );
}


SolidBase::~SolidBase()
{}


void SolidBase::postInputInitialization()
{
  ContinuumBase::postInputInitialization();

  this->getWrapper< array1d< real64 > >( viewKeyStruct::thermalExpansionCoefficientString() ).
    setApplyDefaultValue( m_defaultThermalExpansionCoefficient );
}


void SolidBase::allocateConstitutiveData( dataRepository::Group & parent,
                                          localIndex const numConstitutivePointsPerParentIndex )
{
  ContinuumBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_thermalExpansionCoefficient.resize( 0 );
}


void SolidBase::saveConvergedState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView3d< real64 const, solid::STRESS_USD > newStress = m_newStress;
  arrayView3d< real64, solid::STRESS_USD > oldStress = m_oldStress;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      LvArray::tensorOps::copy< 6 >( oldStress[k][q], newStress[k][q] );
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SolidBase, string const &, Group * const )
} /* namespace constitutive */
} /* namespace geos */
