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
 * @file SlurryFluidBase.cpp
 */

#include "SlurryFluidBase.hpp"

#include "SlurryFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SlurryFluidBase::SlurryFluidBase( string const & name, Group * const parent ):
  SingleFluidBase( name, parent ),
  m_isNewtonianFluid( true )
{

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "List of fluid component names" );

  registerWrapper( viewKeyStruct::defaultComponentDensityString(), &m_defaultComponentDensity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Default value for the component density." );

  registerWrapper( viewKeyStruct::defaultCompressibilityString(), &m_defaultComponentCompressibility ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Default value for the component compressibility." );

  registerWrapper( viewKeyStruct::defaultComponentViscosityString(), &m_defaultComponentViscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default value for the component viscosity." );

  registerWrapper( viewKeyStruct::flowBehaviorIndexString(), &m_nIndices ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flow behavior index" );

  registerWrapper( viewKeyStruct::flowConsistencyIndexString(), &m_Ks ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flow consistency index" );

  registerField( fields::slurryfluid::dDensity_dProppantConcentration{}, &m_dDensity_dProppantConc );
  registerField( fields::slurryfluid::dDensity_dComponentConcentration{}, &m_dDensity_dCompConc );

  registerField( fields::slurryfluid::fluidDensity{}, &m_fluidDensity );
  registerField( fields::slurryfluid::dFluidDensity_dPressure{}, &m_dFluidDens_dPres );
  registerField( fields::slurryfluid::dFluidDensity_dComponentConcentration{}, &m_dFluidDens_dCompConc );

  registerField( fields::slurryfluid::fluidViscosity{}, &m_fluidViscosity );
  registerField( fields::slurryfluid::dFluidViscosity_dPressure{}, &m_dFluidVisc_dPres );
  registerField( fields::slurryfluid::dFluidViscosity_dComponentConcentration{}, &m_dFluidVisc_dCompConc );

  registerField( fields::slurryfluid::componentDensity{}, &m_componentDensity );
  registerField( fields::slurryfluid::dComponentDensity_dPressure{}, &m_dCompDens_dPres );
  registerField( fields::slurryfluid::dComponentDensity_dComponentConcentration{}, &m_dCompDens_dCompConc );

  registerField( fields::slurryfluid::dViscosity_dProppantConcentration{}, &m_dViscosity_dProppantConc );
  registerField( fields::slurryfluid::dViscosity_dComponentConcentration{}, &m_dViscosity_dCompConc );
}

SlurryFluidBase::~SlurryFluidBase() = default;

void SlurryFluidBase::postInputInitialization()
{
  SingleFluidBase::postInputInitialization();

  localIndex const NC = numFluidComponents();

  GEOS_ERROR_IF( m_defaultComponentDensity.size() != NC,
                 getFullName() << ": The number of default density values is not the same as the component number" );

  GEOS_ERROR_IF( m_defaultComponentCompressibility.size() != NC,
                 getFullName() << ": The number of default compressibility values is not the same as the component number" );

  GEOS_ERROR_IF( m_defaultComponentViscosity.size() != NC,
                 getFullName() << ": The number of default viscosity values is not the same as the component number" );

}

localIndex SlurryFluidBase::numFluidComponents() const
{
  return LvArray::integerConversion< localIndex >( m_componentNames.size());
}

void SlurryFluidBase::allocateConstitutiveData( Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  SingleFluidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent.size() );

  localIndex const NC = numFluidComponents();

  m_dDensity_dProppantConc.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );

  m_componentDensity.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );
  m_dCompDens_dPres.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );
  m_dCompDens_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC, NC );

  m_fluidDensity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dFluidDens_dPres.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dFluidDens_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );

  m_fluidViscosity.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dFluidVisc_dPres.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dFluidVisc_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );

  m_dViscosity_dProppantConc.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dCompConc.resize( parent.size(), numConstitutivePointsPerParentIndex, NC );

}


} //namespace constitutive

} //namespace geos
