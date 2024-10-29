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
 *  @file ElasticTransverseIsotropic.cpp
 */

#include "ElasticTransverseIsotropic.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

ElasticTransverseIsotropic::ElasticTransverseIsotropic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultYoungModulusTransverse(),
  m_defaultYoungModulusAxial(),
  m_defaultPoissonRatioTransverse(),
  m_defaultPoissonRatioAxialTransverse(),
  m_defaultShearModulusAxialTransverse(),
  m_c11(),
  m_c13(),
  m_c33(),
  m_c44(),
  m_c66(),
  m_effectiveBulkModulus(),
  m_effectiveShearModulus(),
  m_materialDirection()
{
  registerWrapper( viewKeyStruct::defaultYoungModulusTransverseString(), &m_defaultYoungModulusTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Transverse Young's Modulus" );

  registerWrapper( viewKeyStruct::defaultYoungModulusAxialString(), &m_defaultYoungModulusAxial ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial Young's Modulus" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioTransverseString(), &m_defaultPoissonRatioTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Transverse Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultPoissonRatioAxialTransverseString(), &m_defaultPoissonRatioAxialTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial-Transverse Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultShearModulusAxialTransverseString(), &m_defaultShearModulusAxialTransverse ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Axial-Transverse Shear Modulus" );

  registerWrapper< real64 >( viewKeyStruct::defaultC11String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C11" );

  registerWrapper< real64 >( viewKeyStruct::defaultC13String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C13" );

  registerWrapper< real64 >( viewKeyStruct::defaultC33String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C33" );

  registerWrapper< real64 >( viewKeyStruct::defaultC44String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C44" );

  registerWrapper< real64 >( viewKeyStruct::defaultC66String() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Stiffness Parameter C66" );

  registerWrapper( viewKeyStruct::c11String(), &m_c11 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C11" );

  registerWrapper( viewKeyStruct::c13String(), &m_c13 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C13" );

  registerWrapper( viewKeyStruct::c33String(), &m_c33 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C33" );

  registerWrapper( viewKeyStruct::c44String(), &m_c44 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C44" );

  registerWrapper( viewKeyStruct::c66String(), &m_c66 ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Stiffness Field C66" );

  registerWrapper( viewKeyStruct::effectiveBulkModulusString(), &m_effectiveBulkModulus ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Effective bulk modulus for stress control and wavespeed calculations" );
  
  registerWrapper( viewKeyStruct::effectiveShearModulusString(), &m_effectiveShearModulus ).
    setInputFlag( InputFlags::FALSE).
    setDescription( "Effective shear modulus for stress control and wavespeed calculations");

  registerWrapper( viewKeyStruct::materialDirectionString(), &m_materialDirection ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Material direction" );
}

ElasticTransverseIsotropic::~ElasticTransverseIsotropic()
{}

void ElasticTransverseIsotropic::allocateConstitutiveData( dataRepository::Group & parent, 
                                                           localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_materialDirection.resize( 0, 3 );
}

void ElasticTransverseIsotropic::postInputInitialization()
{
  SolidBase::postInputInitialization();

  real64 & c11  = getReference< real64 >( viewKeyStruct::defaultC11String() );
  real64 & c13  = getReference< real64 >( viewKeyStruct::defaultC13String() );
  real64 & c33  = getReference< real64 >( viewKeyStruct::defaultC33String() );
  real64 & c44  = getReference< real64 >( viewKeyStruct::defaultC44String() );
  real64 & c66  = getReference< real64 >( viewKeyStruct::defaultC66String() );

  real64 & Et = m_defaultYoungModulusTransverse;
  real64 & Ea = m_defaultYoungModulusAxial;
  real64 & Nut = m_defaultPoissonRatioTransverse;
  real64 & Nuat = m_defaultPoissonRatioAxialTransverse;
  real64 & Gat = m_defaultShearModulusAxialTransverse;

  // CC: TODO make sure that if stiffness constants are set or issues with other variable error is through or other logic to ensure effective bulk and shear moduli can be computed
  if( Et > 0.0 && Ea > 0.0 && Gat > 0.0 && Nut > -0.5 && Nut < 0.5 )
  {
    real64 const Nuta = Nuat * ( Et / Ea );
    real64 const delta = ( 1.0 + Nut ) * ( 1.0 - Nut - 2.0 * Nuta * Nuat );

    if( delta > 0.0 && Nuta * Nuat < 1.0 )
    {
      c11 = ( 1.0 - Nuta * Nuat ) * Et / delta;
      c13 = Nuat * ( 1.0 + Nut ) * Et / delta;
      c33 = ( 1.0 - Nut * Nut ) * Ea / delta;
      c44 = 2.0 * Gat;
      c66 = Et / ( 1.0 + Nut );
    }
  } 
  else 
  {
    Et = 4 * c66 * (c11 * c33 - c66 * c33 - c13 * c13) / ( c11 * c33 - c13 * c13 );
    Ea = c33 - c13 * c13 / ( c11 - c66 );
    Gat = c44 / 2.0;
    Nut = 4 * (c11 * c33 - c66 * c33 - c13 * c13 ) / ( c11 * c33 - c13 * c13 ) - 1;
    Nuat = c13 / ( 2 * ( c11 - c66 ) );
  }

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c11String() ).
    setApplyDefaultValue( c11 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c13String() ).
    setApplyDefaultValue( c13 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c33String() ).
    setApplyDefaultValue( c33 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c44String() ).
    setApplyDefaultValue( c44 );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::c66String() ).
    setApplyDefaultValue( c66 );

  real64 Keff = -Et*Ea/(2*Ea*(Nut+Nuat-1) + Et*(2*Nuat-1));
  real64 Geff=  0.6*Keff;

  this->getWrapper< array1d< real64 > >( viewKeyStruct::effectiveBulkModulusString() ).
    setApplyDefaultValue( Keff );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::effectiveShearModulusString() ).
    setApplyDefaultValue( Geff );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticTransverseIsotropic, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
