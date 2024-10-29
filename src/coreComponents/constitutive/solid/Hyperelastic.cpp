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
 *  @file Hyperelastic.cpp
 */

#include "Hyperelastic.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

Hyperelastic::Hyperelastic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultBulkModulus(),
  m_defaultShearModulus(),
  m_bulkModulus(),
  m_shearModulus()
{
  registerWrapper< real64 >( viewKeyStruct::defaultLambdaString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default First Lame Parameter" );

  registerWrapper< real64 >( viewKeyStruct::defaultYoungModulusString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus" );

  registerWrapper< real64 >( viewKeyStruct::defaultPoissonRatioString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultBulkModulusString(), &m_defaultBulkModulus).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default bulk modulus" );

  registerWrapper( viewKeyStruct::defaultShearModulusString(), &m_defaultShearModulus).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default shear modulus" );

  registerWrapper( viewKeyStruct::bulkModulusString(), &m_bulkModulus ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::shearModulusString(), &m_shearModulus ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Shear Modulus Field" );
}

Hyperelastic::~Hyperelastic()
{}

void Hyperelastic::postInputInitialization()
{
  // check what constants the user actually input, and do conversions as needed

  SolidBase::postInputInitialization();

  real64 & nu = getReference< real64 >( viewKeyStruct::defaultPoissonRatioString() );
  real64 & E  = getReference< real64 >( viewKeyStruct::defaultYoungModulusString() );
  real64 & lambda = getReference< real64 >( viewKeyStruct::defaultLambdaString() );
  real64 & K  = m_defaultBulkModulus;
  real64 & G  = m_defaultShearModulus;

  // Poisson ratio range is: -0.5 < nu < 0.5
  // Zero bulk, shear or Young modulus is not accepted to avoid devision to zero

  string errorCheck( "( " );
  int numConstantsSpecified = 0;
  if( nu > -0.5 && nu < 0.5 )
  {
    ++numConstantsSpecified;
    errorCheck += "nu, ";
  }
  if( E > 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "E, ";
  }
  if( K > 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "K, ";
  }
  if( G > 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "G, ";
  }
  if(lambda > 0.0)
  {
    ++numConstantsSpecified;
    errorCheck += "Lambda, ";
  }
  errorCheck += ")";

  GEOS_ERROR_IF( numConstantsSpecified != 2,
                 "A specific pair of elastic constants is required. Either (K,G), (K,E), (G,E), (K,nu), (G,nu) or (E,nu). "<<
                 "You have specified "<<errorCheck );

  if( nu > -0.5 && nu < 0.5 && E > 0.0 )
  {
    K = conversions::youngModAndPoissonRatio::toBulkMod( E, nu );
    G = conversions::youngModAndPoissonRatio::toShearMod( E, nu );
    lambda = conversions::youngModAndPoissonRatio::toFirstLame( E, nu); 
  }
  else if( nu > -0.5 && nu < 0.5 && G > 0.0 )
  {
    E = conversions::shearModAndPoissonRatio::toYoungMod( G, nu );
    K = conversions::shearModAndPoissonRatio::toBulkMod( G, nu );
    lambda = conversions::shearModAndPoissonRatio::toFirstLame( G, nu); 
  }
  else if( nu > -0.5 && nu < 0.5 && K > 0.0 )
  {
    E = conversions::bulkModAndPoissonRatio::toYoungMod( K, nu );
    G = conversions::bulkModAndPoissonRatio::toShearMod( K, nu );
    lambda = conversions::bulkModAndPoissonRatio::toFirstLame( K, nu);
  }
  else if( E > 0.0 && K > 0.0 )
  {
    nu = conversions::bulkModAndYoungMod::toPoissonRatio( K, E );
    G  = conversions::bulkModAndYoungMod::toShearMod( K, E );
    lambda = conversions::bulkModAndYoungMod::toFirstLame( K, E);
  }
  else if( E > 0.0 && G > 0.0 )
  {
    nu = conversions::shearModAndYoungMod::toPoissonRatio( G, E );
    K  = conversions::shearModAndYoungMod::toBulkMod( G, E );
    lambda = conversions::shearModAndYoungMod::toFirstLame( G, E);
  }
  else if( K > 0.0 && G > 0.0 )
  {
    E  = conversions::bulkModAndShearMod::toYoungMod( K, G );
    nu = conversions::bulkModAndShearMod::toPoissonRatio( K, G );
    lambda = conversions::bulkModAndShearMod::toFirstLame( K, G);
  }
  else if(lambda > 0.0 && E > 0.0){
    G = conversions::firstLameAndYoungMod::toShearMod(lambda, E);
    K = conversions::firstLameAndYoungMod::toBulkMod(lambda, E);
    nu = conversions::firstLameAndYoungMod::toPoissonRatio(lambda, E);
  } 
  else if(lambda > 0.0 && G > 0.0)
  {
    E = conversions::lameConstants::toYoungMod( lambda, G );
    K = conversions::lameConstants::toBulkMod( lambda, G );
    nu = conversions::lameConstants::toPoissonRatio( lambda, G );
  } 
  else if (lambda > 0.0 && K > 0.0) 
  {
    E = conversions::firstLameAndBulkMod::toYoungMod(lambda, K);
    G = conversions::firstLameAndBulkMod::toShearMod(lambda, K);
    nu = conversions::firstLameAndBulkMod::toPoissonRatio(lambda, K);
  }
  else if(lambda > 0.0 && nu > 0.-0.5 && nu < 0.5) 
  {
    E = conversions::firstLameAndPoissonRatio::toYoungMod(lambda, nu);
    K = conversions::firstLameAndPoissonRatio::toBulkMod(lambda, nu);
    G = conversions::firstLameAndPoissonRatio::toShearMod(lambda, nu);
  }
  else
  {
    GEOS_ERROR( "Invalid specification for default elastic constants. " << errorCheck << " has been specified." );
  }

  // set results as array default values
  this->getWrapper< array1d< real64 > >( viewKeyStruct::bulkModulusString() ).
    setApplyDefaultValue( m_defaultBulkModulus );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::shearModulusString() ).
    setApplyDefaultValue( m_defaultShearModulus );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, Hyperelastic, string const &, Group * const )
}
} /* namespace geos */
