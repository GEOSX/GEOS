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
 *  @file HyperelasticMMS.cpp
 */

#include "HyperelasticMMS.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

HyperelasticMMS::HyperelasticMMS( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultLambda(),
  m_defaultShearModulus(),
  m_lambda(),
  m_shearModulus()
{
  registerWrapper< real64 >( viewKeyStruct::defaultBulkModulusString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Bulk Modulus Parameter" );

  registerWrapper< real64 >( viewKeyStruct::defaultYoungModulusString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Young's Modulus" );

  registerWrapper< real64 >( viewKeyStruct::defaultPoissonRatioString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultLambdaString(), &m_defaultLambda ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Poisson's Ratio" );

  registerWrapper( viewKeyStruct::defaultShearModulusString(), &m_defaultShearModulus).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default shear modulus" );

  registerWrapper( viewKeyStruct::lambdaString(), &m_lambda ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::shearModulusString(), &m_shearModulus ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Shear Modulus Field" );
}

HyperelasticMMS::~HyperelasticMMS()
{}

void HyperelasticMMS::postInputInitialization()
{
  // check what constants the user actually input, and do conversions as needed

  SolidBase::postInputInitialization();

  real64 & nu = getReference< real64 >( viewKeyStruct::defaultPoissonRatioString() );
  real64 & E  = getReference< real64 >( viewKeyStruct::defaultYoungModulusString() );
  real64 & K  = getReference< real64 >( viewKeyStruct::defaultBulkModulusString() );
  real64 & lambda = m_defaultLambda;
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
  this->getWrapper< array1d< real64 > >( viewKeyStruct::lambdaString() ).
    setApplyDefaultValue( m_defaultLambda );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::shearModulusString() ).
    setApplyDefaultValue( m_defaultShearModulus );

  // This may not be needed since these are computed every call
  this->getWrapper< array2d< real64 > >( viewKeyStruct::wavespeedString() ).
    setApplyDefaultValue( sqrt( ( conversions::lameConstants::toBulkMod( m_defaultLambda, m_defaultShearModulus ) + (4.0/3.0) * m_defaultShearModulus ) / m_defaultDensity ) );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, HyperelasticMMS, string const &, Group * const )
}
} /* namespace geos */
