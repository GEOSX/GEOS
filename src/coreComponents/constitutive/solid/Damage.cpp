
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
 * @file Damage.cpp
 */

#include "Damage.hpp"

namespace geos
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
Damage< BASE >::Damage( string const & name, Group * const parent ):
  BASE( name, parent ),
  m_newDamage(),
  m_oldDamage(),
  m_damageGrad(),
  m_strainEnergyDensity(),
  m_volStrain(),
  m_extDrivingForce(),
  m_lengthScale(),
  m_defaultCriticalFractureEnergy(),
  m_criticalStrainEnergy(),
  m_degradationLowerLimit( 0.0 ),
  m_extDrivingForceFlag( 0 ),
  m_defaultTensileStrength(),
  m_defaultCompressStrength(),
  m_defaultDeltaCoefficient(),
  m_damagePressure(),
  m_biotCoefficient(),
  m_criticalFractureEnergy(),
  m_tensileStrength(),
  m_compressStrength(),
  m_deltaCoefficient()
{
  this->registerWrapper( viewKeyStruct::newDamageString(), &m_newDamage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Material New Damage Variable" );

  this->registerWrapper( viewKeyStruct::oldDamageString(), &m_oldDamage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Material Old Damage Variable" );

  this->registerWrapper( viewKeyStruct::damageGradString(), &m_damageGrad ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Material Damage Gradient" );

  this->registerWrapper( viewKeyStruct::strainEnergyDensityString(), &m_strainEnergyDensity ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Strain Energy Density" );

  this->registerWrapper( viewKeyStruct::volumetricStrainString(), &m_volStrain ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Volumetric strain" );

  this->registerWrapper( viewKeyStruct::extDrivingForceString(), &m_extDrivingForce ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "External Driving Force" );

  this->registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Length scale l in the phase-field equation" );

  this->registerWrapper( viewKeyStruct::defaultCriticalFractureEnergyString(), &m_defaultCriticalFractureEnergy ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default critical fracture energy" );

  this->registerWrapper( viewKeyStruct::criticalFractureEnergyString(), &m_criticalFractureEnergy ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Critical fracture energy" );

  this->registerWrapper( viewKeyStruct::criticalStrainEnergyString(), &m_criticalStrainEnergy ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Critical stress in a 1d tension test" );

  this->registerWrapper( viewKeyStruct::degradationLowerLimitString(), &m_degradationLowerLimit ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The lower limit of the degradation function" );

  this->registerWrapper( viewKeyStruct::extDrivingForceFlagString(), &m_extDrivingForceFlag ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Whether to have external driving force. Can be 0 or 1" );

  this->registerWrapper( viewKeyStruct::defaultTensileStrengthString(), &m_defaultTensileStrength ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default tensile strength from the uniaxial tension test" );

  this->registerWrapper( viewKeyStruct::tensileStrengthString(), &m_tensileStrength ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Tensile strength from the uniaxial tension test" );

  this->registerWrapper( viewKeyStruct::defaultCompressStrengthString(), &m_defaultCompressStrength ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default compressive strength from the uniaxial compression test" );
  
  this->registerWrapper( viewKeyStruct::compressStrengthString(), &m_compressStrength ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Compressive strength from the uniaxial compression test" );

  this->registerWrapper( viewKeyStruct::defaultDeltaCoefficientString(), &m_defaultDeltaCoefficient ).
    setApplyDefaultValue( -1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default coefficient in the calculation of the external driving force" );

  this->registerWrapper( viewKeyStruct::deltaCoefficientString(), &m_deltaCoefficient ).
    setApplyDefaultValue( -1.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Coefficient in the calculation of the external driving force" );

  this->registerWrapper( viewKeyStruct::damagePressureString(), &m_damagePressure ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "A prescribed uniform pressure in the phase-field crack" );

  this->registerWrapper( viewKeyStruct::biotCoefficientString(), &m_biotCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Biot coefficient" );
}


template< typename BASE >
void Damage< BASE >::postProcessInput()
{
  BASE::postProcessInput();

  GEOS_ERROR_IF( m_extDrivingForceFlag != 0 && m_extDrivingForceFlag!= 1,
                 BASE::getDataContext() << ": invalid external driving force flag option - must"
                                           " be 0 or 1" );
  GEOS_ERROR_IF( m_extDrivingForceFlag == 1 && m_defaultTensileStrength <= 0.0,
                 BASE::getDataContext() << ": tensile strength must be input and positive when the"
                                           " external driving force flag is turned on" );
  GEOS_ERROR_IF( m_extDrivingForceFlag == 1 && m_defaultCompressStrength  <= 0.0,
                 BASE::getDataContext() << ": compressive strength must be input and positive when the"
                                           " external driving force flag is turned on" );
  GEOS_ERROR_IF( m_extDrivingForceFlag == 1 && m_defaultDeltaCoefficient < 0.0,
                 BASE::getDataContext() << ": delta coefficient must be input and non-negative when the"
                                           " external driving force flag is turned on" );

  // set results as array default values
  this->template getWrapper< array1d< real64 > >( viewKeyStruct::criticalFractureEnergyString() ).
    setApplyDefaultValue( m_defaultCriticalFractureEnergy );

  this->template getWrapper< array1d< real64 > >( viewKeyStruct::tensileStrengthString() ).
    setApplyDefaultValue( m_defaultTensileStrength );

  this->template getWrapper< array1d< real64 > >( viewKeyStruct::compressStrengthString() ).
    setApplyDefaultValue( m_defaultCompressStrength );
  
  this->template getWrapper< array1d< real64 > >( viewKeyStruct::deltaCoefficientString() ).
    setApplyDefaultValue( m_defaultDeltaCoefficient );
}

template< typename BASE >
void Damage< BASE >::allocateConstitutiveData( dataRepository::Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  m_newDamage.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldDamage.resize( 0, numConstitutivePointsPerParentIndex );
  m_damageGrad.resize( 0, numConstitutivePointsPerParentIndex, 3 );
  m_strainEnergyDensity.resize( 0, numConstitutivePointsPerParentIndex );
  m_volStrain.resize( 0, numConstitutivePointsPerParentIndex );
  m_extDrivingForce.resize( 0, numConstitutivePointsPerParentIndex );
  m_biotCoefficient.resize( parent.size() );
  m_criticalFractureEnergy.resize( parent.size() );
  m_tensileStrength.resize( parent.size() );
  m_compressStrength.resize( parent.size() );
  m_deltaCoefficient.resize( parent.size() );
  BASE::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

template< typename BASE >
void Damage< BASE >::saveConvergedState() const
{
  SolidBase::saveConvergedState(); // TODO: not ideal, as we have separate loops for base and derived data

  localIndex const numE = SolidBase::numElem();
  localIndex const numQ = SolidBase::numQuad();

  arrayView2d< real64 const > newDamage = m_newDamage;
  arrayView2d< real64 > oldDamage = m_oldDamage;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldDamage( k, q ) = newDamage( k, q );
    }
  } );
}

typedef Damage< ElasticIsotropic > DamageElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DamageElasticIsotropic, string const &, Group * const )

}
} /* namespace geos */
