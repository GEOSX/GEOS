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
 * @file CompositionalMultiphaseFluid.cpp
 */

#include "CompositionalMultiphaseFluid.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geos
{

namespace fields
{
namespace multifluid
{
DECLARE_FIELD( kValues,
               "kValues",
               array4dLayoutPhaseComp,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase equilibrium ratios" );
}
}

namespace constitutive
{

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::
CompositionalMultiphaseFluid( string const & name, Group * const parent )
  : MultiFluidBase( name, parent ),
  m_componentProperties( std::make_unique< compositional::ComponentProperties >( m_componentNames, m_componentMolarWeight ) ),
  m_parameters( createModelParameters() )
{
  using InputFlags = dataRepository::InputFlags;

  getWrapperBase( viewKeyStruct::componentNamesString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::componentCriticalPressureString(), &m_componentProperties->m_componentCriticalPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical pressures" );

  registerWrapper( viewKeyStruct::componentCriticalTemperatureString(), &m_componentProperties->m_componentCriticalTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical temperatures" );

  registerWrapper( viewKeyStruct::componentAcentricFactorString(), &m_componentProperties->m_componentAcentricFactor ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component acentric factors" );

  registerWrapper( viewKeyStruct::componentVolumeShiftString(), &m_componentProperties->m_componentVolumeShift ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component volume shifts" );

  registerWrapper( viewKeyStruct::componentBinaryCoeffString(), &m_componentProperties->m_componentBinaryCoeff ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Table of binary interaction coefficients" );

  registerField( fields::multifluid::kValues{}, &m_kValues );

  // Link parameters specific to each model
  m_parameters->registerParameters( this );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
integer CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] = { "water" };
  return PVTProps::PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
string CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::catalogName()
{
  return GEOS_FMT( "Compositional{}Fluid{}",
                   FLASH::catalogName(),
                   PHASE1::Viscosity::catalogName() );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::allocateConstitutiveData( dataRepository::Group & parent,
                                                                                              localIndex const numConstitutivePointsPerParentIndex )
{
  MultiFluidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  // Zero k-Values to force initialisation with Wilson k-Values
  m_kValues.zero();
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::postInputInitialization()
{
  MultiFluidBase::postInputInitialization();

  integer const NC = numFluidComponents();
  integer const NP = numFluidPhases();

  GEOS_THROW_IF_NE_MSG( NP, NUM_PHASES,
                        GEOS_FMT( "{}: invalid number of phases in '{}'. There should be {} phases",
                                  getFullName(), viewKeyStruct::phaseNamesString(), NUM_PHASES ),
                        InputError );

  auto const checkInputSize = [&]( auto const & array, integer const expected, string const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), expected,
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );

  };
  checkInputSize( m_componentProperties->m_componentCriticalPressure, NC, viewKeyStruct::componentCriticalPressureString() );
  checkInputSize( m_componentProperties->m_componentCriticalTemperature, NC, viewKeyStruct::componentCriticalTemperatureString() );
  checkInputSize( m_componentProperties->m_componentAcentricFactor, NC, viewKeyStruct::componentAcentricFactorString() );

  if( m_componentProperties->m_componentVolumeShift.empty() )
  {
    m_componentProperties->m_componentVolumeShift.resize( NC );
    m_componentProperties->m_componentVolumeShift.zero();
  }
  checkInputSize( m_componentProperties->m_componentVolumeShift, NC, viewKeyStruct::componentVolumeShiftString() );

  array2d< real64 > & componentBinaryCoeff = m_componentProperties->m_componentBinaryCoeff;
  if( componentBinaryCoeff.empty() )
  {
    componentBinaryCoeff.resize( NC, NC );
    componentBinaryCoeff.zero();
  }
  checkInputSize( componentBinaryCoeff, NC * NC, viewKeyStruct::componentBinaryCoeffString() );

  // Binary interaction coefficients should be symmetric and have zero diagonal
  GEOS_THROW_IF_NE_MSG( componentBinaryCoeff.size( 0 ), NC,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::componentBinaryCoeffString() ),
                        InputError );
  GEOS_THROW_IF_NE_MSG( componentBinaryCoeff.size( 1 ), NC,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), viewKeyStruct::componentBinaryCoeffString() ),
                        InputError );
  for( integer ic = 0; ic < NC; ++ic )
  {
    GEOS_THROW_IF_GT_MSG( LvArray::math::abs( componentBinaryCoeff( ic, ic )), MultiFluidConstants::epsilon,
                          GEOS_FMT( "{}: {} entry at ({},{}) is {}: should be zero", getFullName(), viewKeyStruct::componentBinaryCoeffString(),
                                    ic, ic, componentBinaryCoeff( ic, ic ) ),
                          InputError );
    for( integer jc = ic + 1; jc < NC; ++jc )
    {
      real64 const difference = LvArray::math::abs( componentBinaryCoeff( ic, jc )-componentBinaryCoeff( jc, ic ));
      GEOS_THROW_IF_GT_MSG( difference, MultiFluidConstants::epsilon,
                            GEOS_FMT( "{}: {} entry at ({},{}) is {} and is different from entry at ({},{}) which is {}",
                                      getFullName(), viewKeyStruct::componentBinaryCoeffString(),
                                      ic, jc, componentBinaryCoeff( ic, jc ), jc, ic, componentBinaryCoeff( jc, ic ) ),
                            InputError );
    }
  }

  m_parameters->postInputInitialization( this, *m_componentProperties );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();

  // Create the fluid models
  createModels();
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::resizeFields( localIndex const size, localIndex const numPts )
{
  MultiFluidBase::resizeFields( size, numPts );

  m_kValues.resize( size, numPts, numFluidPhases()-1, numFluidComponents() );
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::deliverClone( string const & name,
                                                                             Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  return clone;
}

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
typename CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::KernelWrapper
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::createKernelWrapper()
{
  return KernelWrapper( *m_componentProperties,
                        *m_flash,
                        *m_phase1,
                        *m_phase2,
                        *m_phase3,
                        m_componentMolarWeight,
                        m_useMass,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView(),
                        m_kValues.toView() );
}

// Create the fluid models
template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
void CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::createModels()
{
  m_flash = std::make_unique< FLASH >( getName() + '_' + FLASH::catalogName(),
                                       *m_componentProperties,
                                       *m_parameters );

  m_phase1 = std::make_unique< PHASE1 >( GEOS_FMT( "{}_PhaseModel1", getName() ),
                                         *m_componentProperties,
                                         0,
                                         *m_parameters );

  m_phase2 = std::make_unique< PHASE2 >( GEOS_FMT( "{}_PhaseModel2", getName() ),
                                         *m_componentProperties,
                                         1,
                                         *m_parameters );

  m_phase3 = std::make_unique< PHASE3 >( GEOS_FMT( "{}_PhaseModel3", getName() ),
                                         *m_componentProperties,
                                         2,
                                         *m_parameters );
}

// Create the fluid models
template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
std::unique_ptr< compositional::ModelParameters >
CompositionalMultiphaseFluid< FLASH, PHASE1, PHASE2, PHASE3 >::createModelParameters()
{
  std::unique_ptr< compositional::ModelParameters > parameters;
  parameters = FLASH::createParameters( std::move( parameters ));
  parameters = PHASE1::createParameters( std::move( parameters ));
  parameters = PHASE2::createParameters( std::move( parameters ));
  parameters = PHASE3::createParameters( std::move( parameters ));
  return parameters;
}

// Explicit instantiation of the model template.
template class CompositionalMultiphaseFluid<
    compositional::NegativeTwoPhaseFlashModel,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::ConstantViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::ConstantViscosity, compositional::NullModel > >;
template class CompositionalMultiphaseFluid<
    compositional::NegativeTwoPhaseFlashModel,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel >,
    compositional::PhaseModel< compositional::CompositionalDensity, compositional::LohrenzBrayClarkViscosity, compositional::NullModel > >;

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalTwoPhaseConstantViscosity,
                        string const &,
                        dataRepository::Group * const )

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalTwoPhaseLohrenzBrayClarkViscosity,
                        string const &,
                        dataRepository::Group * const )

} // namespace constitutive

} // namespace geos
