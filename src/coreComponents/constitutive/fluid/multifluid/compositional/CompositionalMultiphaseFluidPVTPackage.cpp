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
 * @file CompositionalMultiphaseFluidPVTPackage.cpp
 */

#include "CompositionalMultiphaseFluidPVTPackage.hpp"

#include "codingUtilities/Utilities.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"

#include "pvt/pvt.hpp"

#include <map>
#include <utility>

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CompositionalMultiphaseFluidPVTPackage::CompositionalMultiphaseFluidPVTPackage( string const & name, Group * const parent )
  : MultiFluidBase( name, parent )
{
  getWrapperBase( viewKeyStruct::componentNamesString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::componentMolarWeightString() ).setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString() ).setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::equationsOfStateString(), &m_equationsOfState ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of equation of state types for each phase" );

  registerWrapper( viewKeyStruct::constantPhaseViscosityString(), &m_constantPhaseViscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Viscosity for each phase" );

  registerWrapper( viewKeyStruct::componentCriticalPressureString(), &m_componentCriticalPressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical pressures" );

  registerWrapper( viewKeyStruct::componentCriticalTemperatureString(), &m_componentCriticalTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component critical temperatures" );

  registerWrapper( viewKeyStruct::componentAcentricFactorString(), &m_componentAcentricFactor ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Component acentric factors" );

  registerWrapper( viewKeyStruct::componentVolumeShiftString(), &m_componentVolumeShift ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component volume shifts" );

  registerWrapper( viewKeyStruct::componentBinaryCoeffString(), &m_componentBinaryCoeff ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Table of binary interaction coefficients" );
}

integer CompositionalMultiphaseFluidPVTPackage::getWaterPhaseIndex() const
{
  string const expectedWaterPhaseNames[] = { "water" };
  return PVTProps::PVTFunctionHelpers::findName( m_phaseNames, expectedWaterPhaseNames, viewKeyStruct::phaseNamesString() );
}

void CompositionalMultiphaseFluidPVTPackage::postInputInitialization()
{
  MultiFluidBase::postInputInitialization();

  auto const getPVTPackagePhaseType = [&]( string const & phaseName )
  {
    static map< string, pvt::PHASE_TYPE > const phaseTypes
    {
      { "gas", pvt::PHASE_TYPE::GAS },
      { "oil", pvt::PHASE_TYPE::OIL },
      { "water", pvt::PHASE_TYPE::LIQUID_WATER_RICH }
    };
    return findOption( phaseTypes, phaseName, viewKeyStruct::phaseNamesString(), getFullName() );
  };

  m_phaseTypes.resize( numFluidPhases() );
  std::transform( m_phaseNames.begin(), m_phaseNames.end(), m_phaseTypes.begin(), getPVTPackagePhaseType );

  integer const NC = numFluidComponents();
  integer const NP = numFluidPhases();

  auto const checkInputSize = [&]( auto const & array, integer const expected, string const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), expected,
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );

  };
  checkInputSize( m_equationsOfState, NP, viewKeyStruct::equationsOfStateString() );
  checkInputSize( m_componentCriticalPressure, NC, viewKeyStruct::componentCriticalPressureString() );
  checkInputSize( m_componentCriticalTemperature, NC, viewKeyStruct::componentCriticalTemperatureString() );
  checkInputSize( m_componentAcentricFactor, NC, viewKeyStruct::componentAcentricFactorString() );

  if( m_constantPhaseViscosity.empty() )
  {
    m_constantPhaseViscosity.resize( NP );
    for( integer ip = 0; ip < NP; ++ip )
    {
      m_constantPhaseViscosity[ip] = 0.001;      // Default value = 1 cP
    }
  }
  checkInputSize( m_constantPhaseViscosity, NP, viewKeyStruct::constantPhaseViscosityString() );

  if( m_componentVolumeShift.empty() )
  {
    m_componentVolumeShift.resize( NC );
    m_componentVolumeShift.zero();
  }
  checkInputSize( m_componentVolumeShift, NC, viewKeyStruct::componentVolumeShiftString() );

  if( m_componentBinaryCoeff.empty() )
  {
    m_componentBinaryCoeff.resize( NC, NC );
    m_componentBinaryCoeff.zero();
  }
  checkInputSize( m_componentBinaryCoeff, NC * NC, viewKeyStruct::componentBinaryCoeffString() );
}

void CompositionalMultiphaseFluidPVTPackage::initializePostSubGroups()
{
  MultiFluidBase::initializePostSubGroups();
  createFluid();
}

void CompositionalMultiphaseFluidPVTPackage::createFluid()
{
  auto const getCompositionalEosType = [&]( string const & name )
  {
    static map< string, pvt::EOS_TYPE > const eosTypes =
    {
      { "PR", pvt::EOS_TYPE::PENG_ROBINSON },
      { "SRK", pvt::EOS_TYPE::REDLICH_KWONG_SOAVE }
    };
    return findOption( eosTypes, name, viewKeyStruct::equationsOfStateString(), getFullName() );
  };

  std::vector< pvt::EOS_TYPE > eos( numFluidPhases() );
  std::transform( m_equationsOfState.begin(), m_equationsOfState.end(), eos.begin(), getCompositionalEosType );

  std::vector< pvt::PHASE_TYPE > phases( m_phaseTypes.begin(), m_phaseTypes.end() );
  std::vector< string > const components( m_componentNames.begin(), m_componentNames.end() );
  std::vector< double > const Mw( m_componentMolarWeight.begin(), m_componentMolarWeight.end() );
  std::vector< double > const Tc( m_componentCriticalTemperature.begin(), m_componentCriticalTemperature.end() );
  std::vector< double > const Pc( m_componentCriticalPressure.begin(), m_componentCriticalPressure.end() );
  std::vector< double > const Omega( m_componentAcentricFactor.begin(), m_componentAcentricFactor.end() );

  m_fluid = pvt::MultiphaseSystemBuilder::buildCompositional( pvt::COMPOSITIONAL_FLASH_TYPE::NEGATIVE_OIL_GAS, phases, eos,
                                                              components, Mw, Tc, Pc, Omega );
}

std::unique_ptr< ConstitutiveBase >
CompositionalMultiphaseFluidPVTPackage::deliverClone( string const & name,
                                                      Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = MultiFluidBase::deliverClone( name, parent );
  CompositionalMultiphaseFluidPVTPackage & fluid = dynamicCast< CompositionalMultiphaseFluidPVTPackage & >( *clone );
  fluid.m_phaseTypes = m_phaseTypes;
  fluid.createFluid();
  return clone;
}

CompositionalMultiphaseFluidPVTPackage::KernelWrapper::
  KernelWrapper( pvt::MultiphaseSystem & fluid,
                 arrayView1d< pvt::PHASE_TYPE > const & phaseTypes,
                 arrayView1d< geos::real64 const > const & constantPhaseViscosity,
                 arrayView1d< geos::real64 const > const & componentMolarWeight,
                 bool useMass,
                 PhaseProp::ViewType phaseFraction,
                 PhaseProp::ViewType phaseDensity,
                 PhaseProp::ViewType phaseMassDensity,
                 PhaseProp::ViewType phaseViscosity,
                 PhaseProp::ViewType phaseEnthalpy,
                 PhaseProp::ViewType phaseInternalEnergy,
                 PhaseComp::ViewType phaseCompFraction,
                 FluidProp::ViewType totalDensity )
  : MultiFluidBase::KernelWrapper( componentMolarWeight,
                                   useMass,
                                   std::move( phaseFraction ),
                                   std::move( phaseDensity ),
                                   std::move( phaseMassDensity ),
                                   std::move( phaseViscosity ),
                                   std::move( phaseEnthalpy ),
                                   std::move( phaseInternalEnergy ),
                                   std::move( phaseCompFraction ),
                                   std::move( totalDensity ) ),
  m_fluid( fluid ),
  m_phaseTypes( phaseTypes ),
  m_constantPhaseViscosity( constantPhaseViscosity )
{}

CompositionalMultiphaseFluidPVTPackage::KernelWrapper
CompositionalMultiphaseFluidPVTPackage::createKernelWrapper()
{
  return KernelWrapper( *m_fluid,
                        m_phaseTypes,
                        m_constantPhaseViscosity,
                        m_componentMolarWeight,
                        m_useMass,
                        m_phaseFraction.toView(),
                        m_phaseDensity.toView(),
                        m_phaseMassDensity.toView(),
                        m_phaseViscosity.toView(),
                        m_phaseEnthalpy.toView(),
                        m_phaseInternalEnergy.toView(),
                        m_phaseCompFraction.toView(),
                        m_totalDensity.toView() );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompositionalMultiphaseFluidPVTPackage, string const &, Group * const )

} // namespace constitutive

} // namespace geos
