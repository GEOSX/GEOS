/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalDensity.cpp
 */

#include "CompositionalDensity.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

CompositionalDensity::CompositionalDensity( string const & name,
                                            ComponentProperties const & componentProperties,
                                            integer const phaseIndex,
                                            ModelParameters const & modelParameters )
  : FunctionBase( name, componentProperties )
{
  EquationOfState const * equationOfState = modelParameters.get< EquationOfState >();
  string const eosName = equationOfState->m_equationsOfStateNames[phaseIndex];
  m_equationOfState = EnumStrings< EquationOfStateType >::fromString( eosName );

  Parameters const * densityParameters = modelParameters.get< Parameters >();

  // Calculate the dimensional volume shift
  m_componentDimensionalVolumeShift.resize( componentProperties.getNumberOfComponents());
  calculateDimensionalVolumeShift( componentProperties,
                                   m_equationOfState,
                                   densityParameters->m_componentVolumeShift.toSliceConst(),
                                   m_componentDimensionalVolumeShift );
}

CompositionalDensity::KernelWrapper
CompositionalDensity::createKernelWrapper() const
{
  return KernelWrapper( m_componentDimensionalVolumeShift, m_equationOfState );
}

std::unique_ptr< ModelParameters >
CompositionalDensity::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  auto params = EquationOfState::create( std::move( parameters ) );
  return Parameters::create( std::move( params ) );
}

void CompositionalDensity::calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                                            EquationOfStateType const & equationOfState,
                                                            arraySlice1d< real64 const > componentVolumeShift,
                                                            arraySlice1d< real64 > componentDimensionalVolumeShift )
{
  if( equationOfState == EquationOfStateType::PengRobinson )
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::calculateDimensionalVolumeShift( componentProperties,
                                                                            componentVolumeShift,
                                                                            componentDimensionalVolumeShift );
  }
  else if( equationOfState == EquationOfStateType::SoaveRedlichKwong )
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::calculateDimensionalVolumeShift( componentProperties,
                                                                                 componentVolumeShift,
                                                                                 componentDimensionalVolumeShift );
  }
}

// Compositional density parameters

CompositionalDensity::Parameters::Parameters( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

std::unique_ptr< ModelParameters >
CompositionalDensity::Parameters::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< Parameters >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< Parameters >( std::move( parameters ) );
}

void CompositionalDensity::Parameters::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::componentVolumeShiftString(), &m_componentVolumeShift ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Component volume shifts" );
}

void CompositionalDensity::Parameters::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                                    ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( componentProperties );

  integer const numComponents = fluid->numFluidComponents();

  if( m_componentVolumeShift.empty() )
  {
    m_componentVolumeShift.resize( numComponents );
    m_componentVolumeShift.zero();
  }
  GEOS_THROW_IF_NE_MSG( m_componentVolumeShift.size(), numComponents,
                        GEOS_FMT( "{}: invalid number of values in attribute '{}'", fluid->getFullName(),
                                  viewKeyStruct::componentVolumeShiftString() ),
                        InputError );
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
