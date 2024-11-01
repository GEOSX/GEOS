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

  // Calculate the dimensional volume shift
  m_componentDimensionalVolumeShift.resize( componentProperties.getNumberOfComponents());
  calculateDimensionalVolumeShift( componentProperties,
                                   m_equationOfState,
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
  return EquationOfState::create( std::move( parameters ) );
}

void CompositionalDensity::calculateDimensionalVolumeShift( ComponentProperties const & componentProperties,
                                                            EquationOfStateType const & equationOfState,
                                                            arraySlice1d< real64 > componentDimensionalVolumeShift )
{
  if( equationOfState == EquationOfStateType::PengRobinson )
  {
    CubicEOSPhaseModel< PengRobinsonEOS >::calculateDimensionalVolumeShift( componentProperties,
                                                                            componentDimensionalVolumeShift );
  }
  else if( equationOfState == EquationOfStateType::SoaveRedlichKwong )
  {
    CubicEOSPhaseModel< SoaveRedlichKwongEOS >::calculateDimensionalVolumeShift( componentProperties,
                                                                                 componentDimensionalVolumeShift );
  }
}

} // namespace compositional

} // namespace constitutive

} // namespace geos
