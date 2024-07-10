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
 * @file NegativeTwoPhaseFlashModel.cpp
 */

#include "NegativeTwoPhaseFlashModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

// Naming conventions
string NegativeTwoPhaseFlashModel::catalogName()
{
  return "TwoPhase";
}

NegativeTwoPhaseFlashModel::NegativeTwoPhaseFlashModel( string const & name,
                                                        ComponentProperties const & componentProperties,
                                                        ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties )
{
  m_parameters = modelParameters.get< EquationOfState >();
}

NegativeTwoPhaseFlashModel::KernelWrapper
NegativeTwoPhaseFlashModel::createKernelWrapper() const
{
  constexpr integer liquidIndex = 0;
  constexpr integer vapourIndex = 1;
  EquationOfStateType const liquidEos =  EnumStrings< EquationOfStateType >::fromString( m_parameters->m_equationsOfStateNames[liquidIndex] );
  EquationOfStateType const vapourEos =  EnumStrings< EquationOfStateType >::fromString( m_parameters->m_equationsOfStateNames[vapourIndex] );
  return KernelWrapper( m_componentProperties.getNumberOfComponents(), liquidIndex, vapourIndex, liquidEos, vapourEos );
}

NegativeTwoPhaseFlashModelUpdate::NegativeTwoPhaseFlashModelUpdate(
  integer const numComponents,
  integer const liquidIndex,
  integer const vapourIndex,
  EquationOfStateType const liquidEos,
  EquationOfStateType const vapourEos ):
  m_numComponents( numComponents ),
  m_liquidIndex( liquidIndex ),
  m_vapourIndex( vapourIndex ),
  m_liquidEos( liquidEos ),
  m_vapourEos( vapourEos )
{}

std::unique_ptr< ModelParameters >
NegativeTwoPhaseFlashModel::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  return EquationOfState::create( std::move( parameters ) );
}

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
