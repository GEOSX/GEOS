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
 * @file EquationOfState.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_EQUATIONOFSTATE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_EQUATIONOFSTATE_HPP_

#include "ModelParameters.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "dataRepository/InputFlags.hpp"
#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

enum class EquationOfStateType : integer
{
  PengRobinson,
  SoaveRedlichKwong
};

ENUM_STRINGS( EquationOfStateType,
              "pr",
              "srk" );

class EquationOfState : public ModelParameters
{
public:
  EquationOfState( std::unique_ptr< ModelParameters > parameters ):
    ModelParameters( std::move( parameters ) )
  {}

  ~EquationOfState() override = default;

  static std::unique_ptr< ModelParameters > create( std::unique_ptr< ModelParameters > parameters )
  {
    if( parameters && parameters->get< EquationOfState >() != nullptr )
    {
      return parameters;
    }
    return std::make_unique< EquationOfState >( std::move( parameters ) );
  }

  string_array m_equationsOfStateNames;

protected:
  void registerParametersImpl( MultiFluidBase * fluid ) override
  {
    fluid->registerWrapper( viewKeyStruct::equationsOfStateString(), &m_equationsOfStateNames ).
      setInputFlag( dataRepository::InputFlags::REQUIRED ).
      setDescription( "List of equation of state types for each phase. Valid options:\n* " +
                      EnumStrings< EquationOfStateType >::concat( "\n* " ) );
  }

  void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override
  {
    GEOS_UNUSED_VAR( componentProperties );

    integer const numPhase = fluid->numFluidPhases();

    GEOS_THROW_IF_NE_MSG( m_equationsOfStateNames.size(), numPhase,
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", fluid->getFullName(),
                                    viewKeyStruct::equationsOfStateString() ),
                          InputError );

    // If any value is invalid conversion will throw
    for( string const & eos : m_equationsOfStateNames )
    {
      EnumStrings< EquationOfStateType >::fromString( eos );
    }
  }

  struct viewKeyStruct
  {
    static constexpr char const * equationsOfStateString() { return "equationsOfState"; }
  };
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_EQUATIONOFSTATE_HPP_
