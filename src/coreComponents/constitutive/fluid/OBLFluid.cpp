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
 * @file OBLFluid.cpp
 */

#include "OBLFluid.hpp"
#include "functions/FunctionManager.hpp"
#include "dataRepository/InputFlags.hpp"
#include "dataRepository/RestartFlags.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

MultivariableTableFunction const *
makeOBLOperatorsTable( string const & OBLOperatorsTableFile,
                       string const & OBLFluidName,
                       FunctionManager & functionManager )
{
  string const tableName = OBLFluidName + "OBL_table";
  if( functionManager.hasGroup< MultivariableTableFunction >( tableName ) )
  {
    return functionManager.getGroupPointer< MultivariableTableFunction >( tableName );
  }
  else
  {
    MultivariableTableFunction * const table = dynamicCast< MultivariableTableFunction * >( functionManager.createChild( "MultivariableTableFunction", tableName ) );
    table->initializeFunctionFromFile ( OBLOperatorsTableFile );
    return table;
  }
}

template< typename INDEX_T = __uint128_t >
PythonFunction< INDEX_T > *
makePythonFunction( string const & OBLFluidName, FunctionManager & functionManager )
{
  string const pythonFunctionName = OBLFluidName + "PythonFunction";
  if( functionManager.hasGroup< PythonFunction< INDEX_T > >( pythonFunctionName ))
  {
    return functionManager.getGroupPointer< PythonFunction< INDEX_T > >( pythonFunctionName );
  }
  else
  {
    PythonFunction< INDEX_T > * function = dynamicCast< PythonFunction< INDEX_T > * >( functionManager.createChild( "PythonFunction", pythonFunctionName ));
    return function;
  }
}

OBLFluid::OBLFluid( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent ),
  m_OBLOperatorsTable( nullptr ),
  m_pythonFunction( nullptr )
{
  this->registerWrapper( viewKeyStruct::interpolatorModeString(), &m_interpolatorModeString ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "OBL interpolator mode: static or adaptive" );

  this->registerWrapper( viewKeyStruct::interpolatorTypeString(), &m_interpolatorTypeString ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "OBL interpolator type: multilinear or linear" );

  this->registerWrapper( viewKeyStruct::oblOperatorsTableFileString(), &m_OBLOperatorsTableFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "File containing OBL operator values for static mode interpolation" );
}

void OBLFluid::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  // set interpolator mode
  GEOS_THROW_IF( m_interpolatorModeString.empty() ||
                 (m_interpolatorModeString != "static" &&
                  m_interpolatorModeString != "adaptive"),
                 GEOS_FMT( "{}: Invalid interpolator mode: {}",
                           getFullName(),
                           m_interpolatorModeString ),
                 InputError );

  if( m_interpolatorModeString == "static" )
    m_interpolatorMode = OBLInterpolatorMode::Static;
  else if( m_interpolatorModeString == "adaptive" )
    m_interpolatorMode = OBLInterpolatorMode::Adaptive;

  // set interpolator type
  /*GEOS_WARNING_IF( m_interpolatorTypeString.empty() ||
                  (m_interpolatorTypeString != "multilinear" &&
                   m_interpolatorTypeString != "linear"),
                  GEOS_FMT( "{}: Invalid interpolator type, using multilinear interpolator",
                            getFullName() ) );*/

  GEOS_WARNING_IF( m_interpolatorTypeString == "linear",
                   GEOS_FMT( "{}: Linear interpolator type is not supported yet, using multilinear interpolator",
                             getFullName()) );

  if( !m_interpolatorTypeString.empty())
  {
    if( m_interpolatorTypeString == "multilinear" )
      m_interpolatorType = OBLInterpolatorType::Multilinear;
    else if( m_interpolatorTypeString == "linear" )
      m_interpolatorType = OBLInterpolatorType::Multilinear;
    else
      m_interpolatorType = OBLInterpolatorType::Multilinear;
  }
  else
    m_interpolatorType = OBLInterpolatorType::Multilinear;

  // set table file
  GEOS_THROW_IF( m_OBLOperatorsTableFile.empty() &&
                 m_interpolatorMode == OBLInterpolatorMode::Static,
                 GEOS_FMT( "{}: Invalid operator table file",
                           getFullName() ),
                 InputError );


  if( m_interpolatorMode == OBLInterpolatorMode::Static )
    m_OBLOperatorsTable = makeOBLOperatorsTable( m_OBLOperatorsTableFile, getName(), FunctionManager::getInstance());
  else if( m_interpolatorMode == OBLInterpolatorMode::Adaptive )
    m_pythonFunction = makePythonFunction< longIndex >( getName(), FunctionManager::getInstance());

  // raise intialization flag
  m_isInitialized = true;
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, OBLFluid, string const &, Group * const )

} // namespace constitutive

} // namespace geos
