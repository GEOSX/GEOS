/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file KValueFlashParameters.hpp
 */

#include "KValueFlashParameters.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"
#ifdef GEOS_USE_MATHPRESSO
#include "functions/SymbolicFunction.hpp"
#include "functions/CompositeFunction.hpp"
#endif

#include "common/Units.hpp"
#include "common/format/table/TableFormatter.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< integer NUM_PHASE >
KValueFlashParameters< NUM_PHASE >::KValueFlashParameters( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

template< integer NUM_PHASE >
std::unique_ptr< ModelParameters > KValueFlashParameters< NUM_PHASE >::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< KValueFlashParameters >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< KValueFlashParameters >( std::move( parameters ) );
}

template< integer NUM_PHASE >
void KValueFlashParameters< NUM_PHASE >::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::pressureCoordinatesString(), &m_pressureCoordinates ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "List of pressure values for interpolation of k-values." );

  fluid->registerWrapper( viewKeyStruct::temperatureCoordinatesString(), &m_temperatureCoordinates ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "List of temperature values for interpolation of k-values." );

  fluid->registerWrapper( viewKeyStruct::kValueTablesString(), &m_kValueTables ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "List of k-value tables for each phase." );
}

template< integer NUM_PHASE >
void KValueFlashParameters< NUM_PHASE >::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                                      ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( componentProperties );

  integer const numFluidPhase = fluid->numFluidPhases();
  integer const numFluidComponent = fluid->numFluidComponents();

  GEOS_THROW_IF_NE_MSG( numPhases, numFluidPhase,
                        GEOS_FMT( "{}: invalid number of phases for the fluid.", fluid->getFullName() ),
                        InputError );

  if( !m_pressureCoordinates.empty())
  {
    // If coordinates are provided then there must be at least 2
    GEOS_THROW_IF_LT_MSG( m_pressureCoordinates.size(), 2,
                          GEOS_FMT( "{}: invalid number of pressure coordinates provided in {}. "
                                    "At least 2 values must be provided", fluid->getFullName(), viewKeyStruct::pressureCoordinatesString() ),
                          InputError );

    // Values must be increasing
    GEOS_THROW_IF( isIncreasing( m_pressureCoordinates.toSliceConst()),
                   GEOS_FMT( "{}: invalid number of pressure coordinates provided in {}. "
                             "Values must be strictly increasing.", fluid->getFullName(), viewKeyStruct::pressureCoordinatesString() ),
                   InputError );
  }

  if( !m_temperatureCoordinates.empty())
  {
    // If coordinates are provided then there must be at least 2
    GEOS_THROW_IF_LT_MSG( m_temperatureCoordinates.size(), 2,
                          GEOS_FMT( "{}: invalid number of temperature coordinates provided in {}. "
                                    "At least 2 values must be provided", fluid->getFullName(), viewKeyStruct::temperatureCoordinatesString() ),
                          InputError );

    // Values must be increasing
    GEOS_THROW_IF( isIncreasing( m_temperatureCoordinates.toSliceConst()),
                   GEOS_FMT( "{}: invalid number of temperature coordinates provided in {}. "
                             "Values must be strictly increasing.", fluid->getFullName(), viewKeyStruct::temperatureCoordinatesString() ),
                   InputError );
  }

  integer const numTables = m_kValueTables.size();

  GEOS_THROW_IF_NE_MSG( numTables, (numPhases-1)*numFluidComponent,
                        GEOS_FMT( "{}: invalid number of k-value tables provided.", fluid->getFullName() ),
                        InputError );

  // Check that the tables exist and are 2D
  FunctionManager & functionManager = FunctionManager::getInstance();
  for( integer tableIndex = 0; tableIndex < numTables; ++tableIndex )
  {
    string const functionName = m_kValueTables[tableIndex];
    FunctionBase * function = functionManager.getGroupPointer< FunctionBase >( functionName );
    GEOS_THROW_IF( function == nullptr,
                   GEOS_FMT( "Function with name {} not found. ", functionName ),
                   InputError );

    function->initializeFunction();

    integer numDims = 0;
    if( TableFunction const * tableFunction = dynamicCast< TableFunction const * >( function ))
    {
      numDims = tableFunction->numDimensions();
    }
#ifdef GEOS_USE_MATHPRESSO
    else if( SymbolicFunction const * symbolicFunction = dynamicCast< SymbolicFunction const * >( function ))
    {
      numDims = symbolicFunction->getWrapper< string_array >( "variableNames" ).reference().size();
    }
    else if( CompositeFunction const * compositeFunction = dynamicCast< CompositeFunction const * >( function ))
    {
      numDims = compositeFunction->getWrapper< string_array >( "variableNames" ).reference().size();
    }
#endif
    GEOS_THROW_IF_NE_MSG( numDims, 2,
                          GEOS_FMT( "Function with name {} must have a dimension of 2. ", functionName ),
                          InputError );
  }

  generateHyperCube( numFluidComponent );
  validateKValues( fluid );
}

template< integer NUM_PHASE >
void KValueFlashParameters< NUM_PHASE >::createTables( string const & name,
                                                       string & pressureTableName,
                                                       string & temperatureTableName ) const
{
  FunctionManager & functionManager = FunctionManager::getInstance();
  TableFunction * tableFunction = nullptr;

  // Pressure lookup table
  pressureTableName = GEOS_FMT( "{}_pressure_index_table", name );
  if( !functionManager.hasGroup< TableFunction >( pressureTableName ) )
  {
    array1d< real64 > indices( m_pressureValues[0].size() );
    for( integer i = 0; i < m_pressureValues[0].size(); ++i )
    {
      indices[i] = static_cast< real64 >(i);
    }
    tableFunction = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", pressureTableName ) );
    tableFunction->setTableCoordinates( m_pressureValues, { units::Pressure } );
    tableFunction->setTableValues( indices, units::Dimensionless );
    tableFunction->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    tableFunction->initializeFunction();
  }

  // Temperature lookup table
  temperatureTableName = GEOS_FMT( "{}_temperature_index_table", name );
  if( !functionManager.hasGroup< TableFunction >( temperatureTableName ) )
  {
    array1d< real64 > indices( m_temperatureValues[0].size() );
    for( integer i = 0; i < m_temperatureValues[0].size(); ++i )
    {
      indices[i] = static_cast< real64 >(i);
    }
    tableFunction = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", temperatureTableName ) );
    tableFunction->setTableCoordinates( m_temperatureValues, { units::Temperature } );
    tableFunction->setTableValues( indices, units::Dimensionless );
    tableFunction->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    tableFunction->initializeFunction();
  }
}

template< integer NUM_PHASE >
void KValueFlashParameters< NUM_PHASE >::generateHyperCube( integer const numComps )
{
  FunctionManager & functionManager = FunctionManager::getInstance();
  TableFunction const * tableFunction = nullptr;

  real64 minPressure = LvArray::NumericLimits< real64 >::max;
  real64 maxPressure = LvArray::NumericLimits< real64 >::min;
  real64 minTemperature = LvArray::NumericLimits< real64 >::max;
  real64 maxTemperature = LvArray::NumericLimits< real64 >::min;
  integer numPressurePoints = 0;
  integer numTemperaturePoints = 0;

  // If the pressure coordinates or temperature coordinates are not provided, read the maximum
  // pressure and temperature intervals from the provided tables. These will be uniformly
  // discretised.
  if( m_pressureCoordinates.empty() || m_temperatureCoordinates.empty() )
  {
    for( integer phaseIndex = 0; phaseIndex < numPhases-1; ++phaseIndex )
    {
      for( integer compIndex = 0; compIndex < numComps; ++compIndex )
      {
        integer const tableIndex = numComps*phaseIndex + compIndex;
        string const tableName = m_kValueTables[tableIndex];
        tableFunction = functionManager.getGroupPointer< TableFunction >( tableName );

        if( tableFunction == nullptr )
        {
          continue;
        }

        ArrayOfArraysView< real64 const > coordinates = tableFunction->getCoordinates();

        integer const np = coordinates[0].size();
        numPressurePoints = LvArray::math::max( numPressurePoints, np );
        minPressure = LvArray::math::min( minPressure, coordinates[0][0] );
        maxPressure = LvArray::math::max( maxPressure, coordinates[0][np-1] );

        integer const nt = coordinates[1].size();
        numTemperaturePoints = LvArray::math::max( numTemperaturePoints, nt );
        minTemperature = LvArray::math::min( minTemperature, coordinates[1][0] );
        maxTemperature = LvArray::math::max( maxTemperature, coordinates[1][nt-1] );
      }
    }

    if( numPressurePoints == 1 )
    {
      numPressurePoints = 2;
      maxPressure = minPressure + 1.0e7;
    }
    if( numTemperaturePoints == 1 )
    {
      numTemperaturePoints = 2;
      maxTemperature = minTemperature + 10.0;
    }
  }

  if( !m_pressureCoordinates.empty() )
  {
    numPressurePoints = m_pressureCoordinates.size();
  }
  if( !m_temperatureCoordinates.empty() )
  {
    numTemperaturePoints = m_temperatureCoordinates.size();
  }

  GEOS_THROW_IF_EQ_MSG( numPressurePoints, 0,
                        GEOS_FMT( "Failed to calculate number of pressure points for k-value interpolation. "
                                  "Provide values for {}.", KValueFlashParameters::viewKeyStruct::pressureCoordinatesString() ),
                        InputError );
  GEOS_THROW_IF_EQ_MSG( numTemperaturePoints, 0,
                        GEOS_FMT( "Failed to calculate number of temperature points for k-value interpolation. "
                                  "Provide values for {}.", KValueFlashParameters::viewKeyStruct::temperatureCoordinatesString() ),
                        InputError );

  // Create the pressure index lookup table
  m_pressureValues.resize( 1 );
  if( m_pressureCoordinates.empty() )
  {
    real64 const dp = (maxPressure - minPressure) / (numPressurePoints-1);
    m_pressureValues[0].resize( numPressurePoints );
    for( localIndex i = 0; i < numPressurePoints; ++i )
    {
      m_pressureValues[0][i] = minPressure + i * dp;
    }
  }
  else
  {
    m_pressureValues[0].resize( numPressurePoints );
    for( localIndex i = 0; i < numPressurePoints; ++i )
    {
      m_pressureValues[0][i] = m_pressureCoordinates[i];
    }
  }

  // Create the temperature index lookup table
  m_temperatureValues.resize( 1 );
  if( m_temperatureCoordinates.empty() )
  {
    real64 const dt = (maxTemperature - minTemperature) / (numTemperaturePoints-1);
    m_temperatureValues[0].resize( numTemperaturePoints );
    for( localIndex i = 0; i < numTemperaturePoints; ++i )
    {
      m_temperatureValues[0][i] = minTemperature + i * dt;
    }
  }
  else
  {
    m_temperatureValues[0].resize( numTemperaturePoints );
    for( localIndex i = 0; i < numTemperaturePoints; ++i )
    {
      m_temperatureValues[0][i] = m_temperatureCoordinates[i];
    }
  }

  // Create the "hypercube" of k-values
  real64 lookupValue[2]{};
  m_kValueHyperCube.resize( numPhases-1, numComps, numPressurePoints, numTemperaturePoints );
  for( integer phaseIndex = 0; phaseIndex < numPhases-1; ++phaseIndex )
  {
    for( integer compIndex = 0; compIndex < numComps; ++compIndex )
    {
      integer const tableIndex = numComps*phaseIndex + compIndex;
      string const tableName = m_kValueTables[tableIndex];
      FunctionBase const * function = functionManager.getGroupPointer< FunctionBase >( tableName );
      for( integer pressureIndex = 0; pressureIndex < numPressurePoints; ++pressureIndex )
      {
        lookupValue[0] = m_pressureValues[0][pressureIndex];
        for( integer temperatureIndex = 0; temperatureIndex < numTemperaturePoints; ++temperatureIndex )
        {
          lookupValue[1] = m_temperatureValues[0][temperatureIndex];
          m_kValueHyperCube( phaseIndex, compIndex, pressureIndex, temperatureIndex ) = function->evaluate( lookupValue );
        }
      }
    }
  }
}

template< integer NUM_PHASE >
bool KValueFlashParameters< NUM_PHASE >::isIncreasing( arraySlice1d< real64 const > const & array )
{
  localIndex const size = array.size();
  GEOS_ASSERT( 2 <= size );
  bool increasing = true;
  real64 constexpr epsilon = MultiFluidConstants::epsilon;
  for( localIndex i = 1; i < size; ++i )
  {
    increasing = increasing && (epsilon < array[i] - array[i-1]);
  }
  return increasing;
}

template< integer NUM_PHASE >
bool KValueFlashParameters< NUM_PHASE >::validateKValues( MultiFluidBase const * fluid ) const
{
  integer const numComps = m_kValueHyperCube.size( 1 );
  integer const numPressures = m_kValueHyperCube.size( 2 );
  integer const numTemperatures = m_kValueHyperCube.size( 3 );

  arrayView1d< string const > const phaseNames = fluid->phaseNames();
  arrayView1d< string const > const componentNames = fluid->componentNames();

  bool hasAtLeastOneNegative = false;
  bool hasAtLeastOneOneSided = false;

  integer const numTableColumns = numComps+3;
  TableData tableData;
  std::vector< string > tableRow( numTableColumns );
  for( integer phaseIndex = 0; phaseIndex < numPhases-1; ++phaseIndex )
  {
    for( integer pressureIndex = 0; pressureIndex < numPressures; ++pressureIndex )
    {
      for( integer temperatureIndex = 0; temperatureIndex < numTemperatures; ++temperatureIndex )
      {
        bool hasNegative = false;
        bool allMoreThanUnity = true;
        bool allLessThanUnity = true;
        for( integer compIndex = 0; compIndex < numComps; ++compIndex )
        {
          real64 const kValue = m_kValueHyperCube( phaseIndex, compIndex, pressureIndex, temperatureIndex );
          hasNegative = hasNegative || (kValue < 0.0);
          allMoreThanUnity = allMoreThanUnity && (1.0 < kValue);
          allLessThanUnity = allLessThanUnity && (kValue < 1.0);
        }
        hasAtLeastOneNegative = hasAtLeastOneNegative || hasNegative;
        hasAtLeastOneOneSided = hasAtLeastOneOneSided || (allMoreThanUnity || allLessThanUnity);
        if( (allMoreThanUnity || allLessThanUnity || hasNegative) && tableData.getTableDataRows().size() < 5 )
        {
          tableRow[0] = phaseNames[phaseIndex+1];
          tableRow[1] = GEOS_FMT( "{0:.3e}", m_pressureValues[0][pressureIndex] );
          tableRow[2] = GEOS_FMT( "{0:.2f}", m_temperatureValues[0][temperatureIndex] );
          for( integer compIndex = 0; compIndex < numComps; ++compIndex )
          {
            real64 const kValue = m_kValueHyperCube( phaseIndex, compIndex, pressureIndex, temperatureIndex );
            tableRow[3+compIndex] = GEOS_FMT( "{0:.3e}", kValue );
          }
          tableData.addRow( tableRow );
        }
      }
    }
  }

  if( !tableData.getTableDataRows().empty())
  {
    std::vector< TableLayout::ColumnParam > columnParameters;
    columnParameters.emplace_back( "Phase", TableLayout::Alignment::left );
    columnParameters.emplace_back( "Pressure", TableLayout::Alignment::right );
    columnParameters.emplace_back( "Temperature", TableLayout::Alignment::right );
    for( integer compIndex = 0; compIndex < numComps; ++compIndex )
    {
      columnParameters.emplace_back( componentNames[compIndex], TableLayout::Alignment::right );
    }
    TableLayout const tableLayout( columnParameters );
    TableTextFormatter const tableText( tableLayout );

    string message;
    if( hasAtLeastOneNegative )
    {
      message += "The provided tables of k-values have values which are not positive. ";
    }
    if( hasAtLeastOneOneSided )
    {
      message += "The provided tables of k-values have a pressure and temperature at "
                 "which all k-values are greater than unity or all k-values are less than unity. ";
    }

    GEOS_WARNING( GEOS_FMT( "{}: {}\n{}",
                            fluid->getFullName(), message, tableText.toString( tableData ) ));

    GEOS_THROW_IF( hasAtLeastOneNegative, GEOS_FMT( "{}: negative k-value found. ", fluid->getFullName() ),
                   InputError );
  }

  return true;
}

// Instantiate
template class KValueFlashParameters< 2 >;
template class KValueFlashParameters< 3 >;

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
