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
 * @file KValueFlashModel.cpp
 */

#include "KValueFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/KValueFlashParameters.hpp"

#include "functions/FunctionManager.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

// Naming conventions
namespace
{
template< integer NUM_PHASE >
struct KValueFlashName
{
  static constexpr char const * name = "KValue";
};

template<>
struct KValueFlashName< 3 >
{
  static constexpr char const * name = "ThreePhaseKValue";
};
}

template< integer NUM_PHASE >
string KValueFlashModel< NUM_PHASE >::catalogName()
{
  return KValueFlashName< NUM_PHASE >::name;
}

template< integer NUM_PHASE >
KValueFlashModel< NUM_PHASE >::KValueFlashModel( string const & name,
                                                 ComponentProperties const & componentProperties,
                                                 ModelParameters const & modelParameters ):
  FunctionBase( name, componentProperties )
{
  m_parameters = modelParameters.get< KValueFlashParameters< NUM_PHASE > >();
  integer const numComps = m_componentProperties.getNumberOfComponents();
  m_presentComponents.resize( numComps );
  for( integer ic =0; ic < numComps; ++ic )
  {
    m_presentComponents[ic] = ic;
  }
}

template< integer NUM_PHASE >
typename KValueFlashModel< NUM_PHASE >::KernelWrapper
KValueFlashModel< NUM_PHASE >::createKernelWrapper() const
{
  string pressureTableName;
  string temperatureTableName;
  m_parameters->createTables( functionName(),
                              pressureTableName,
                              temperatureTableName );

  FunctionManager const & functionManager = FunctionManager::getInstance();
  TableFunction const * pressureTable = functionManager.getGroupPointer< TableFunction >( pressureTableName );
  TableFunction const * temperatureTable = functionManager.getGroupPointer< TableFunction >( temperatureTableName );

  return KernelWrapper( m_componentProperties.getNumberOfComponents(),
                        *pressureTable,
                        *temperatureTable,
                        m_presentComponents,
                        m_parameters->m_kValueHyperCube );
}

template< integer NUM_PHASE >
KValueFlashModelUpdate< NUM_PHASE >::KValueFlashModelUpdate( integer const numComponents,
                                                             TableFunction const & pressureTable,
                                                             TableFunction const & temperatureTable,
                                                             arrayView1d< integer const > const & presentComponents,
                                                             arrayView4d< real64 const > const & kValues ):
  m_numComponents( numComponents ),
  m_numPressurePoints( pressureTable.getCoordinates()[0].size()),
  m_numTemperaturePoints( temperatureTable.getCoordinates()[0].size()),
  m_pressureTable( pressureTable.createKernelWrapper()),
  m_temperatureTable( temperatureTable.createKernelWrapper()),
  m_presentComponents( presentComponents ),
  m_kValues( kValues )
{}

template< integer NUM_PHASE >
std::unique_ptr< ModelParameters >
KValueFlashModel< NUM_PHASE >::createParameters( std::unique_ptr< ModelParameters > parameters )
{
  return KValueFlashParameters< NUM_PHASE >::create( std::move( parameters ) );
}

// Instantiate
template class KValueFlashModel< 2 >;
template class KValueFlashModel< 3 >;

} // end namespace compositional

} // namespace constitutive

} // end namespace geos
