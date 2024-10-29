/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file Gas.cpp
 */

#include "Gas.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

Gas::Gas( string const & name, Group * const parent ):
  ContinuumBase( name, parent ),
  m_referencePressure(),
  m_referenceTemperature( 300.0 ),
  m_jacobian( 0 ),
  m_temperature( 0 )
{
  // register default values
  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setApplyDefaultValue( m_referenceTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference temperature" );

  // register fields
  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::LEVEL_0).
    setDescription( "Jacobian" );

  registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setApplyDefaultValue( m_referenceTemperature ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Temperature" );
}


Gas::~Gas()
{}


void Gas::allocateConstitutiveData( dataRepository::Group & parent,
                                    localIndex const numConstitutivePointsPerParentIndex )
{
  ContinuumBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_jacobian.resize( 0 );
  m_temperature.resize( 0 );
}


void Gas::postInputInitialization()
{
  ContinuumBase::postInputInitialization();

  GEOS_THROW_IF( m_referencePressure < 0.0, "Reference pressure must be greater than 0.0", InputError );
  GEOS_THROW_IF( m_referenceTemperature < 0.0, "Reference temperature must be greater than 0.0", InputError );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, Gas, std::string const &, Group * const )
}
} /* namespace geos */
