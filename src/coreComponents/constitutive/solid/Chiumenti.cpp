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
 *  @file Chiumenti.cpp
 */

#include "Chiumenti.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

Chiumenti::Chiumenti( string const & name, Group * const parent ):
  HyperelasticMMS( name, parent ),
  m_damage(),
  m_lengthScale(),
  m_strengthScale(),
  m_criticalLength(),
  m_failureStrength(),
  m_energyReleaseRate()
{
  // register default values
  registerWrapper( viewKeyStruct::criticalLengthString(), &m_criticalLength ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Critical length" );

  registerWrapper( viewKeyStruct::failureStrengthString(), &m_failureStrength ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Failure strength" );

  registerWrapper( viewKeyStruct::energyReleaseRateString(), &m_energyReleaseRate ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Energy release rate" );

  // register fields
  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::lengthScaleString(), &m_lengthScale ).
    setApplyDefaultValue( DBL_MIN ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "Array of quadrature point damage values" );

  registerWrapper( viewKeyStruct::strengthScaleString(), &m_strengthScale ).
    setApplyDefaultValue( 1.0 ).
    setPlotLevel( PlotLevel::LEVEL_0).
    setDescription( "Strength scale" );
}


Chiumenti::~Chiumenti()
{}


void Chiumenti::allocateConstitutiveData( dataRepository::Group & parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  HyperelasticMMS::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
  m_lengthScale.resize( 0 );
  m_strengthScale.resize( 0 );
}


void Chiumenti::postInputInitialization()
{
  HyperelasticMMS::postInputInitialization();

  // CC: TODO double check model inputs
  GEOS_THROW_IF( m_criticalLength < 0.0, "Critical length must be a positive number.", InputError );
  GEOS_THROW_IF( m_failureStrength < 0.0, "Failure strength must be a positive number.", InputError );
  GEOS_THROW_IF( m_energyReleaseRate < 0.0, "Energy release rate must be a positive number.", InputError );
}


void Chiumenti::saveConvergedState() const
{
  HyperelasticMMS::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, Chiumenti, std::string const &, Group * const )
}
} /* namespace geos */
