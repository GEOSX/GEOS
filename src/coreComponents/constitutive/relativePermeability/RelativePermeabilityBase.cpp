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
 * @file RelativePermeabilityBase.cpp
 */

#include "RelativePermeabilityBase.hpp"
#include "RelativePermeabilityFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

RelativePermeabilityBase::RelativePermeabilityBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::phaseTypesString(), &m_phaseTypes ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseOrderString(), &m_phaseOrder ).
    setSizedFromParent( 0 );

  registerField( fields::relperm::phaseRelPerm{}, &m_phaseRelPerm );
  registerField( fields::relperm::dPhaseRelPerm_dPhaseVolFraction{}, &m_dPhaseRelPerm_dPhaseVolFrac );

  registerField( fields::relperm::phaseTrappedVolFraction{}, &m_phaseTrappedVolFrac );

  registerField( fields::relperm::phaseRelPerm_n{}, &m_phaseRelPerm_n );

}

void RelativePermeabilityBase::postInputInitialization()
{
  ConstitutiveBase::postInputInitialization();

  integer const numPhases = numFluidPhases();
  GEOS_THROW_IF_LT_MSG( numPhases, 2,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );
  GEOS_THROW_IF_GT_MSG( numPhases, MAX_NUM_PHASES,
                        GEOS_FMT( "{}: invalid number of phases", getFullName() ),
                        InputError );

  m_phaseTypes.resize( numPhases );
  m_phaseOrder.resizeDefault( MAX_NUM_PHASES, -1 );

  auto const toPhaseType = [&]( string const & lookup )
  {
    static unordered_map< string, integer > const phaseDict =
    {
      { "gas", PhaseType::GAS },
      { "oil", PhaseType::OIL },
      { "water", PhaseType::WATER }
    };
    return findOption( phaseDict, lookup, viewKeyStruct::phaseNamesString(), getFullName() );
  };

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    m_phaseTypes[ip] = toPhaseType( m_phaseNames[ip] );
    m_phaseOrder[m_phaseTypes[ip]] = ip;
  }

  // call to correctly set member array tertiary sizes on the 'main' material object
  resizeFields( 0, 0 );

  // set labels on array wrappers for plottable fields
  setLabels();
}

void RelativePermeabilityBase::resizeFields( localIndex const size, localIndex const numPts )
{
  integer const numPhases = numFluidPhases();

  m_phaseRelPerm.resize( size, numPts, numPhases );
  m_phaseRelPerm_n.resize( size, numPts, numPhases );
  m_dPhaseRelPerm_dPhaseVolFrac.resize( size, numPts, numPhases, numPhases );
  //phase trapped for stats
  m_phaseTrappedVolFrac.resize( size, numPts, numPhases );
  m_phaseTrappedVolFrac.zero();
}

void RelativePermeabilityBase::setLabels()
{
  getField< fields::relperm::phaseRelPerm >().
    setDimLabels( 2, m_phaseNames );
  getField< fields::relperm::phaseRelPerm_n >().
    setDimLabels( 2, m_phaseNames );
  getField< fields::relperm::phaseTrappedVolFraction >().
    setDimLabels( 2, m_phaseNames );
}

void RelativePermeabilityBase::saveConvergedState( ) const
{
  m_phaseRelPerm_n.setValues< parallelDevicePolicy<> >( m_phaseRelPerm.toViewConst() );
}

void RelativePermeabilityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                         localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}

/// for use in RelpermDriver to browse the drainage curves
/// by setting the MaxHistoricalNonWettingSat to Snwmin and MinWettingSat to Sw
std::tuple< integer, integer > RelativePermeabilityBase::wettingAndNonWettingPhaseIndices() const
{
  using PT = PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];
  integer const ipOil = m_phaseOrder[PT::OIL];
  integer const ipGas = m_phaseOrder[PT::GAS];

  integer ipWetting = -1, ipNonWetting = -1;

  if( ipWater >= 0 && ipOil >= 0 && ipGas >= 0 )
  {
    ipWetting = ipWater;
    ipNonWetting = ipGas;
  }
  else if( ipWater < 0 )
  {
    ipWetting = ipOil;
    ipNonWetting = ipGas;
  }
  else if( ipOil < 0 )
  {
    ipWetting = ipWater;
    ipNonWetting = ipGas;
  }
  else if( ipGas < 0 )
  {
    ipWetting = ipWater;
    ipNonWetting = ipOil;
  }

  //maybe a bit too pythonic
  return std::make_tuple( ipWetting, ipNonWetting );
}

} // namespace constitutive

} // namespace geos
