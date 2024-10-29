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
 * @file ApertureTableContact.cpp
 */

#include "../deprecated/ApertureTableContact.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

ApertureTableContact::ApertureTableContact( string const & name,
                                            Group * const parent )
  : ContactBase( name, parent ),
  m_apertureTolerance( 1.0e-99 ),
  m_apertureTable( nullptr )
{
  registerWrapper( viewKeyStruct::apertureToleranceString(), &m_apertureTolerance ).
    setApplyDefaultValue( 1.0e-9 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to be used to avoid floating point errors in expressions involving aperture. "
                    "For example in the case of dividing by the actual aperture (not the effective aperture "
                    "that results from the aperture function) this value may be used to avoid the 1/0 error. "
                    "Note that this value may have some physical significance in its usage, as it may be used "
                    "to smooth out highly nonlinear behavior associated with 1/0 in addition to avoiding the "
                    "1/0 error." );

  registerWrapper( viewKeyStruct::apertureTableNameString(), &m_apertureTableName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the aperture table" );
}

ApertureTableContact::~ApertureTableContact()
{}

void ApertureTableContact::postInputInitialization()
{
  FunctionManager const & functionManager = FunctionManager::getInstance();

  GEOS_THROW_IF( m_apertureTableName.empty(),
                 getFullName() << ": the aperture table name " << m_apertureTableName << " is empty",
                 InputError );

  GEOS_THROW_IF( !functionManager.hasGroup( m_apertureTableName ),
                 getFullName() << ": the aperture table named " << m_apertureTableName << " could not be found",
                 InputError );
}

void ApertureTableContact::initializePreSubGroups()
{
  FunctionManager & functionManager = FunctionManager::getInstance();
  TableFunction & apertureTable = functionManager.getGroup< TableFunction >( m_apertureTableName );
  validateApertureTable( apertureTable );

  ArrayOfArraysView< real64 > coords = apertureTable.getCoordinates();
  arraySlice1d< real64 const > apertureValues = coords[0];
  array1d< real64 > & effectiveApertureValues = apertureTable.getValues();

  localIndex const n = apertureValues.size()-1;
  real64 const slope = ( effectiveApertureValues[n] - effectiveApertureValues[n-1] ) / ( apertureValues[n] - apertureValues[n-1] );
  real64 const apertureTransition = ( effectiveApertureValues[n] - slope * apertureValues[n] ) / ( 1.0 - slope );

  coords.emplaceBack( 0, apertureTransition );
  effectiveApertureValues.emplace_back( apertureTransition );
  coords.emplaceBack( 0, apertureTransition*10e9 );
  effectiveApertureValues.emplace_back( apertureTransition*10e9 );
  apertureTable.reInitializeFunction();

  m_apertureTable = &apertureTable;
}

void ApertureTableContact::validateApertureTable( TableFunction const & apertureTable ) const
{
  ArrayOfArraysView< real64 const > const coords = apertureTable.getCoordinates();
  arrayView1d< real64 const > const & effectiveApertureValues = apertureTable.getValues();

  GEOS_THROW_IF( coords.size() > 1,
                 getFullName() << ": Aperture limiter table cannot be greater than a 1D table.",
                 InputError );

  arraySlice1d< real64 const > apertureValues = coords[0];
  localIndex const size = apertureValues.size();

  GEOS_THROW_IF( coords( 0, size-1 ) > 0.0 || coords( 0, size-1 ) < 0.0,
                 getFullName() << ": Invalid aperture limiter table. Last coordinate must be zero!",
                 InputError );

  GEOS_THROW_IF( apertureValues.size() < 2,
                 getFullName() << ": Invalid aperture limiter table. Must have more than two points specified",
                 InputError );

  localIndex const n = apertureValues.size()-1;
  real64 const slope = ( effectiveApertureValues[n] - effectiveApertureValues[n-1] ) / ( apertureValues[n] - apertureValues[n-1] );

  GEOS_THROW_IF( slope >= 1.0,
                 getFullName() << ": Invalid aperture table. The slope of the last two points >= 1 is invalid.",
                 InputError );
}

ApertureTableContactUpdates ApertureTableContact::createKernelWrapper() const
{
  return ApertureTableContactUpdates( m_penaltyStiffness,
                                      *m_apertureTable );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ApertureTableContact, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
