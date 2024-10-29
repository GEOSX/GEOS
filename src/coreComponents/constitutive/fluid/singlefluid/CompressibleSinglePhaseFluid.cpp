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
 * @file CompressibleSinglePhaseFluid.cpp
 */

#include "CompressibleSinglePhaseFluid.hpp"

#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CompressibleSinglePhaseFluid::CompressibleSinglePhaseFluid( string const & name, Group * const parent ):
  SingleFluidBase( name, parent ),
  m_densityModelType( ExponentApproximationType::Linear ),
  m_viscosityModelType( ExponentApproximationType::Linear )
{
  registerWrapper( viewKeyStruct::defaultDensityString(), &m_defaultDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value for density." );

  registerWrapper( viewKeyStruct::defaultViscosityString(), &m_defaultViscosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value for viscosity." );

  registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid compressibility" );

  registerWrapper( viewKeyStruct::viscosibilityString(), &m_viscosibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fluid viscosity exponential coefficient" );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::referenceDensityString(), &m_referenceDensity ).
    setApplyDefaultValue( 1000.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid density" );

  registerWrapper( viewKeyStruct::referenceViscosityString(), &m_referenceViscosity ).
    setApplyDefaultValue( 0.001 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference fluid viscosity" );

  registerWrapper( viewKeyStruct::densityModelTypeString(), &m_densityModelType ).
    setApplyDefaultValue( m_densityModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of density model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::viscosityModelTypeString(), &m_viscosityModelType ).
    setApplyDefaultValue( m_viscosityModelType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Type of viscosity model. Valid options:\n* " + EnumStrings< ExponentApproximationType >::concat( "\n* " ) );

}

CompressibleSinglePhaseFluid::~CompressibleSinglePhaseFluid() = default;

void CompressibleSinglePhaseFluid::allocateConstitutiveData( dataRepository::Group & parent,
                                                             localIndex const numConstitutivePointsPerParentIndex )
{
  SingleFluidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  getField< fields::singlefluid::density >().setApplyDefaultValue( m_defaultDensity );
  getField< fields::singlefluid::viscosity >().setApplyDefaultValue( m_defaultViscosity );

  m_density.setValues< serialPolicy >( m_referenceDensity );
  m_viscosity.setValues< serialPolicy >( m_referenceViscosity );
}

void CompressibleSinglePhaseFluid::postInputInitialization()
{
  SingleFluidBase::postInputInitialization();

  auto const checkNonnegative = [&]( real64 const value, auto const & attribute )
  {
    GEOS_THROW_IF_LT_MSG( value, 0.0,
                          GEOS_FMT( "{}: invalid value of attribute '{}'", getFullName(), attribute ),
                          InputError );
  };
  checkNonnegative( m_compressibility, viewKeyStruct::compressibilityString() );
  checkNonnegative( m_viscosibility, viewKeyStruct::viscosibilityString() );

  auto const checkPositive = [&]( real64 const value, auto const & attribute )
  {
    GEOS_THROW_IF_LE_MSG( value, 0.0,
                          GEOS_FMT( "{}: invalid value of attribute '{}'", getFullName(), attribute ),
                          InputError );
  };
  checkPositive( m_referenceDensity, viewKeyStruct::referenceDensityString() );
  checkPositive( m_referenceViscosity, viewKeyStruct::referenceViscosityString() );

  // Due to the way update wrapper is currently implemented, we can only support one model type
  auto const checkModelType = [&]( ExponentApproximationType const value, auto const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( value, ExponentApproximationType::Linear,
                          GEOS_FMT( "{}: invalid model type in attribute '{}' (only linear currently supported)", getFullName(), attribute ),
                          InputError );
  };
  checkModelType( m_densityModelType, viewKeyStruct::densityModelTypeString() );
  checkModelType( m_viscosityModelType, viewKeyStruct::viscosityModelTypeString() );

  // Set default values for derivatives (cannot be done in base class)
  // TODO: reconsider the necessity of this

  real64 dRho_dP;
  real64 dVisc_dP;
  createKernelWrapper().compute( m_referencePressure, m_referenceDensity, dRho_dP, m_referenceViscosity, dVisc_dP );
  getField< fields::singlefluid::dDensity_dPressure >().setDefaultValue( dRho_dP );
  getField< fields::singlefluid::dViscosity_dPressure >().setDefaultValue( dVisc_dP );
}

CompressibleSinglePhaseFluid::KernelWrapper
CompressibleSinglePhaseFluid::createKernelWrapper()
{
  return KernelWrapper( KernelWrapper::DensRelationType( m_referencePressure, m_referenceDensity, m_compressibility ),
                        KernelWrapper::ViscRelationType( m_referencePressure, m_referenceViscosity, m_viscosibility ),
                        m_density,
                        m_dDensity_dPressure,
                        m_viscosity,
                        m_dViscosity_dPressure );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleSinglePhaseFluid, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
