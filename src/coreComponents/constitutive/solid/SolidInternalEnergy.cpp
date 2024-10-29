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
 * @file SolidInternalEnergy.cpp
 */

#include "SolidInternalEnergy.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SolidInternalEnergy::SolidInternalEnergy( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_internalEnergy(),
  m_dInternalEnergy_dTemperature(),
  m_referenceVolumetricHeatCapacity(),
  m_dVolumetricHeatCapacity_dTemperature(),
  m_referenceTemperature(),
  m_referenceInternalEnergy()
{
  registerWrapper( viewKeyStruct::internalEnergyString(), &m_internalEnergy ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Internal energy of the solid per unit volume [J/m^3]" );

  registerWrapper( viewKeyStruct::oldInternalEnergyString(), &m_internalEnergy_n ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Internal energy of the solid per unit volume at the previous time-step [J/m^3]" );

  registerWrapper( viewKeyStruct::dInternalEnergy_dTemperatureString(), &m_dInternalEnergy_dTemperature ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Derivative of the solid internal energy w.r.t. temperature [J/(m^3.K)]" );

  registerWrapper( viewKeyStruct::referenceVolumetricHeatCapacityString(), &m_referenceVolumetricHeatCapacity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference solid volumetric heat capacity [J/(kg.K)]" );

  registerWrapper( viewKeyStruct::dVolumetricHeatCapacity_dTemperatureString(), &m_dVolumetricHeatCapacity_dTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Derivative of the solid volumetric heat capacity w.r.t. temperature [J/(m^3.K^2)]" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference temperature [K]" );

  registerWrapper( viewKeyStruct::referenceInternalEnergyString(), &m_referenceInternalEnergy ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Internal energy at the reference temperature [J/kg]" );
}

void SolidInternalEnergy::allocateConstitutiveData( Group & parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  m_internalEnergy.resize( 0, 1 );
  m_dInternalEnergy_dTemperature.resize( 0, 1 );
  m_internalEnergy_n.resize( 0, 1 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void SolidInternalEnergy::saveConvergedState() const
{
  arrayView2d< real64 const > internalEnergy   = m_internalEnergy;
  arrayView2d< real64 >       internalEnergy_n = m_internalEnergy_n;

  forAll< parallelDevicePolicy<> >( internalEnergy.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    internalEnergy_n[k][0] = internalEnergy[k][0];
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SolidInternalEnergy, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
