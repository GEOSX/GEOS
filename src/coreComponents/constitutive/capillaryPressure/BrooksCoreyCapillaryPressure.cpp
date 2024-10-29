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
 * @file BrooksCoreyCapillaryPressure.cpp
 */

#include "BrooksCoreyCapillaryPressure.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

BrooksCoreyCapillaryPressure::BrooksCoreyCapillaryPressure( string const & name,
                                                            Group * const parent )
  : CapillaryPressureBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::phaseCapPressureExponentInvString(), &m_phaseCapPressureExponentInv ).
    setApplyDefaultValue( 2.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Inverse of capillary power law exponent for each phase" );

  registerWrapper( viewKeyStruct::phaseEntryPressureString(), &m_phaseEntryPressure ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Entry pressure value for each phase" );

  registerWrapper( viewKeyStruct::capPressureEpsilonString(), &m_capPressureEpsilon ).
    setApplyDefaultValue( 1e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription(
    "Wetting-phase saturation at which the max cap. pressure is attained; used to avoid infinite cap. pressure values for saturations close to zero" );

  registerWrapper( viewKeyStruct::volFracScaleString(), &m_volFracScale ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Factor used to scale the phase capillary pressure, defined as: one minus the sum of the phase minimum volume fractions." );
}

void BrooksCoreyCapillaryPressure::postInputInitialization()
{
  CapillaryPressureBase::postInputInitialization();

  auto const checkInputSize = [&]( auto const & array, auto const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), m_phaseNames.size(),
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );
  };
  checkInputSize( m_phaseMinVolumeFraction, viewKeyStruct::phaseMinVolumeFractionString() );
  checkInputSize( m_phaseCapPressureExponentInv, viewKeyStruct::phaseCapPressureExponentInvString() );
  checkInputSize( m_phaseEntryPressure, viewKeyStruct::phaseEntryPressureString() );

  m_volFracScale = 1.0;
  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    auto const errorMsg = [&]( auto const & attribute )
    {
      return GEOS_FMT( "{}: invalid value at {}[{}]", getFullName(), attribute, ip );
    };

    GEOS_THROW_IF_LT_MSG( m_phaseMinVolumeFraction[ip], 0.0,
                          errorMsg( viewKeyStruct::phaseMinVolumeFractionString() ),
                          InputError );
    GEOS_THROW_IF_GT_MSG( m_phaseMinVolumeFraction[ip], 1.0,
                          errorMsg( viewKeyStruct::phaseMinVolumeFractionString() ),
                          InputError );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];

    if( m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE )
    {
      GEOS_THROW_IF_LE_MSG( m_phaseCapPressureExponentInv[ip], 0.0,
                            errorMsg( viewKeyStruct::phaseCapPressureExponentInvString() ),
                            InputError );
      GEOS_THROW_IF_LT_MSG( m_phaseEntryPressure[ip], 0.0,
                            errorMsg( viewKeyStruct::phaseEntryPressureString() ),
                            InputError );
      GEOS_THROW_IF_LT_MSG( m_capPressureEpsilon, 0.0,
                            errorMsg( viewKeyStruct::capPressureEpsilonString() ),
                            InputError );
      GEOS_THROW_IF_GT_MSG( m_capPressureEpsilon, 0.2,
                            errorMsg( viewKeyStruct::capPressureEpsilonString() ),
                            InputError );
    }
  }

  GEOS_THROW_IF_LT_MSG( m_volFracScale, 0.0,
                        GEOS_FMT( "{}: sum of min volume fractions exceeds 1.0", getFullName() ),
                        InputError );
}

BrooksCoreyCapillaryPressure::KernelWrapper
BrooksCoreyCapillaryPressure::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_phaseCapPressureExponentInv,
                        m_phaseEntryPressure,
                        m_capPressureEpsilon,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyCapillaryPressure, string const &, Group * const )
} // namespace constitutive

} // namespace geos
