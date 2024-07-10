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
 * @file BrooksCoreyBakerRelativePermeability.cpp
 */

#include "BrooksCoreyBakerRelativePermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

BrooksCoreyBakerRelativePermeability::BrooksCoreyBakerRelativePermeability( string const & name,
                                                                            Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::waterOilRelPermExponentString(), &m_waterOilRelPermExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rel perm power law exponent for the pair (water phase, oil phase) at residual gas saturation\n"
                    "The expected format is \"{ waterExp, oilExp }\", in that order" );

  registerWrapper( viewKeyStruct::waterOilRelPermMaxValueString(), &m_waterOilRelPermMaxValue ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum rel perm value for the pair (water phase, oil phase) at residual gas saturation\n"
                    "The expected format is \"{ waterMax, oilMax }\", in that order" );

  registerWrapper( viewKeyStruct::gasOilRelPermExponentString(), &m_gasOilRelPermExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rel perm power law exponent for the pair (gas phase, oil phase) at residual water saturation\n"
                    "The expected format is \"{ gasExp, oilExp }\", in that order" );

  registerWrapper( viewKeyStruct::gasOilRelPermMaxValueString(), &m_gasOilRelPermMaxValue ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum rel perm value for the pair (gas phase, oil phase) at residual water saturation\n"
                    "The expected format is \"{ gasMax, oilMax }\", in that order" );

  registerWrapper( viewKeyStruct::volFracScaleString(), &m_volFracScale ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Factor used to scale the phase capillary pressure, defined as: one minus the sum of the phase minimum volume fractions." );

}

void BrooksCoreyBakerRelativePermeability::postInputInitialization()
{
  RelativePermeabilityBase::postInputInitialization();

  GEOS_THROW_IF( m_phaseOrder[PhaseType::OIL] < 0,
                 GEOS_FMT( "{}: reference oil phase has not been defined and must be included in model", getFullName() ),
                 InputError );

  auto const checkInputSize = [&]( auto const & array, localIndex const expected, auto const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), expected,
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );
  };
  checkInputSize( m_phaseMinVolumeFraction, numFluidPhases(), viewKeyStruct::phaseMinVolumeFractionString() );

  if( m_phaseOrder[PhaseType::WATER] >= 0 )
  {
    checkInputSize( m_waterOilRelPermExponent, 2, viewKeyStruct::waterOilRelPermExponentString() );
    checkInputSize( m_waterOilRelPermMaxValue, 2, viewKeyStruct::waterOilRelPermMaxValueString() );
  }

  if( m_phaseOrder[PhaseType::GAS] >=0 )
  {
    checkInputSize( m_gasOilRelPermExponent, 2, viewKeyStruct::gasOilRelPermExponentString() );
    checkInputSize( m_gasOilRelPermMaxValue, 2, viewKeyStruct::gasOilRelPermMaxValueString() );
  }

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
  }

  GEOS_THROW_IF_LT_MSG( m_volFracScale, 0.0,
                        GEOS_FMT( "{}: sum of min volume fractions exceeds 1.0", getFullName() ),
                        InputError );


  for( integer ip = 0; ip < 2; ++ip )
  {
    auto const errorMsg = [&]( auto const & attribute )
    {
      return GEOS_FMT( "{}: invalid value at {}[{}]", getFullName(), attribute, ip );
    };
    if( m_phaseOrder[PhaseType::WATER] >= 0 )
    {
      GEOS_THROW_IF_LT_MSG( m_waterOilRelPermExponent[ip], 0.0,
                            errorMsg( viewKeyStruct::waterOilRelPermExponentString() ),
                            InputError );
      GEOS_THROW_IF_LT_MSG( m_waterOilRelPermMaxValue[ip], 0.0,
                            errorMsg( viewKeyStruct::waterOilRelPermMaxValueString() ),
                            InputError );
      GEOS_THROW_IF_GT_MSG( m_waterOilRelPermMaxValue[ip], 1.0,
                            errorMsg( viewKeyStruct::waterOilRelPermMaxValueString() ),
                            InputError );
    }

    if( m_phaseOrder[PhaseType::GAS] >= 0 )
    {
      GEOS_THROW_IF_LT_MSG( m_gasOilRelPermExponent[ip], 0.0,
                            errorMsg( viewKeyStruct::gasOilRelPermExponentString() ),
                            InputError );
      GEOS_THROW_IF_LT_MSG( m_gasOilRelPermMaxValue[ip], 0.0,
                            errorMsg( viewKeyStruct::gasOilRelPermMaxValueString() ),
                            InputError );
      GEOS_THROW_IF_GT_MSG( m_gasOilRelPermMaxValue[ip], 1.0,
                            errorMsg( viewKeyStruct::gasOilRelPermMaxValueString() ),
                            InputError );
    }
  }

  if( m_phaseOrder[PhaseType::WATER] >= 0 && m_phaseOrder[PhaseType::GAS] >= 0 )
  {
    real64 const mean = 0.5 * ( m_gasOilRelPermMaxValue[GasOilPairPhaseType::OIL]
                                + m_waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL] );
    m_gasOilRelPermMaxValue[GasOilPairPhaseType::OIL]     = mean;
    m_waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL] = mean;
  }
}

BrooksCoreyBakerRelativePermeability::KernelWrapper
BrooksCoreyBakerRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_waterOilRelPermExponent,
                        m_waterOilRelPermMaxValue,
                        m_gasOilRelPermExponent,
                        m_gasOilRelPermMaxValue,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac,
                        m_phaseTrappedVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyBakerRelativePermeability, string const &, Group * const )
} // namespace constitutive

} // namespace geos
