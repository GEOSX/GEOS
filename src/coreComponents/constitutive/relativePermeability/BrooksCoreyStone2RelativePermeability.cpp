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
 * @file BrooksCoreyStone2RelativePermeability.cpp
 */

#include "BrooksCoreyStone2RelativePermeability.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

BrooksCoreyStone2RelativePermeability::BrooksCoreyStone2RelativePermeability( string const & name,
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



void BrooksCoreyStone2RelativePermeability::resizeFields( localIndex const size, localIndex const numPts )
{
  RelativePermeabilityBase::resizeFields( size, numPts );

  integer const numPhases = numFluidPhases();
  integer const numDir = m_waterOilRelPermExponent.size(0);


  m_phaseRelPerm.resize( size, numPts, numPhases, numDir );
  m_phaseRelPerm_n.resize( size, numPts, numPhases, numDir );
  m_dPhaseRelPerm_dPhaseVolFrac.resize( size, numPts, numPhases, numPhases, numDir );
  //phase trapped for stats
  m_phaseTrappedVolFrac.resize( size, numPts, numPhases );
  m_phaseTrappedVolFrac.zero();


}




void BrooksCoreyStone2RelativePermeability::postInputInitialization()
{
  RelativePermeabilityBase::postInputInitialization();

  integer const numDir = 3;

  m_volFracScale.resize( numDir /*ndims*/ );

  GEOS_THROW_IF( m_phaseOrder[PhaseType::OIL] < 0,
                 GEOS_FMT( "{}: reference oil phase has not been defined and must be included in model", getFullName() ),
                 InputError );

  auto const checkInputSize = [&]( auto const & array, localIndex const expected, auto const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), expected,
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );
  };

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


  for( integer ip = 0; ip < numFluidPhases(); ++ip )
  {
    auto const errorMsg = [&]( auto const & attribute ) {
      return GEOS_FMT( "{}: invalid value at {}[{}]", getFullName(), attribute, ip );
    };
    for( int dir=0; dir<numDir; ++dir )
    {
      checkInputSize( m_phaseMinVolumeFraction[dir], numFluidPhases(), viewKeyStruct::phaseMinVolumeFractionString() );
      m_volFracScale[dir] = 1.0;
      GEOS_THROW_IF_LT_MSG( m_phaseMinVolumeFraction[dir][ip], 0.0,
                            errorMsg( viewKeyStruct::phaseMinVolumeFractionString()),
                            InputError );
      GEOS_THROW_IF_GT_MSG( m_phaseMinVolumeFraction[dir][ip], 1.0,
                            errorMsg( viewKeyStruct::phaseMinVolumeFractionString()),
                            InputError );
      m_volFracScale[dir] -= m_phaseMinVolumeFraction[dir][ip];
    }


  }


  for( int dir=0; dir<numDir; ++dir )
  {
    GEOS_THROW_IF_LT_MSG( m_volFracScale[dir], 0.0,
                          GEOS_FMT( "{}: sum of min volume fractions exceeds 1.0", getFullName()),
                          InputError );
    for( integer ip = 0; ip < 2; ++ip )
    {
      auto const errorMsg = [&]( auto const & attribute ) {
        return GEOS_FMT( "{}: invalid value at {}[{}]", getFullName(), attribute, ip );
      };
      if( m_phaseOrder[PhaseType::WATER] >= 0 )
      {
        GEOS_THROW_IF_LT_MSG( m_waterOilRelPermExponent[dir][ip], 0.0,
                              errorMsg( viewKeyStruct::waterOilRelPermExponentString()),
                              InputError );
        GEOS_THROW_IF_LT_MSG( m_waterOilRelPermMaxValue[dir][ip], 0.0,
                              errorMsg( viewKeyStruct::waterOilRelPermMaxValueString()),
                              InputError );
        GEOS_THROW_IF_GT_MSG( m_waterOilRelPermMaxValue[dir][ip], 1.0,
                              errorMsg( viewKeyStruct::waterOilRelPermMaxValueString()),
                              InputError );
      }

      if( m_phaseOrder[PhaseType::GAS] >= 0 )
      {
        GEOS_THROW_IF_LT_MSG( m_gasOilRelPermExponent[dir][ip], 0.0,
                              errorMsg( viewKeyStruct::gasOilRelPermExponentString()),
                              InputError );
        GEOS_THROW_IF_LT_MSG( m_gasOilRelPermMaxValue[dir][ip], 0.0,
                              errorMsg( viewKeyStruct::gasOilRelPermMaxValueString()),
                              InputError );
        GEOS_THROW_IF_GT_MSG( m_gasOilRelPermMaxValue[dir][ip], 1.0,
                              errorMsg( viewKeyStruct::gasOilRelPermMaxValueString()),
                              InputError );
      }
    }

    if( m_phaseOrder[PhaseType::WATER] >= 0 && m_phaseOrder[PhaseType::GAS] >= 0 )
    {
      real64 const mean = 0.5 * (m_gasOilRelPermMaxValue[dir][GasOilPairPhaseType::OIL]
                                 + m_waterOilRelPermMaxValue[dir][WaterOilPairPhaseType::OIL]);
      m_gasOilRelPermMaxValue[dir][GasOilPairPhaseType::OIL] = mean;
      m_waterOilRelPermMaxValue[dir][WaterOilPairPhaseType::OIL] = mean;
    }
  }
}

BrooksCoreyStone2RelativePermeability::KernelWrapper
BrooksCoreyStone2RelativePermeability::createKernelWrapper()
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

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyStone2RelativePermeability, string const &, Group * const )
}     // namespace constitutive

} // namespace geos
