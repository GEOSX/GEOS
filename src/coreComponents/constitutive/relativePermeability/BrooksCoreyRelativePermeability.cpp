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
 * @file BrooksCoreyRelativePermeability.cpp
 */

#include "BrooksCoreyRelativePermeability.hpp"

#include <cmath>

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

//START_SPHINX_INCLUDE_00

BrooksCoreyRelativePermeability::BrooksCoreyRelativePermeability( string const & name,
                                                                  Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::phaseRelPermExponentString(), &m_phaseRelPermExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum relative permeability power law exponent for each phase" );


  registerWrapper( viewKeyStruct::phaseRelPermMaxValueString(), &m_phaseRelPermMaxValue ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum relative permeability value for each phase" );

  registerWrapper( viewKeyStruct::volFracScaleString(), &m_volFracScale ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Factor used to scale the phase relative permeability, defined as: one minus the sum of the phase minimum volume fractions." );

}

void BrooksCoreyRelativePermeability::postInputInitialization()
{
  RelativePermeabilityBase::postInputInitialization();

  auto const checkInputSize = [&]( auto const & array, auto const & attribute )
  {
    GEOS_THROW_IF_NE_MSG( array.size(), m_phaseNames.size(),
                          GEOS_FMT( "{}: invalid number of values in attribute '{}'", getFullName(), attribute ),
                          InputError );
  };
  checkInputSize( m_phaseMinVolumeFraction, viewKeyStruct::phaseMinVolumeFractionString() );
  checkInputSize( m_phaseRelPermExponent, viewKeyStruct::phaseRelPermExponentString() );
  checkInputSize( m_phaseRelPermMaxValue, viewKeyStruct::phaseRelPermMaxValueString() );

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

    GEOS_THROW_IF_LT_MSG( m_phaseRelPermExponent[ip], 0.0,
                          errorMsg( viewKeyStruct::phaseRelPermExponentString() ),
                          InputError );
    GEOS_THROW_IF_LT_MSG( m_phaseRelPermMaxValue[ip], 0.0,
                          errorMsg( viewKeyStruct::phaseRelPermMaxValueString() ),
                          InputError );
    GEOS_THROW_IF_GT_MSG( m_phaseRelPermMaxValue[ip], 1.0,
                          errorMsg( viewKeyStruct::phaseRelPermMaxValueString() ),
                          InputError );
  }

  GEOS_THROW_IF_LT_MSG( m_volFracScale, 0.0,
                        GEOS_FMT( "{}: sum of min volume fractions exceeds 1.0", getFullName() ),
                        InputError );
}

BrooksCoreyRelativePermeability::KernelWrapper
BrooksCoreyRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_phaseRelPermExponent,
                        m_phaseRelPermMaxValue,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac,
                        m_phaseTrappedVolFrac );
}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyRelativePermeability, string const &, Group * const )

} // namespace constitutive

} // namespace geos
