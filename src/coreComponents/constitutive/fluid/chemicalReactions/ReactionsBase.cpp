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
 * @file EquilibriumReaction.cpp
 */

#include "constitutive/fluid/chemicalReactions/ReactionsBase.hpp"

#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{

ReactionsBase::ReactionsBase( string const & name, integer const numPrimarySpecies, integer const numSecSpecies ):
  m_name( name ),
  m_numPrimarySpecies( numPrimarySpecies ),
  m_numSecSpecies( numSecSpecies )
{}

void ReactionsBase::KernelWrapper::computeLog10ActCoefBDotModel( real64 const temperature,
                                                                 real64 const ionicStrength,
                                                                 arraySlice1d< real64 > const & log10PrimaryActCoeff,
                                                                 arraySlice1d< real64 > const & dLog10PrimaryActCoeff_dIonicStrength,
                                                                 arraySlice1d< real64 > const & log10SecActCoeff,
                                                                 arraySlice1d< real64 > const & dLog10SecActCoeff_dIonicStrength ) const
{
  // Compute log10(ActivityCoefficient) for basis and dependent species along with their
  // derivatives with respect to Ionic strength using the B-Dot Model
  // which is the same as the Extended Debye-Huckel model in GEOS.
  // localIndex const NBasis = m_numPrimarySpecies;
  // localIndex const NDependent = m_numSecSpecies;

  real64 const TK = temperature + 273.15;

  for( localIndex i = 0; i < m_numPrimarySpecies; ++i )
  {
    log10PrimaryActCoeff[i] = m_WATEQBDot * ionicStrength - m_DebyeHuckelA * m_chargePrimary[i] * m_chargePrimary[i] * sqrt( ionicStrength ) /
                              (1.0 + m_ionSizePrimary[i] * m_DebyeHuckelB * sqrt( ionicStrength ));
    dLog10PrimaryActCoeff_dIonicStrength[i] = m_WATEQBDot - m_DebyeHuckelA * m_chargePrimary[i] * m_chargePrimary[i] *
                                              (0.5 / sqrt( ionicStrength ) / (1.0 + m_ionSizePrimary[i] * m_DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * m_ionSizePrimary[i] * m_DebyeHuckelB /
                                               (1.0 + m_ionSizePrimary[i] * m_DebyeHuckelB * sqrt( ionicStrength )) /
                                               (1.0 + m_ionSizePrimary[i] * m_DebyeHuckelB * sqrt( ionicStrength )));
  }
  for( localIndex i = 0; i < m_numSecSpecies; ++i )
  {
    log10SecActCoeff[i] = m_WATEQBDot * ionicStrength - m_DebyeHuckelA * m_chargeSec[i] * m_chargeSec[i] * sqrt( ionicStrength ) /
                          (1.0 + m_ionSizeSec[i] * m_DebyeHuckelB * sqrt( ionicStrength ));
    dLog10SecActCoeff_dIonicStrength[i] = m_WATEQBDot - m_DebyeHuckelA * m_chargeSec[i] * m_chargeSec[i] *
                                          (0.5 / sqrt( ionicStrength ) / (1.0 + m_ionSizeSec[i] * m_DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * m_ionSizeSec[i] * m_DebyeHuckelB /
                                           (1.0 + m_ionSizeSec[i] * m_DebyeHuckelB * sqrt( ionicStrength )) /
                                           (1.0 + m_ionSizeSec[i] * m_DebyeHuckelB * sqrt( ionicStrength )));
  }
}

void ReactionsBase::KernelWrapper::computeIonicStrength( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                                         real64 & ionicStrength ) const
{
  //get ionic strength
  ionicStrength = 0.0;
  // Primary species
  for( localIndex i = 0; i < m_numPrimarySpecies; ++i )
  {
    ionicStrength += 0.5 * m_chargePrimary[i] * m_chargePrimary[i] * primarySpeciesConcentration[i];
  }
  // Secondary species
  for( int j = 0; j < m_numSecSpecies; ++j )
  {
    ionicStrength += 0.5 * m_chargeSec[j] * m_chargeSec[j] * secondarySpeciesConcentration[j];
  }
}

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geosx
