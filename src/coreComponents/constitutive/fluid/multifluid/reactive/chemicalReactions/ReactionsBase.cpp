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
 * @file EquilibriumReaction.cpp
 */

#include "constitutive/fluid/multifluid/reactive/chemicalReactions/ReactionsBase.hpp"

#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{

ReactionsBase::ReactionsBase( string const & name, integer const numPrimarySpecies, integer const numSecSpecies ):
  m_name( name ),
  m_numPrimarySpecies( numPrimarySpecies ),
  m_numSecondarySpecies( numSecSpecies )
{
  // Activity coefficient related constants
  m_DebyeHuckelA = 0.5465;
  m_DebyeHuckelB = 0.3346;
  m_WATEQBDot = 0.0438;

  m_ionSizePrimary.resize( m_numPrimarySpecies );
  m_ionSizePrimary[0] = 9.00;
  m_ionSizePrimary[1] = 4.00;
  m_ionSizePrimary[2] = 6.00;
  m_ionSizePrimary[3] = 4.00;
  m_ionSizePrimary[4] = 3.00;
  m_ionSizePrimary[5] = 8.00;
  m_ionSizePrimary[6] = 4.00;

  m_ionSizeSec.resize( m_numSecondarySpecies );
  m_ionSizeSec[0] = 3.50;
  m_ionSizeSec[1] = 3.00;
  m_ionSizeSec[2] = 4.50;
  m_ionSizeSec[3] = 3.00;
  m_ionSizeSec[4] = 4.00;
  m_ionSizeSec[5] = 3.00;
  m_ionSizeSec[6] = 3.00;
  m_ionSizeSec[7] = 4.00;
  m_ionSizeSec[8] = 3.00;
  m_ionSizeSec[9] = 3.00;
  m_ionSizeSec[10] = 4.00;

  m_chargePrimary.resize( m_numPrimarySpecies );
  m_chargePrimary[0] = 1;
  m_chargePrimary[1] = -1;
  m_chargePrimary[2] = 2;
  m_chargePrimary[3] = -2;
  m_chargePrimary[4] = -1;
  m_chargePrimary[5] = 2;
  m_chargePrimary[6] = 1;

  m_chargeSec.resize( m_numSecondarySpecies );
  m_chargeSec[0] = -1;
  m_chargeSec[1] = 0;
  m_chargeSec[2] = -2;
  m_chargeSec[3] = 0;
  m_chargeSec[4] = 1;
  m_chargeSec[5] = 0;
  m_chargeSec[6] = 0;
  m_chargeSec[7] = 1;
  m_chargeSec[8] = 0;
  m_chargeSec[9] = 0;
  m_chargeSec[10] = -1;

}

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
  // localIndex const NDependent = m_numSecondarySpecies;

  GEOS_UNUSED_VAR( temperature );

  for( localIndex i = 0; i < m_numPrimarySpecies; ++i )
  {
    log10PrimaryActCoeff[i] = m_WATEQBDot * ionicStrength - m_DebyeHuckelA * m_chargePrimary[i] * m_chargePrimary[i] * sqrt( ionicStrength ) /
                              (1.0 + m_ionSizePrimary[i] * m_DebyeHuckelB * sqrt( ionicStrength ));
    dLog10PrimaryActCoeff_dIonicStrength[i] = m_WATEQBDot - m_DebyeHuckelA * m_chargePrimary[i] * m_chargePrimary[i] *
                                              (0.5 / sqrt( ionicStrength ) / (1.0 + m_ionSizePrimary[i] * m_DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * m_ionSizePrimary[i] * m_DebyeHuckelB /
                                               (1.0 + m_ionSizePrimary[i] * m_DebyeHuckelB * sqrt( ionicStrength )) /
                                               (1.0 + m_ionSizePrimary[i] * m_DebyeHuckelB * sqrt( ionicStrength )));
//    log10PrimaryActCoeff[i] = 0;
//    dLog10PrimaryActCoeff_dIonicStrength[i] = 0;
  }
  for( localIndex i = 0; i < m_numSecondarySpecies; ++i )
  {
    log10SecActCoeff[i] = m_WATEQBDot * ionicStrength - m_DebyeHuckelA * m_chargeSec[i] * m_chargeSec[i] * sqrt( ionicStrength ) /
                          (1.0 + m_ionSizeSec[i] * m_DebyeHuckelB * sqrt( ionicStrength ));
    dLog10SecActCoeff_dIonicStrength[i] = m_WATEQBDot - m_DebyeHuckelA * m_chargeSec[i] * m_chargeSec[i] *
                                          (0.5 / sqrt( ionicStrength ) / (1.0 + m_ionSizeSec[i] * m_DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * m_ionSizeSec[i] * m_DebyeHuckelB /
                                           (1.0 + m_ionSizeSec[i] * m_DebyeHuckelB * sqrt( ionicStrength )) /
                                           (1.0 + m_ionSizeSec[i] * m_DebyeHuckelB * sqrt( ionicStrength )));
//    log10SecActCoeff[i] = 0;
//    dLog10SecActCoeff_dIonicStrength[i] = 0;
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
  for( int j = 0; j < m_numSecondarySpecies; ++j )
  {
    ionicStrength += 0.5 * m_chargeSec[j] * m_chargeSec[j] * secondarySpeciesConcentration[j];
  }
}

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geos
