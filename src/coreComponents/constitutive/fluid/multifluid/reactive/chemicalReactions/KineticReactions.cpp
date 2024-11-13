/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file KineticReactions.cpp
 */

#include "constitutive/fluid/multifluid/reactive/chemicalReactions/KineticReactions.hpp"
#include "functions/FunctionManager.hpp"
#include "common/Units.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{

KineticReactions::KineticReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies, integer const numKineticReactions ):
  ReactionsBase( name, numPrimarySpecies, numSecSpecies ),
  m_numKineticReactions( numKineticReactions )
{

  // Stochiometric Matrix for the kinetic reactions (in terms of primary species only)
  // First index: 0 = Ca(OH)2 dissolution, 1 = CaCO3 dissolution
  // Second index: 0 = H+, 1 = HCO3-, 2 = Ca+2, 3 = SO4-2, 4 = Cl-, 5 = Mg+2, 6 = Na+1
  m_stoichMatrix.resize( m_numKineticReactions, m_numPrimarySpecies );
  m_stoichMatrix[0][0] = -2;
  m_stoichMatrix[0][2] = 1;
  m_stoichMatrix[1][0] = -1;
  m_stoichMatrix[1][1] = 1;
  m_stoichMatrix[1][2] = 1;

  // Equilibrium Constant
  m_log10EqConst.resize( m_numKineticReactions );
  m_log10EqConst[0] = 20.19;
  m_log10EqConst[1] = 1.32;

  // Rate Constant
  // have to check the values as the functional form of the rate equation here differs from the implementation in GEOS
  m_reactionRateConstant.resize( m_numKineticReactions );
  m_reactionRateConstant[0] = 9.95e-1;
  m_reactionRateConstant[1] = 9.95e-3;
  // Here we should either read the database or the input values.

  m_specificSurfaceArea = 1.0;
}

KineticReactions::KernelWrapper KineticReactions::createKernelWrapper() const
{
  return KernelWrapper( m_numPrimarySpecies,
                        m_numSecondarySpecies,
                        m_numKineticReactions,
                        m_log10EqConst,
                        m_stoichMatrix,
                        m_chargePrimary,
                        m_chargeSec,
                        m_ionSizePrimary,
                        m_ionSizeSec,
                        m_DebyeHuckelA,
                        m_DebyeHuckelB,
                        m_WATEQBDot,
                        m_reactionRateConstant,
                        m_specificSurfaceArea );
}

// function to  the reaction rate. Includes impact of temperature, concentration, surface area, volume fraction and porosity
void KineticReactions::KernelWrapper::computeReactionRates( real64 const & temperature,
                                                            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                            arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                                            arraySlice1d< real64, compflow::USD_COMP - 1 > const & reactionRates ) const
{
  /// 1. Create local vectors
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > log10PrimaryActCoeff( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumSecondarySpecies > log10SecActCoeff( m_numSecondarySpecies );

  real64 ionicStrength = 0.0;
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > dIonicStrength_dPrimaryConcentration( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > dLog10PrimaryActCoeff_dIonicStrength( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumSecondarySpecies > dLog10SecActCoeff_dIonicStrength( m_numSecondarySpecies );

  /// 2. Compute activity coefficients
  computeIonicStrength( primarySpeciesConcentration,
                        secondarySpeciesConcentration,
                        ionicStrength );

  computeLog10ActCoefBDotModel( temperature,
                                ionicStrength,
                                log10PrimaryActCoeff,
                                dLog10PrimaryActCoeff_dIonicStrength,
                                log10SecActCoeff,
                                dLog10SecActCoeff_dIonicStrength );


  /// 3. Compute reaction rates
  for( int iRxn = 0; iRxn < m_numKineticReactions; iRxn++ )
  {
    real64 saturationIndex = -m_log10EqConst[iRxn];

    for( int iPri = 0; iPri < m_numPrimarySpecies; ++iPri )
    {
      saturationIndex += m_stoichMatrix[iRxn][iPri] * log10( primarySpeciesConcentration[iPri] );
      saturationIndex += m_stoichMatrix[iRxn][iPri] * log10PrimaryActCoeff[iPri];
    }

    reactionRates[iRxn] = m_specificSurfaceArea * (1.0 - pow( 10, saturationIndex ) ) * m_reactionRateConstant[iRxn];
  }
}

} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geos

/*
   for( localIndex ir = 0; ir < NReaction; ++ir )
   {

    for( localIndex ic = 0; ic < kineticReaction.stochs.size(); ++ic )
    {
      dSIndex[id] = kineticReaction.stochs[ic];

    }


    for( localIndex ic = 0; ic < NBasis; ++ic )
    {

      dKineticReactionRate_dConc[ir][ic] = S * kineticReaction.rateConst * rateTemp * pow( 10.0, SIndex ) * log( 10.0 ) * dSIndex[ic];

    }

   }


   for( localIndex ir = 0; ir < NReaction; ++ir )
   {

    for( localIndex i = 0; i < (kineticReaction.stochs).size(); ++i )
    {

      localIndex ic = basisSpeciesIndices[i];

      kineticSpeciesReactionRate[ic] += -(kineticReaction.stochs)[i] * kineticReactionRate[ir];

      for( localIndex id = 0; id < NBasis; ++id )
      {

        dKineticSpeciesReactionRate_dConc[ic][id] += -(kineticReaction.stochs)[i] * dKineticReactionRate_dConc[ir][id];

      }

    }

   }

   }
   // arraySlice1d< real64 const > const & concentration,
   // arraySlice1d< real64 const > const & surfaceArea0,
   // arraySlice1d< real64 const > const & volumeFraction0,
   // arraySlice1d< real64 const > const & volumeFraction,
   // real64 const & porosity0,
   // real64 const & porosity,
   // arraySlice1d< real64 > const & kineticReactionRate )
   // {
   // for( localIndex ir = 0; ir < kineticReactionArray.size(); ++ir )
   // {
   //   const KineticReaction & kineticReaction = kineticReactionArray[ir];
   //   const array1d< localIndex > & basisSpeciesIndices = kineticReaction.basisSpeciesIndices;
   //   // calculate saturation index
   //   real64 SIndex = -kineticReaction.logK;
   //   for( localIndex ic = 0; ic < kineticReaction.stochs.size(); ++ic )
   //   {
   //     SIndex += kineticReaction.stochs[ic] * concentration[basisSpeciesIndices[ic] ];     // Check that the input "concentration" is
   //                                                                                         // actually ln(activity
   // coefficient*concentration)
   //   }
   //   // surface area is assumed to scale with the volume fraction. Check whether this is the volume fraction of the mineral
   //   // dissolving/precipitating. Not sure why porosity is included.
   //   real64 S = surfaceArea0[ir] * pow( volumeFraction[ir] / volumeFraction0[ir], 2.0/3.0 ) * pow( porosity / porosity0, 2.0/3.0 );
   //   // computing the rate at the correct temperature. Looks like EQ36 database has it at 298.15 K
   //   real64 rateTemp = exp( -kineticReaction.E / RConst * (1.0 / units::convertCToK( temperature ) - 1.0 / 298.15));
   //   real64 SS = (pow( 10.0, SIndex ) - 1.0);
   //   kineticReactionRate[ir] = S * kineticReaction.rateConst * rateTemp * SS;
   //}


 */
