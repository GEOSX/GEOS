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
 * @file EquilibriumReactions.cpp
 */

#include "EquilibriumReactions.hpp"

#include "functions/FunctionManager.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace chemicalReactions
{

EquilibriumReactions::EquilibriumReactions( string const & name, integer const numPrimarySpecies, integer const numSecSpecies ):
  ReactionsBase( name, numPrimarySpecies, numSecSpecies )
{
  // Here we should either read the database or the input values.

  // Hardcoding values for now



  // Stochiometric Matrix
  // First index: 0 = OH-, 1 = CO2, 2 = CO3-2, 3 = H2CO3, 4 = CaHCO3+, 5 = CaCO3, 6 = CaSO4, 7 = CaCl+, 8 = CaCl2, 9 = MgSO4, 10 = NaSO4-
  // Second index: 0 = H+, 1 = HCO3-, 2 = Ca+2, 3 = SO4-2, 4 = Cl-, 5 = Mg+2, 6 = Na+1
  m_stoichMatrix.resize( m_numSecondarySpecies, m_numPrimarySpecies );
  m_stoichMatrix[0][0] = -1;
  m_stoichMatrix[1][0] = 1;
  m_stoichMatrix[1][1] = 1;
  m_stoichMatrix[2][0] = -1;
  m_stoichMatrix[2][1] = 1;
  m_stoichMatrix[3][0] = 1;
  m_stoichMatrix[3][1] = 1;
  m_stoichMatrix[4][1] = 1;
  m_stoichMatrix[4][2] = 1;
  m_stoichMatrix[5][0] = -1;
  m_stoichMatrix[5][1] = 1;
  m_stoichMatrix[5][2] = 1;
  m_stoichMatrix[6][2] = 1;
  m_stoichMatrix[6][3] = 1;
  m_stoichMatrix[7][2] = 1;
  m_stoichMatrix[7][4] = 1;
  m_stoichMatrix[8][2] = 1;
  m_stoichMatrix[8][4] = 2;
  m_stoichMatrix[9][5] = 1;
  m_stoichMatrix[9][3] = 1;
  m_stoichMatrix[10][6] = 1;
  m_stoichMatrix[10][3] = 1;

  // Equilibrium Constant
  m_log10EqConst.resize( m_numSecondarySpecies );
  m_log10EqConst[0] = 13.99;
  m_log10EqConst[1] = -6.36;
  m_log10EqConst[2] = 10.33;
  m_log10EqConst[3] = -3.77;
  m_log10EqConst[4] = -1.09;
  m_log10EqConst[5] = 7.07;
  m_log10EqConst[6] = -2.16;
  m_log10EqConst[7] = 0.67;
  m_log10EqConst[8] = 0.60;
  m_log10EqConst[9] = -2.43;
  m_log10EqConst[10] = -0.82;
}

EquilibriumReactions::KernelWrapper EquilibriumReactions::createKernelWrapper() const
{
  return KernelWrapper( m_numPrimarySpecies,
                        m_numSecondarySpecies,
                        m_log10EqConst,
                        m_stoichMatrix,
                        m_chargePrimary,
                        m_chargeSec,
                        m_ionSizePrimary,
                        m_ionSizeSec,
                        m_DebyeHuckelA,
                        m_DebyeHuckelB,
                        m_WATEQBDot );
}

void EquilibriumReactions::KernelWrapper::assembleEquilibriumReactionSystem( real64 const temperature,
                                                                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                                                             arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                             arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration,
                                                                             arraySlice2d< real64 > const & matrix,
                                                                             arraySlice1d< real64 > const & rhs ) const
{

  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > log10PrimaryActCoeff( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumSecondarySpecies > log10SecActCoeff( m_numSecondarySpecies );
  stackArray2d< real64, ReactionsBase::maxNumSecondarySpecies * ReactionsBase::maxNumPrimarySpecies > dLog10SecConc_dLog10PrimaryConc( m_numSecondarySpecies, m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > totalConcentration( m_numPrimarySpecies );
  stackArray2d< real64, ReactionsBase::maxNumPrimarySpecies * ReactionsBase::maxNumPrimarySpecies > dTotalConc_dLog10PrimaryConc( m_numPrimarySpecies, m_numPrimarySpecies );

  real64 ionicStrength = 0.0;
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > dIonicStrength_dPrimaryConcentration( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > dLog10PrimaryActCoeff_dIonicStrength( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumSecondarySpecies > dLog10SecActCoeff_dIonicStrength( m_numSecondarySpecies );

  /// activity coefficients
  computeIonicStrength( primarySpeciesConcentration,
                        secondarySpeciesConcentration,
                        ionicStrength );

  computeLog10ActCoefBDotModel( temperature,
                                ionicStrength,
                                log10PrimaryActCoeff,
                                dLog10PrimaryActCoeff_dIonicStrength,
                                log10SecActCoeff,
                                dLog10SecActCoeff_dIonicStrength );

  computeSeondarySpeciesConcAndDerivative( temperature,
                                           log10PrimaryActCoeff,
                                           dLog10PrimaryActCoeff_dIonicStrength,
                                           log10SecActCoeff,
                                           dLog10SecActCoeff_dIonicStrength,
                                           primarySpeciesConcentration,
                                           secondarySpeciesConcentration,
                                           dLog10SecConc_dLog10PrimaryConc );

  computeTotalConcAndDerivative( temperature,
                                 primarySpeciesConcentration,
                                 secondarySpeciesConcentration,
                                 dLog10SecConc_dLog10PrimaryConc,
                                 totalConcentration,
                                 dTotalConc_dLog10PrimaryConc );

  for( int i=0; i<m_numPrimarySpecies; i++ )
  {
    rhs[i] = 1 - totalConcentration[i] / primarySpeciesTotalConcentration[i];
    rhs[i] = -rhs[i];
    for( int j=0; j<m_numPrimarySpecies; j++ )
    {
      matrix( i, j ) = -dTotalConc_dLog10PrimaryConc( i, j ) / primarySpeciesTotalConcentration[i];
    }
  }
}

void EquilibriumReactions::KernelWrapper::updateConcentrations( real64 const temperature,
                                                                arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                                                arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConcentration ) const

{
  stackArray2d< real64, ReactionsBase::maxNumPrimarySpecies * ReactionsBase::maxNumPrimarySpecies > matrix( m_numPrimarySpecies, m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > rhs( m_numPrimarySpecies );
  stackArray1d< real64, ReactionsBase::maxNumPrimarySpecies > solution( m_numPrimarySpecies );

  setInitialGuess( primarySpeciesTotalConcentration, primarySpeciesConcentration );

  bool converged = false;
  for( int iteration = 0; iteration < m_maxNumIterations; iteration++ )
  {

    for( int i = 0; i< m_numPrimarySpecies; i++ )
    {
      rhs[i] = 0.0;
      solution[i] = 0.0;
      for( int j = 0; j< m_numPrimarySpecies; j++ )
      {
        matrix( i, j ) = 0.0;
      }
    }

    assembleEquilibriumReactionSystem( temperature,
                                       primarySpeciesTotalConcentration,
                                       primarySpeciesConcentration,
                                       secondarySpeciesConcentration,
                                       matrix,
                                       rhs );

    real64 const residualNorm = BlasLapackLA::vectorNorm2( rhs.toSliceConst() );

    if( residualNorm < m_newtonTol && iteration >= 1 )
    {
      converged = true;
      break;
    }

    BlasLapackLA::solveLinearSystem( matrix, rhs, solution );

    updatePrimarySpeciesConcentrations( solution, primarySpeciesConcentration );
  }
  GEOS_ERROR_IF( !converged, "Equilibrium reactions did not converge." );
}

// function to compute the derivative of the concentration of secondary species with respect to the concentration of the primary species.
void EquilibriumReactions::KernelWrapper::computeSeondarySpeciesConcAndDerivative( real64 const temperature,
                                                                                   arraySlice1d< real64 const > const & log10PrimaryActCoeff,
                                                                                   arraySlice1d< real64 const > const & dLog10PrimaryActCoeff_dIonicStrength,
                                                                                   arraySlice1d< real64 const > const & log10SecActCoeff,
                                                                                   arraySlice1d< real64 const > const & dLog10SecActCoeff_dIonicStrength,
                                                                                   arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                                   arraySlice1d< real64, compflow::USD_COMP - 1 > const & secondarySpeciesConectration,
                                                                                   arraySlice2d< real64 > const & dLog10SecConc_dLog10PrimaryConc ) const
{
  GEOS_UNUSED_VAR( temperature );

  // Compute d(concentration of dependent species)/d(concentration of basis species)
  for( int iSec = 0; iSec < m_numSecondarySpecies; iSec++ )
  {
    real64 log10SecConc = -m_log10EqConst[iSec] - log10SecActCoeff[iSec];

    for( int jPri = 0; jPri < m_numPrimarySpecies; jPri++ )
    {
      real64 const dIonicStrength_dPrimaryConc = log( 10 ) * 0.5 * m_chargePrimary[jPri] * m_chargePrimary[jPri];

      log10SecConc += m_stoichMatrix[iSec][jPri] * ( log10( primarySpeciesConcentration[jPri] ) + log10PrimaryActCoeff[jPri] );
      dLog10SecConc_dLog10PrimaryConc[iSec][jPri] += m_stoichMatrix[iSec][jPri] - dLog10SecActCoeff_dIonicStrength[iSec] * primarySpeciesConcentration[jPri] *
                                                     dIonicStrength_dPrimaryConc;
      for( int kDerivative = 0; kDerivative < m_numPrimarySpecies; kDerivative++ )
      {
        // add contribution to the derivtive from all primary activity coefficients
        dLog10SecConc_dLog10PrimaryConc[iSec][jPri] += m_stoichMatrix[iSec][kDerivative] * dLog10PrimaryActCoeff_dIonicStrength[kDerivative] * primarySpeciesConcentration[jPri] *
                                                       dIonicStrength_dPrimaryConc;
      }

    }
    secondarySpeciesConectration[iSec] = pow( 10, log10SecConc );
  }

}

void EquilibriumReactions::KernelWrapper::computeTotalConcAndDerivative( real64 const temperature,
                                                                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesConcentration,
                                                                         arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & secondarySpeciesConectration,
                                                                         arraySlice2d< real64 const > const & dLog10SecConc_dLog10PrimaryConc,
                                                                         arraySlice1d< real64 > const & totalConc,
                                                                         arraySlice2d< real64 > const & dTotalConc_dLog10PrimaryConc ) const


{
  GEOS_UNUSED_VAR( temperature );

  // This function computes the total concentration and its derivative with respect to log10(basis species concentrations).
  for( int iPri = 0; iPri < m_numPrimarySpecies; iPri++ )
  {
    totalConc[iPri] = primarySpeciesConcentration[iPri];
    // d(total concentration)/d(log10(concentration))
    dTotalConc_dLog10PrimaryConc[iPri][iPri] = log( 10.0 ) * primarySpeciesConcentration[iPri];
    // contribution from all dependent species
    for( int jSec = 0; jSec < m_numSecondarySpecies; jSec++ )
    {
      totalConc[iPri] += m_stoichMatrix[jSec][iPri] * secondarySpeciesConectration[jSec];
      for( int kDerivative = 0; kDerivative < m_numPrimarySpecies; kDerivative++ )
      {
        // add contribution to the derivtive from dependent species via the chain rule
        dTotalConc_dLog10PrimaryConc[iPri][kDerivative] += m_stoichMatrix[jSec][iPri] * log( 10.0 ) *
                                                           secondarySpeciesConectration[jSec] * dLog10SecConc_dLog10PrimaryConc[jSec][kDerivative];
      }
    }
  }
}

void EquilibriumReactions::KernelWrapper::
  updatePrimarySpeciesConcentrations( arraySlice1d< real64 const > const solution,
                                      arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration ) const
{
  for( integer i = 0; i < m_numPrimarySpecies; i++ )
  {
    primarySpeciesConcentration[i] *= pow( 10, solution[i] );
  }
}

void EquilibriumReactions::KernelWrapper::setInitialGuess( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & primarySpeciesTotalConcentration,
                                                           arraySlice1d< real64, compflow::USD_COMP - 1 > const & primarySpeciesConcentration ) const
{
  for( integer i = 0; i < m_numPrimarySpecies; i++ )
  {
    primarySpeciesConcentration[i] = primarySpeciesTotalConcentration[i];
  }
  real64 const hPlusConcentration = 2*primarySpeciesConcentration[2]-2*primarySpeciesConcentration[3]-primarySpeciesConcentration[4]+2*primarySpeciesConcentration[5]+primarySpeciesConcentration[6];
  if( hPlusConcentration < 0 )
  {
    primarySpeciesConcentration[0] = -hPlusConcentration;
  }
  else
  {
    primarySpeciesConcentration[0] = 1e-7;
  }

}


} // end namespace chemicalReactions

} // namespace constitutive

} // end namespace geos
