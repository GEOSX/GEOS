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
 * @file CO2SolubilityDuanSun.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2SolubilityDuanSun.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2EOSSolver.hpp"

#include "common/Units.hpp"

namespace geos
{
namespace constitutive
{
namespace PVTProps
{

namespace
{

inline constexpr real64 P_Pa_f = 1e+5;
inline constexpr real64 P_c    = 73.773 * P_Pa_f;
inline constexpr real64 T_c    = 304.1282;
inline constexpr real64 Rgas   = constants::gasConstant;
inline constexpr real64 V_c    = Rgas*T_c/P_c;

// these coefficients are in Table (A1) of Duan and Sun (2003)
inline constexpr real64 acoef[15] =
{ 8.99288497e-2, -4.94783127e-1, 4.77922245e-2, 1.03808883e-2, -2.82516861e-2, 9.49887563e-2, 5.20600880e-4,
  -2.93540971e-4, -1.77265112e-3, -2.51101973e-5, 8.93353441e-5, 7.88998563e-5, -1.66727022e-2, 1.398, 2.96e-2 };

real64 co2EOS( real64 const & T, real64 const & P, real64 const & V_r )
{
  // reduced pressure
  real64 const P_r = P*P_Pa_f/P_c;
  // reduced temperature
  real64 const T_r = units::convertCToK( T )/T_c;

  // CO2 equation of state
  // see equation (A1) in Duan and Sun (2003)
  real64 const f_Z = 1.0
                     + ( acoef[0] + acoef[1]/(T_r * T_r) + acoef[2]/(T_r * T_r * T_r) )/V_r
                     + ( acoef[3] + acoef[4]/(T_r * T_r) + acoef[5]/(T_r * T_r * T_r) )/(V_r*V_r)
                     + ( acoef[6] + acoef[7]/(T_r * T_r) + acoef[8]/(T_r * T_r * T_r) )/(V_r*V_r*V_r*V_r)
                     + ( acoef[9] + acoef[10]/(T_r * T_r) + acoef[11]/(T_r * T_r * T_r) )/(V_r*V_r*V_r*V_r*V_r)
                     + acoef[12]/(T_r * T_r * T_r)/(V_r * V_r) * (acoef[13] + acoef[14]/(V_r * V_r)) * exp( -acoef[14]/(V_r * V_r)) - P_r * V_r / T_r;

  return f_Z;
}

real64 PWater( real64 const & T )
{
  // these coefficients are defined in Table (B1) of Duan and Sun (2003)
  static constexpr real64 ccoef[5] = { -38.640844, 5.8948420, 59.876516, 26.654627, 10.637097 };

  // H2O critical pressure (bars)
  real64 const P_c_w = 220.85;
  // H2O critical temperature (K)
  real64 const T_c_w = 647.29;
  real64 const tt = ( units::convertCToK( T )-T_c_w )/T_c_w;
  // Empirical model for water pressure of equation (B1) of Duan and Sun (2003)
  real64 const x = ( P_c_w*units::convertCToK( T )/T_c_w )
                   * (1
                      + ccoef[0]*pow( -tt, 1.9 )
                      + ccoef[1]*tt
                      + ccoef[2]*tt*tt
                      + ccoef[3]*tt*tt*tt
                      + ccoef[4]*tt*tt*tt*tt);

  return x;
}

real64 logF( real64 const & T, real64 const & P, real64 const & V_r )
{
  // reduced pressure
  real64 const P_r = P*P_Pa_f/P_c;
  // reduced temperature
  real64 const T_r = units::convertCToK( T ) / T_c;
  real64 const Z   = P_r * V_r/T_r;

  real64 const inv_T_r = 1.0/T_r;
  real64 const inv_T_r2 = inv_T_r*inv_T_r;
  real64 const inv_T_r3 = inv_T_r2*inv_T_r;
  real64 const inv_V_r = 1.0/V_r;
  real64 const inv_V_r2 = inv_V_r*inv_V_r;
  real64 const inv_V_r3 = inv_V_r2*inv_V_r;
  real64 const inv_V_r4 = inv_V_r3*inv_V_r;
  real64 const inv_V_r5 = inv_V_r4*inv_V_r;

  // fugacity coefficient of CO2, equation (A6) of Duan and Sun (2003)
  real64 const log_f = Z - 1 - log( Z ) +
                       ( acoef[0] + acoef[1]*inv_T_r2 + acoef[2]*inv_T_r3 )*inv_V_r
                       + ( acoef[3] + acoef[4]*inv_T_r2 + acoef[5]*inv_T_r3 )*0.5*inv_V_r2
                       + ( acoef[6] + acoef[7]*inv_T_r2 + acoef[8]*inv_T_r3 )*0.25*inv_V_r4
                       + ( acoef[9] + acoef[10]*inv_T_r2 + acoef[11]*inv_T_r3 )*0.2*inv_V_r5
                       + acoef[12]*0.5*inv_T_r3/acoef[14] * ( acoef[13] + 1.0 - (acoef[13] + 1.0 + acoef[14]*inv_V_r2) * exp( -acoef[14]*inv_V_r2 ) );
  //This causes a divide by zero FPE when using clang14 on ruby
  // real64 const log_f = Z - 1 - log( Z ) +
  //                      ( acoef[0] + acoef[1]/T_r/T_r + acoef[2]/T_r/T_r/T_r )/V_r
  //                      + ( acoef[3] + acoef[4]/T_r/T_r + acoef[5]/T_r/T_r/T_r )/2.0/V_r/V_r
  //                      + ( acoef[6] + acoef[7]/T_r/T_r + acoef[8]/T_r/T_r/T_r )/4.0/V_r/V_r/V_r/V_r
  //                      + ( acoef[9] + acoef[10]/T_r/T_r + acoef[11]/T_r/T_r/T_r )/5.0/V_r/V_r/V_r/V_r/V_r
  //                      + acoef[12]/2.0/T_r/T_r/T_r/acoef[14] * ( acoef[13] + 1.0 - (acoef[13] + 1.0 + acoef[14]/V_r/V_r) * exp(
  // -acoef[14]/V_r/V_r ) );


  return log_f;
}

real64 Par( real64 const & T, real64 const & P, real64 const (&cc)[11] )
{
  // "equation for the parameters", see equation (7) of Duan and Sun (2003)
  real64 x = cc[0]
             + cc[1]*T
             + cc[2]/T
             + cc[3]*T*T
             + cc[4]/(630.0-T)
             + cc[5]*P
             + cc[6]*P *log( T )
             + cc[7]*P/T
             + cc[8]*P/(630.0-T)
             + cc[9]*P*P/(630.0-T)/(630.0-T)
             + cc[10]*T *log( P );

  return x;
}

real64 CO2SolubilityFunction( string const & name,
                              real64 const & tolerance,
                              real64 const & T,
                              real64 const & P,
                              real64 (* f)( real64 const & x1, real64 const & x2, real64 const & x3 ) )
{
  // compute the initial guess for Newton's method
  real64 const initialReducedVolume = 0.75*Rgas*units::convertCToK( T )/(P*P_Pa_f)*(1/V_c);

  // define the local solver parameters
  // for now, this is hard-coded, but we may want to let the user access the parameters at some point
  integer const maxNumNewtonIter = 500;
  integer const maxNumBacktrackIter = 8;
  real64 const maxAbsUpdate = 1e12;
  real64 const minAbsDeriv = 0;
  real64 const allowedMinValue = 0.05; // value chosen to match previous implementation
  real64 const presMultiplierForReporting = 1e5; // this is because P is in hectopascal in this function

  // solve the CO2 equation of state for this pair of (pres, temp)
  // return the reduced volume
  return CO2EOSSolver::solve( name,
                              maxNumNewtonIter,
                              maxNumBacktrackIter,
                              tolerance,
                              minAbsDeriv,
                              maxAbsUpdate,
                              allowedMinValue,
                              initialReducedVolume,
                              T,
                              P,
                              presMultiplierForReporting,
                              f );
}

void calculateCO2Solubility( string const & functionName,
                             real64 const & tolerance,
                             PTTableCoordinates const & tableCoords,
                             real64 const & salinity,
                             array1d< real64 > const & values )
{
  // Interaction parameters, see Table 2 of Duan and Sun (2003)
  static constexpr real64 mu[11] =
  { 28.9447706, -0.0354581768, -4770.67077, 1.02782768e-5, 33.8126098, 9.04037140e-3,
    -1.14934031e-3, -0.307405726, -0.0907301486, 9.32713393e-4, 0 };
  static constexpr real64 lambda[11] = { -0.411370585, 6.07632013e-4, 97.5347708, 0, 0, 0, 0, -0.0237622469, 0.0170656236, 0, 1.41335834e-5 };
  static constexpr real64 zeta[11] = { 3.36389723e-4, -1.98298980e-5, 0, 0, 0, 0, 0, 2.12220830e-3, -5.24873303e-3, 0, 0 };

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const P = tableCoords.getPressure( i ) / P_Pa_f;

    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const T = tableCoords.getTemperature( j );

      // compute reduced volume by solving the CO2 equation of state
      real64 const V_r = CO2SolubilityFunction( functionName, tolerance, T, P, &co2EOS );

      // compute equation (6) of Duan and Sun (2003)
      real64 const logK = Par( units::convertCToK( T ), P, mu )
                          - logF( T, P, V_r )
                          + 2*Par( units::convertCToK( T ), P, lambda ) * salinity
                          + Par( units::convertCToK( T ), P, zeta ) * salinity * salinity;
      real64 const expLogK = exp( logK );

      // mole fraction of CO2 in vapor phase, equation (4) of Duan and Sun (2003)
      real64 const Pw = PWater( T );
      real64 const y_CO2 = (P - Pw)/P;
      values[j*nPressures+i] = y_CO2 * P / expLogK;

      GEOS_WARNING_IF( expLogK <= 1e-10,
                       GEOS_FMT( "CO2Solubility: exp(logK) = {} is too small (logK = {}, P = {}, T = {}, V_r = {}), resulting solubility value is {}",
                                 expLogK, logK, P, T, V_r, values[j*nPressures+i] ));
    }
  }
}

} // end namespace

void CO2SolubilityDuanSun::populateSolubilityTables( string const & functionName,
                                                     PTTableCoordinates const & tableCoords,
                                                     real64 const & salinity,
                                                     real64 const & tolerance,
                                                     array1d< real64 > const & co2SolubilityValues,
                                                     array1d< real64 > const & h2oSolubilityValues )
{
  h2oSolubilityValues.zero();
  calculateCO2Solubility( functionName,
                          tolerance,
                          tableCoords,
                          salinity,
                          co2SolubilityValues );
}

} // end namespace PVTProps
} // namespace constitutive
} // end namespace geos
