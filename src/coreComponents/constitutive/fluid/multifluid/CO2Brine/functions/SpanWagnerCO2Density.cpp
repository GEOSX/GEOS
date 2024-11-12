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
 * @file SpanWagnerCO2Density.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/SpanWagnerCO2Density.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/CO2EOSSolver.hpp"
#include "functions/FunctionManager.hpp"
#include "common/Units.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

namespace
{

real64 co2HelmholtzEnergy( real64 const & T,
                           real64 const & P,
                           real64 const & rho )
{
  // all the coefficients below are defined in Table (31) of Span and Wagner (1996)
  // the variable names used below follow the notation of the Table
  constexpr real64 n[] =
  {0.38856823203161, 2.938547594274, -5.5867188534934, -0.76753199592477, 0.31729005580416, 0.54803315897767, 0.12279411220335, 2.165896154322,
   1.5841735109724, -0.23132705405503, 0.058116916431436, -0.55369137205382, 0.48946615909422, -0.024275739843501, 0.062494790501678, -0.12175860225246,
   -0.37055685270086, -0.016775879700426, -0.11960736637987, -0.045619362508778, 0.035612789270346, -0.0074427727132052, -0.0017395704902432,
   -0.021810121289527, 0.024332166559236, -0.037440133423463, 0.14338715756878, -0.13491969083286, -0.02315122505348, 0.012363125492901, 0.002105832197294,
   -0.00033958519026368, 0.0055993651771592, -0.00030335118055646, -213.65488688320, 26641.569149272, -24027.212204557, -283.41603423999, 212.47284400179,
   -0.66642276540751, 0.72608632349897, 0.055068668612842};

  constexpr real64 d[] = {1, 1, 1, 1, 2, 2, 3, 1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5, 6, 7, 8, 10, 4, 8, 2, 2, 2, 3, 3};

  constexpr real64 t[] =
  {0, 0.75, 1, 2, 0.75, 2, 0.75, 1.5, 1.5, 2.5, 0, 1.5, 2, 0, 1, 2, 3, 6, 3, 6, 8, 6, 0, 7, 12, 16, 22, 24, 16, 24, 8, 2, 28, 14, 1, 0, 1, 3, 3};

  constexpr real64 c[] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 6};

  constexpr real64 alpha[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 25, 25, 15, 20};

  constexpr real64 beta[] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 325, 300, 300, 275, 275, 0.3, 0.3, 0.3};

  constexpr real64 gamma0[] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.16, 1.19, 1.19, 1.25, 1.22};

  constexpr real64 epslon[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1., 1., 1., 1., 1.};

  constexpr real64 a[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.5, 3.5, 3.};

  constexpr real64 b[] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.875, 0.925, 0.875};

  constexpr real64 A[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.7, 0.7, 0.7};

  constexpr real64 B[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.3, 0.3, 1.0};

  constexpr real64 C[] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10.0, 10.0, 12.5};

  constexpr real64 D[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 275, 275, 275};

  constexpr real64 T_c = 304.1282;
  constexpr real64 rho_c = 467.6;

  constexpr real64 R = 188.9241;

  real64 const tau = T_c / T;
  real64 const delta = rho / rho_c;

  real64 theta, Delta, Phi, dPhi_delta, dDelta_delta, dDelta_delta_b;

  // residual part of the energy equation
  // this is defined in Table (32) of Span and Wagner (1996)
  real64 phi_r_delta = 0.0;

  for( int i=0; i<7; i++ )
  {
    phi_r_delta += n[i] * d[i] * pow( delta, d[i]-1.0 ) * pow( tau, t[i] );
  }

  for( int i=7; i<34; i++ )
  {
    phi_r_delta += n[i] * exp( -pow( delta, c[i] )) * pow( delta, d[i] - 1.0 ) * pow( tau, t[i] ) * (d[i] - c[i] * pow( delta, c[i] ));
  }

  for( int i=34; i<39; i++ )
  {
    phi_r_delta += n[i] *
                   pow( delta,
                        d[i] ) *
                   pow( tau,
                        t[i] ) *
                   exp( -alpha[i] *
                        pow( delta - epslon[i], 2.0 ) - beta[i] * pow( tau - gamma0[i], 2.0 )) * (d[i] / delta - 2.0 * alpha[i] * (delta - epslon[i]));
  }

  for( int i=39; i<42; i++ )
  {
    theta = 1.0 - tau + A[i] * pow( pow( delta - 1.0, 2.0 ), 1.0 / (2.0 * beta[i]));
    Delta = pow( theta, 2.0 ) + B[i] * pow( pow( delta - 1.0, 2.0 ), a[i] );

    Phi = exp( -C[i] * pow( delta - 1.0, 2.0 ) - D[i] * pow( tau - 1.0, 2.0 ));

    dPhi_delta = -2.0 * C[i] * (delta - 1.0) * Phi;

    dDelta_delta = (delta - 1.0) *
                   (A[i] * theta * 2.0 / beta[i] *
                    pow( pow( delta - 1.0, 2.0 ), 1.0 / (2.0 * beta[i]) - 1.0 ) + 2.0 * B[i] * a[i] * pow( pow( delta - 1.0, 2.0 ), a[i] - 1.0 ));

    dDelta_delta_b = b[i] * pow( Delta, b[i] - 1.0 ) * dDelta_delta;

    phi_r_delta += n[i] * (pow( Delta, b[i] ) * (Phi + delta *dPhi_delta) + dDelta_delta_b * delta * Phi);
  }

  // Helmholtz energy equation: see equation (2.2) from Span and Wagner (1996)
  return rho * (1.0 + delta * phi_r_delta) / (P / (R * T)) - 1.0;
}

real64 spanWagnerCO2DensityFunction( string const & name,
                                     real64 const & tolerance,
                                     real64 const & T,
                                     real64 const & P,
                                     real64 (* f)( real64 const & x1, real64 const & x2, real64 const & x3 ) )
{
  // define local variables needed to compute the initial guess
  constexpr real64 P_Pa_f = 1e+5;
  constexpr real64 P_c = 73.773 * P_Pa_f;
  constexpr real64 T_c = 304.1282;
  constexpr real64 rho_c = 467.6;
  constexpr real64 R = 188.9241;

  constexpr real64 vpa[] = { -7.0602087, 1.9391218, -1.6463597, -3.2995634 };
  constexpr real64 vpt[] = { 1.0, 1.5, 2.0, 4.0 };

  constexpr real64 lda[] = { 1.9245108, -0.62385555, -0.32731127, 0.39245142 };
  constexpr real64 ldt[] = { 0.340, 0.5, 1.6666666667, 1.833333333 };

  // compute the initial guess for density using the method ported from NUFT
  real64 initialRho = 0.0;
  if( T < T_c )
  {
    real64 sum = 0;
    for( integer i=0; i < 4; i++ )
    {
      sum += vpa[i] * pow( 1.0 - T / T_c, vpt[i] );
    }
    real64 const psat = P_c * exp( T_c / T * sum );

    sum = 0;
    for( integer i=0; i < 4; i++ )
    {
      sum += lda[i] * pow( 1.0 - T / T_c, ldt[i] );
    }
    real64 const denl = rho_c * exp( sum );
    initialRho = ( P >= psat ) ?  denl : P / (R * T);
  }
  else
  {
    initialRho = 1.0 / (30.0 + R * T / P);
  }

  // define the local solver parameters
  // for now, this is hard-coded, but we may want to let the user access the parameters at some point
  integer const maxNumNewtonIter = 500;
  integer const maxNumBacktrackIter = 4;
  real64 const maxAbsUpdate = 500;
  real64 const minAbsDeriv = 1e-12; // needed to avoid divisions by zero in the calculation of the Newton update
  real64 const allowedMinValue = -1e12; // use this negative value to disable the chopping on the min value
  real64 const presMultiplierForReporting = 1;

  // solve the Helmholtz energy equation for this pair of (pres, temp)
  // return the density
  return CO2EOSSolver::solve( name,
                              maxNumNewtonIter,
                              maxNumBacktrackIter,
                              tolerance,
                              minAbsDeriv,
                              maxAbsUpdate,
                              allowedMinValue,
                              initialRho,
                              T,
                              P,
                              presMultiplierForReporting,
                              f );
}


TableFunction const * makeDensityTable( string_array const & inputParams,
                                        string const & functionName,
                                        FunctionManager & functionManager )
{
  string const tableName = functionName + "_table";

  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    TableFunction * const densityTable = functionManager.getGroupPointer< TableFunction >( tableName );
    densityTable->initializeFunction();
    densityTable->setDimUnits( PTTableCoordinates::coordsUnits );
    densityTable->setValueUnits( units::Density );
    return densityTable;
  }
  else
  {
    PTTableCoordinates tableCoords;
    PVTFunctionHelpers::initializePropertyTable( inputParams, tableCoords );

    real64 tolerance = 1e-10;
    try
    {
      if( inputParams.size() >= 9 )
      {
        tolerance = stod( inputParams[8] );
      }
    }
    catch( const std::invalid_argument & e )
    {
      GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value: {}", functionName, e.what() ), InputError );
    }

    array1d< real64 > densities( tableCoords.nPressures() * tableCoords.nTemperatures() );
    SpanWagnerCO2Density::calculateCO2Density( functionName, tolerance, tableCoords, densities );

    TableFunction * const densityTable = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );
    densityTable->setTableCoordinates( tableCoords.getCoords(), tableCoords.coordsUnits );
    densityTable->setTableValues( densities, units::Density );
    densityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return densityTable;
  }
}

} // namespace

void SpanWagnerCO2Density::calculateCO2Density( string const & functionName,
                                                real64 const & tolerance,
                                                PTTableCoordinates const & tableCoords,
                                                array1d< real64 > const & densities )
{

  constexpr real64 TK_f = constants::zeroDegreesCelsiusInKelvin;

  localIndex const nPressures = tableCoords.nPressures();
  localIndex const nTemperatures = tableCoords.nTemperatures();

  for( localIndex i = 0; i < nPressures; ++i )
  {
    real64 const PPa = tableCoords.getPressure( i );
    for( localIndex j = 0; j < nTemperatures; ++j )
    {
      real64 const TK = tableCoords.getTemperature( j ) + TK_f;
      densities[j*nPressures+i] = spanWagnerCO2DensityFunction( functionName, tolerance, TK, PPa, &co2HelmholtzEnergy );
    }
  }
}

SpanWagnerCO2Density::SpanWagnerCO2Density( string const & name,
                                            string_array const & inputParams,
                                            string_array const & componentNames,
                                            array1d< real64 > const & componentMolarWeight,
                                            bool const printTable ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  string const expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames, "componentNames" );

  m_CO2DensityTable = makeDensityTable( inputParams, m_functionName, FunctionManager::getInstance() );
  if( printTable )
    m_CO2DensityTable->print( m_CO2DensityTable->getName() );
}

void SpanWagnerCO2Density::checkTablesParameters( real64 const pressure,
                                                  real64 const temperature ) const
{
  m_CO2DensityTable->checkCoord( pressure, 0 );
  m_CO2DensityTable->checkCoord( temperature, 1 );
}

SpanWagnerCO2Density::KernelWrapper
SpanWagnerCO2Density::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_CO2DensityTable,
                        m_CO2Index );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, SpanWagnerCO2Density, string const &, string_array const &, string_array const &, array1d< real64 > const &, bool const )

} // namespace PVTProps

} // namespace constitutive

} // namespace geos
