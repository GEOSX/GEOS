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
 * @file PhillipsBrineDensity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_SINGLEFLUID_FUNCTIONS_BRINECO2DENSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_SINGLEFLUID_FUNCTIONS_BRINECO2DENSITY_HPP_

#include "PVTFunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/functions/PVTFunctionHelpers.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

class PhillipsBrineDensityUpdate final : public PVTFunctionBaseUpdate
{
public:

  PhillipsBrineDensityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                              TableFunction const & brineDensityTable,
                              integer const CO2Index,
                              integer const waterIndex )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_brineDensityTable( brineDensityTable.createKernelWrapper() ),
    m_CO2Index( CO2Index ),
    m_waterIndex( waterIndex )
  {}

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                real64 & value,
                arraySlice1d< real64, USD3 > const & dValue,
                bool useMass ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    PVTFunctionBaseUpdate::move( space, touch );
    m_brineDensityTable.move( space, touch );
  }

protected:

  /// Table with brine density tabulated as a function (P,T,sal)
  TableFunction::KernelWrapper m_brineDensityTable;

  /// Index of the CO2 component
  integer m_CO2Index;

  /// Index of the water component
  integer m_waterIndex;

};

class PhillipsBrineDensity : public PVTFunctionBase
{
public:

  PhillipsBrineDensity( string const & name,
                        string_array const & inputParams,
                        string_array const & componentNames,
                        array1d< real64 > const & componentMolarWeight,
                        TableFunction::OutputOptions const pvtOutputOpts );

  static string catalogName() { return "PhillipsBrineDensity"; }

  virtual string getCatalogName() const final { return catalogName(); }

  /**
   * @copydoc PVTFunctionBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  void checkTablesParameters( real64 pressure, real64 temperature ) const override final;

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PhillipsBrineDensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;


private:

  /// Table with brine density tabulated as a function of (P,T,sal)
  TableFunction const * m_brineDensityTable;

  /// Index of the CO2 phase
  integer m_CO2Index;

  /// Index of the water phase
  integer m_waterIndex;

};

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PhillipsBrineDensityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                          real64 & value,
                                          arraySlice1d< real64, USD3 > const & dValue,
                                          bool useMass ) const
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

  // this method implements the method proposed by E. Garcia (2001)

  // these coefficients come from equation (2) from Garcia (2001)
  constexpr real64 a = 37.51;
  constexpr real64 b = -9.585e-2;
  constexpr real64 c = 8.740e-4;
  constexpr real64 d = -5.044e-7;

  real64 const input[2] = { pressure, temperature };
  real64 densityDeriv[2]{};
  real64 const density = m_brineDensityTable.compute( input, densityDeriv );

  // equation (2) from Garcia (2001)
  real64 const squaredTemp = temperature * temperature;
  real64 const V = (  a
                      + b * temperature
                      + c * squaredTemp
                      + d * squaredTemp * temperature ) * 1e-6;
  real64 const dV_dTemp = ( b
                            + 2 * c * temperature
                            + 3 * d * squaredTemp ) * 1e-6;

  // CO2 concentration
  real64 const wMwInv = 1.0 / m_componentMolarWeight[m_waterIndex];
  real64 const oneMinusCO2PhaseCompInv = 1.0 / ( 1.0 - phaseComposition[m_CO2Index] );
  real64 const oneMinusCO2PhaseCompInvSquared = oneMinusCO2PhaseCompInv * oneMinusCO2PhaseCompInv;
  real64 const coef = wMwInv * phaseComposition[m_CO2Index] * oneMinusCO2PhaseCompInv;
  real64 const dCoef_dPres = wMwInv * dPhaseComposition[m_CO2Index][Deriv::dP] * oneMinusCO2PhaseCompInvSquared;
  real64 const dCoef_dTemp = wMwInv * dPhaseComposition[m_CO2Index][Deriv::dT] * oneMinusCO2PhaseCompInvSquared;
  real64 dCoef_dComp[2]{};
  dCoef_dComp[m_CO2Index] = wMwInv * dPhaseComposition[m_CO2Index][Deriv::dC+m_CO2Index] * oneMinusCO2PhaseCompInvSquared;
  dCoef_dComp[m_waterIndex] = wMwInv * dPhaseComposition[m_CO2Index][Deriv::dC+m_waterIndex] * oneMinusCO2PhaseCompInvSquared;

  real64 const conc = coef * density;
  real64 const dConc_dPres = dCoef_dPres * density + coef * densityDeriv[0];
  real64 const dConc_dTemp = dCoef_dTemp * density + coef * densityDeriv[1];
  real64 dConc_dComp[2]{};
  dConc_dComp[m_CO2Index] = dCoef_dComp[m_CO2Index] * density;
  dConc_dComp[m_waterIndex] = dCoef_dComp[m_waterIndex] * density;

  // CO2 concentration times density times vol
  real64 const concDensVol = conc * density * V;
  real64 const dConcDensVol_dPres = ( dConc_dPres * density + conc * densityDeriv[0] ) * V;
  real64 const dConcDensVol_dTemp = ( dConc_dTemp * density + conc * densityDeriv[1] ) * V
                                    + conc * density * dV_dTemp;
  real64 dConcDensVol_dComp[2]{};
  dConcDensVol_dComp[m_CO2Index] = dConc_dComp[m_CO2Index] * density * V;
  dConcDensVol_dComp[m_waterIndex] = dConc_dComp[m_waterIndex] * density * V;

  // Brine density
  // equation (1) from Garcia (2001)
  value = density
          + m_componentMolarWeight[m_CO2Index] * conc
          - concDensVol;
  dValue[Deriv::dP] = densityDeriv[0]
                      + m_componentMolarWeight[m_CO2Index] * dConc_dPres
                      - dConcDensVol_dPres;
  dValue[Deriv::dT] = densityDeriv[1]
                      + m_componentMolarWeight[m_CO2Index] * dConc_dTemp
                      - dConcDensVol_dTemp;
  dValue[Deriv::dC+m_CO2Index] = m_componentMolarWeight[m_CO2Index] * dConc_dComp[m_CO2Index]
                                 - dConcDensVol_dComp[m_CO2Index];
  dValue[Deriv::dC+m_waterIndex] = m_componentMolarWeight[m_CO2Index] * dConc_dComp[m_waterIndex]
                                   - dConcDensVol_dComp[m_waterIndex];
  if( !useMass )
  {
    divideByPhaseMolarWeight( phaseComposition, dPhaseComposition, value, dValue );
  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINECO2DENSITY_HPP_
