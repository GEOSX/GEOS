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
 * @file BrineEnthalpy.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_BRINEENTHALPY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_BRINEENTHALPY_HPP_

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

class BrineEnthalpyUpdate final : public PVTFunctionBaseUpdate
{
public:

  BrineEnthalpyUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                       TableFunction const & CO2EnthalpyTable,
                       TableFunction const & brineEnthalpyTable,
                       integer const CO2Index,
                       integer const waterIndex )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_CO2EnthalpyTable( CO2EnthalpyTable.createKernelWrapper() ),
    m_brineEnthalpyTable( brineEnthalpyTable.createKernelWrapper() ),
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
    m_CO2EnthalpyTable.move( space, touch );
    m_brineEnthalpyTable.move( space, touch );
  }

protected:

  /// Table with CO2 enthalpy tabulated as a function of (P,T)
  TableFunction::KernelWrapper m_CO2EnthalpyTable;

  /// Table with brine enthalpy tabulated as a function of (T)
  TableFunction::KernelWrapper m_brineEnthalpyTable;

  /// Index of the CO2 component
  integer m_CO2Index;

  /// Index of the water component
  integer m_waterIndex;

};

class BrineEnthalpy : public PVTFunctionBase
{
public:

  BrineEnthalpy( string const & name,
                 string_array const & inputParams,
                 string_array const & componentNames,
                 array1d< real64 > const & componentMolarWeight,
                 bool const printTable );

  static string catalogName() { return "BrineEnthalpy"; }

  virtual string getCatalogName() const final { return catalogName(); }

  /**
   * @copydoc PVTFunctionBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  void checkTablesParameters( real64 pressure, real64 temperature ) const override final;

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::ENTHALPY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrineEnthalpyUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;


private:

  /// Table with CO2 enthalpy tabulated as a function of (P,T)
  TableFunction const * m_CO2EnthalpyTable;

  /// Table with brine enthalpy tabulated as a function of (T)
  TableFunction const * m_brineEnthalpyTable;

  /// Index of the CO2 phase
  integer m_CO2Index;

  /// Index of the water phase
  integer m_waterIndex;

};

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
void BrineEnthalpyUpdate::compute( real64 const & pressure,
                                   real64 const & temperature,
                                   arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                   arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                   real64 & value,
                                   arraySlice1d< real64, USD3 > const & dValue,
                                   bool useMass ) const
{
  using Deriv = multifluid::DerivativeOffset;

  real64 const input[2] = { pressure, temperature };
  real64 brineEnthalpy_dTemperature = 0.0;
  real64 dvalue_dC = 0.0;
  real64 CO2EnthalpyDeriv[2]{};

  real64 const brineEnthalpy = m_brineEnthalpyTable.compute( &temperature, &brineEnthalpy_dTemperature );
  real64 const CO2Enthalpy = m_CO2EnthalpyTable.compute( input, CO2EnthalpyDeriv );

  //assume there are only CO2 and brine here.

  real64 const C = phaseComposition[m_waterIndex];

  if( useMass )
  {
    value = (1.0 - C ) * CO2Enthalpy
            + C * brineEnthalpy;

    dvalue_dC = brineEnthalpy - CO2Enthalpy;
    dValue[Deriv::dP] = (1.0 - C ) * CO2EnthalpyDeriv[0]
                        + dvalue_dC * dPhaseComposition[m_waterIndex][Deriv::dP];
    dValue[Deriv::dT] =
      LvArray::math::max( 0.0,
                          (1.0 - C ) * CO2EnthalpyDeriv[1]
                          + C * brineEnthalpy_dTemperature
                          + dvalue_dC * dPhaseComposition[m_waterIndex][Deriv::dT] );
  }
  else
  {
    real64 const waterMWInv = 1.0 / m_componentMolarWeight[m_waterIndex];
    real64 const CO2MWInv   = 1.0 / m_componentMolarWeight[m_CO2Index];

    value = (1.0 - C ) * CO2Enthalpy * CO2MWInv + C * brineEnthalpy * waterMWInv;

    dvalue_dC = brineEnthalpy * waterMWInv - CO2Enthalpy * CO2MWInv;
    dValue[Deriv::dP] = (1.0 - C ) * CO2EnthalpyDeriv[0] * CO2MWInv
                        + dvalue_dC * dPhaseComposition[m_waterIndex][Deriv::dP];
    dValue[Deriv::dT] =
      LvArray::math::max( 0.0,
                          (1.0 - C ) * CO2EnthalpyDeriv[1] * CO2MWInv
                          + C * brineEnthalpy_dTemperature * waterMWInv
                          + dvalue_dC * dPhaseComposition[m_waterIndex][Deriv::dT] );
  }

  dValue[Deriv::dC+m_CO2Index]   = dvalue_dC * dPhaseComposition[m_waterIndex][Deriv::dC+m_CO2Index];
  dValue[Deriv::dC+m_waterIndex] = dvalue_dC * dPhaseComposition[m_waterIndex][Deriv::dC+m_waterIndex];
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEENTHALPY_HPP_
