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
 * @file CO2Enthalpy.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2ENTHALPY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_CO2ENTHALPY_HPP_

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

class CO2EnthalpyUpdate final : public PVTFunctionBaseUpdate
{
public:

  CO2EnthalpyUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                     TableFunction const & CO2EnthalpyTable,
                     integer const CO2Index )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_CO2EnthalpyTable( CO2EnthalpyTable.createKernelWrapper() ),
    m_CO2Index( CO2Index )
  {}

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice2d< real64 const, USD2 > const & dPhaseComposition_dGlobalCompFraction,
                real64 & value,
                arraySlice1d< real64, USD3 > const & dValue,
                bool useMass ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    PVTFunctionBaseUpdate::move( space, touch );
    m_CO2EnthalpyTable.move( space, touch );
  }

protected:

  /// Table with CO2 enthalpy tabulated as a function of (P,T)
  TableFunction::KernelWrapper m_CO2EnthalpyTable;

  /// Index of the CO2 component
  integer m_CO2Index;
};

class CO2Enthalpy : public PVTFunctionBase
{
public:

  CO2Enthalpy( string const & name,
               string_array const & inputParams,
               string_array const & componentNames,
               array1d< real64 > const & componentMolarWeight,
               bool const printTable );

  static string catalogName() { return "CO2Enthalpy"; }

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
  using KernelWrapper = CO2EnthalpyUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  static void calculateCO2Enthalpy( PTTableCoordinates const & tableCoords,
                                    array1d< real64 > const & densities,
                                    array1d< real64 > const & enthalpies );


private:

  /// Table with CO2 enthalpy tabulated as a function of (P,T)
  TableFunction const * m_CO2EnthalpyTable;

  /// Index of the CO2 phase
  integer m_CO2Index;
};

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
void CO2EnthalpyUpdate::compute( real64 const & pressure,
                                 real64 const & temperature,
                                 arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                 arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                 real64 & value,
                                 arraySlice1d< real64, USD3 > const & dValue,
                                 bool useMass ) const
{
  GEOS_UNUSED_VAR( phaseComposition, dPhaseComposition );

  using Deriv = multifluid::DerivativeOffset;

  real64 const input[2] = { pressure, temperature };
  real64 CO2EnthalpyDeriv[2]{};

  value = m_CO2EnthalpyTable.compute( input, CO2EnthalpyDeriv );

  LvArray::forValuesInSlice( dValue, []( real64 & val ){ val = 0.0; } );
  dValue[Deriv::dP] = CO2EnthalpyDeriv[0];
  dValue[Deriv::dT] = CO2EnthalpyDeriv[1];

  if( !useMass )
  {
    real64 const CO2MWInv = 1.0 / m_componentMolarWeight[m_CO2Index];
    value *= CO2MWInv;
    dValue[Deriv::dP] *= CO2MWInv;
    dValue[Deriv::dT] *= CO2MWInv;
  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2ENTHALPY_HPP_
