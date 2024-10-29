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
 * @file PhillipsBrineViscosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_PHILLIPSBRINEVISCOSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_PHILLIPSBRINEVISCOSITY_HPP_

#include "PVTFunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

class PhillipsBrineViscosityUpdate final : public PVTFunctionBaseUpdate
{
public:

  PhillipsBrineViscosityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                                TableFunction const & waterViscosityTable,
                                real64 const coef0,
                                real64 const coef1 )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_waterViscosityTable( waterViscosityTable.createKernelWrapper() ),
    m_coef0( coef0 ),
    m_coef1( coef1 )
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
    m_waterViscosityTable.move( space, touch );
  }

protected:

  /// Table with water viscosity tabulated as a function of temperature
  TableFunction::KernelWrapper m_waterViscosityTable;

  real64 m_coef0;

  real64 m_coef1;

};


class PhillipsBrineViscosity : public PVTFunctionBase
{
public:

  PhillipsBrineViscosity( string const & name,
                          string_array const & inputPara,
                          string_array const & componentNames,
                          array1d< real64 > const & componentMolarWeight,
                          bool const printTable );

  virtual ~PhillipsBrineViscosity() override = default;

  static string catalogName() { return "PhillipsBrineViscosity"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  /**
   * @copydoc PVTFunctionBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  void checkTablesParameters( real64 pressure, real64 temperature ) const override final;

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PhillipsBrineViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  void makeCoefficients( string_array const & inputPara );

  /// Table with water viscosity tabulated as a function (T)
  TableFunction const * m_waterViscosityTable;

  real64 m_coef0;

  real64 m_coef1;

};

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PhillipsBrineViscosityUpdate::compute( real64 const & pressure,
                                            real64 const & temperature,
                                            arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                            arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                            real64 & value,
                                            arraySlice1d< real64, USD3 > const & dValue,
                                            bool useMass ) const
{
  GEOS_UNUSED_VAR( pressure,
                   phaseComposition,
                   dPhaseComposition,
                   useMass );

  using Deriv = multifluid::DerivativeOffset;

  // compute the viscosity of pure water as a function of temperature
  real64 dPureWaterVisc_dTemperature;
  real64 const pureWaterVisc = m_waterViscosityTable.compute( &temperature, &dPureWaterVisc_dTemperature );

  // then compute the brine viscosity, accounting for the presence of salt

  real64 const viscMultiplier = m_coef0 + m_coef1 * temperature;
  real64 const dViscMultiplier_dTemperature = m_coef1;
  value = pureWaterVisc * viscMultiplier;
  LvArray::forValuesInSlice( dValue, []( real64 & val ){ val = 0.0; } );
  dValue[Deriv::dP] = 0.0;
  dValue[Deriv::dT] = dPureWaterVisc_dTemperature * viscMultiplier + pureWaterVisc * dViscMultiplier_dTemperature;

}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PHILLIPSBRINEVISCOSITY_HPP_
