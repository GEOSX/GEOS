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
 * @file FenghourCO2Viscosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_FENGHOURCO2VISCOSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_FENGHOURCO2VISCOSITY_HPP_

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

class FenghourCO2ViscosityUpdate final : public PVTFunctionBaseUpdate
{
public:

  FenghourCO2ViscosityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                              TableFunction const & CO2ViscosityTable )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_CO2ViscosityTable( CO2ViscosityTable.createKernelWrapper() )
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
    m_CO2ViscosityTable.move( space, touch );
  }

protected:

  /// Table with viscosity tabulated as a function (P,T)
  TableFunction::KernelWrapper m_CO2ViscosityTable;

};

class FenghourCO2Viscosity : public PVTFunctionBase
{
public:

  FenghourCO2Viscosity( string const & name,
                        string_array const & inputParams,
                        string_array const & componentNames,
                        array1d< real64 > const & componentMolarWeight,
                        TableFunction::OutputOptions const pvtOutputOpts );

  virtual ~FenghourCO2Viscosity() override = default;

  static string catalogName() { return "FenghourCO2Viscosity"; }

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
  using KernelWrapper = FenghourCO2ViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  /// Table with CO2 viscosity tabulated as a function of (P,T)
  TableFunction const * m_CO2ViscosityTable;

};

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void FenghourCO2ViscosityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                          arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                          real64 & value,
                                          arraySlice1d< real64, USD3 > const & dValue,
                                          bool useMass ) const
{
  GEOS_UNUSED_VAR( phaseComposition,
                   dPhaseComposition,
                   useMass );

  using Deriv = constitutive::multifluid::DerivativeOffset;

  real64 const input[2] = { pressure, temperature };
  real64 densityDeriv[2]{};
  value = m_CO2ViscosityTable.compute( input, densityDeriv );

  LvArray::forValuesInSlice( dValue, []( real64 & val ){ val = 0.0; } );
  dValue[Deriv::dP] = densityDeriv[0];
  dValue[Deriv::dT] = densityDeriv[1];

}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITY_HPP_
