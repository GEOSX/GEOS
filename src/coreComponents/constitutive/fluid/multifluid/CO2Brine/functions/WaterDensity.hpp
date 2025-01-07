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
 * @file WaterDensity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_WATERDENSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_WATERDENSITY_HPP_

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

class WaterDensityUpdate final : public PVTFunctionBaseUpdate
{
public:

  WaterDensityUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                      TableFunction const & waterDensityTable )
    : PVTFunctionBaseUpdate( componentMolarWeight ),
    m_waterDensityTable( waterDensityTable.createKernelWrapper() )
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
    m_waterDensityTable.move( space, touch );
  }

protected:

  /// Table with brine density tabulated as a function (P,T,sal)
  TableFunction::KernelWrapper m_waterDensityTable;

};

class WaterDensity : public PVTFunctionBase
{
public:

  WaterDensity( string const & name,
                string_array const & inputParams,
                string_array const & componentNames,
                array1d< real64 > const & componentMolarWeight,
                TableFunction::OutputOptions const pvtOutputOpts );

  static string catalogName() { return "WaterDensity"; }
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
  using KernelWrapper = WaterDensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  /// Table with brine density tabulated as a function of (P,T,sal)
  TableFunction const * m_waterDensityTable;
};

template< int USD1, int USD2, int USD3 >
GEOS_HOST_DEVICE
void WaterDensityUpdate::compute( real64 const & pressure,
                                  real64 const & temperature,
                                  arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                  arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                  real64 & value,
                                  arraySlice1d< real64, USD3 > const & dValue,
                                  bool useMass ) const
{
  GEOS_UNUSED_VAR( phaseComposition, dPhaseComposition, useMass );

  using Deriv = constitutive::multifluid::DerivativeOffset;

  real64 const input[2] = { pressure, temperature };
  real64 densityDeriv[2]{};
  value = m_waterDensityTable.compute( input, densityDeriv );

  dValue[Deriv::dP] = densityDeriv[0];
  dValue[Deriv::dT] = densityDeriv[1];
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_WATERDENSITY_HPP_
