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
 * @file NoOpPVTFunction.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NOOPPVTFUNCTION_HPP_
#define GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NOOPPVTFUNCTION_HPP_

#include "PVTFunctionBase.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

class NoOpPVTFunctionUpdate final : public PVTFunctionBaseUpdate
{
public:

  NoOpPVTFunctionUpdate( arrayView1d< real64 const > const & componentMolarWeight )
    : PVTFunctionBaseUpdate( componentMolarWeight )
  {}

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                real64 & value,
                arraySlice1d< real64, USD3 > const & dValue,
                bool useMass ) const
  {
    GEOS_UNUSED_VAR( pressure, temperature,
                     phaseComposition, dPhaseComposition,
                     value, dValue,
                     useMass );
  }

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    PVTFunctionBaseUpdate::move( space, touch );
  }

};


class NoOpPVTFunction : public PVTFunctionBase
{
public:

  NoOpPVTFunction( string const & name,
                   string_array const & inputPara,
                   string_array const & componentNames,
                   array1d< real64 > const & componentMolarWeight,
                   bool const printTable )
    : PVTFunctionBase( name,
                       componentNames,
                       componentMolarWeight )
  {
    GEOS_UNUSED_VAR( inputPara, printTable );
  }

  virtual ~NoOpPVTFunction() override = default;

  static string catalogName() { return "NoOpPVTFunction"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  /**
   * @copydoc PVTFunctionBase::checkTablesParameters( real64 pressure, real64 temperature )
   */
  virtual void checkTablesParameters( real64 GEOS_UNUSED_PARAM( pressure ),
                                      real64 GEOS_UNUSED_PARAM( temperature ) ) const override final
  {}

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::UNKNOWN;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = NoOpPVTFunctionUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_componentMolarWeight );
  };

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_NOOPPVTFUNCTION_HPP_
