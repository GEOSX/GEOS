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
 * @file PVTFunctionBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_PVTFUNCTIONBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_PVTFUNCTIONBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

enum class PVTFunctionType { UNKNOWN, DENSITY, VISCOSITY, ENTHALPY, INTERNAL_ENERGY };

class PVTFunctionBaseUpdate
{
public:

  explicit PVTFunctionBaseUpdate( arrayView1d< real64 const > const & componentMolarWeight )
    : m_componentMolarWeight( componentMolarWeight )
  {}

  /**
   * @brief Move the KernelWrapper to the given execution space, optionally touching it.
   * @param space the space to move the KernelWrapper to
   * @param touch whether the KernelWrapper should be touched in the new space or not
   * @note This function exists to enable holding KernelWrapper objects in an ArrayView
   *       and have their contents properly moved between memory spaces.
   */
  virtual void move( LvArray::MemorySpace const space, bool const touch )
  {
    m_componentMolarWeight.move( space, touch );
  }

protected:

  template< int USD >
  GEOS_HOST_DEVICE
  real64 computePhaseMolarWeight( arraySlice1d< real64 const, USD > const & phaseComposition ) const
  {
    integer const numComp = phaseComposition.size();
    real64 MT = 0.0;
    for( integer i = 0; i < numComp; i++ )
    {
      MT += phaseComposition[i] * m_componentMolarWeight[i];
    }
    return MT;
  }

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void divideByPhaseMolarWeight( arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                 arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                                 real64 & value, arraySlice1d< real64, USD3 > const & dValue ) const
  {
    integer const numComp = phaseComposition.size();
    integer const numDerivs = dValue.size();

    real64 const MT = computePhaseMolarWeight( phaseComposition );

    value /= MT;

    for( int der = 0; der < numDerivs; der++ )
    {
      real64 dMT = 0.0;
      for( int ic = 0; ic < numComp; ic++ )
      {
        dMT += dPhaseComposition[ic][der] * m_componentMolarWeight[ic];
      }
      dValue[der] = ( dValue[der] - value * dMT ) / MT; // value is already divided by MT
    }
  }

  /// Array storing the component molar weights
  arrayView1d< real64 const > m_componentMolarWeight;

};

class PVTFunctionBase
{

public:

  PVTFunctionBase( string const & name,
                   array1d< string > const & componentNames,
                   array1d< real64 > const & componentMolarWeight )
    :
    m_functionName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~PVTFunctionBase() = default;

  using CatalogInterface = dataRepository::CatalogInterface< PVTFunctionBase,
                                                             string const &,
                                                             array1d< string > const &,
                                                             array1d< string > const &,
                                                             array1d< real64 > const &,
                                                             bool const >;
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string getCatalogName() const = 0;

  string const & functionName() const { return m_functionName; }

  virtual PVTFunctionType functionType() const = 0;

  /**
   * @brief Check if the input values are in the expected pressure & temperature tables bounds
   * @param pressure input pressure to check
   * @param temperature input temperature to check (in C)
   * @throw a SimulationError if one of the input values is out of bound.
   */
  virtual void checkTablesParameters( real64 pressure, real64 temperature ) const = 0;

protected:

  /// Name of the PVT function
  string m_functionName;

  /// Array storing the name of the components
  array1d< string > m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONBASE_HPP_
