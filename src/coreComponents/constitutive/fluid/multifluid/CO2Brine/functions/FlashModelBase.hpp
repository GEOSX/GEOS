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
 * @file FlashModelBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_FLASHMODELBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_FUNCTIONS_FLASHMODELBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"

namespace geos
{

namespace constitutive
{

namespace PVTProps
{

class FlashModelBaseUpdate
{
public:

  /**
   * @brief Constructor.
   * @param componentMolarWeight component molar weights
   */
  explicit FlashModelBaseUpdate( arrayView1d< real64 const > const & componentMolarWeight )
    :
    m_componentMolarWeight( componentMolarWeight )
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

  /// Array storing the component molar weights
  arrayView1d< real64 const > m_componentMolarWeight;

};


class FlashModelBase
{
public:

  FlashModelBase( string const & name,
                  string_array const & componentNames,
                  array1d< real64 > const & componentMolarWeight ):
    m_modelName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~FlashModelBase() = default;

  virtual string getCatalogName() const = 0;

  /**
   * @brief Check if the input values are in the expected pressure & temperature tables bounds
   * @param pressure input pressure to check
   * @param temperature input temperature to check
   * @throw a SimulationError if one of the input values is out of bound.
   */
  virtual void checkTablesParameters( real64 pressure, real64 temperature ) const = 0;

  string const & flashModelName() const { return m_modelName; }

protected:

  /// Name the solubility model
  string m_modelName;

  /// Array storing the name of the components
  string_array m_componentNames;

  /// Array storing the component molar weights
  array1d< real64 > m_componentMolarWeight;

};

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FLASHMODELBASE_HPP_
