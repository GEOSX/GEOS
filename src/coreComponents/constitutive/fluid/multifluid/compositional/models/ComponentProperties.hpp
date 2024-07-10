/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ComponentProperties.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPONENTPROPERTIES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPONENTPROPERTIES_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @brief Class holding standard component properties for a compositional fluid model.
 */
class ComponentProperties final
{
public:
  ComponentProperties( string_array const & componentNames,
                       array1d< real64 > const & componentMolarWeight ):
    m_componentNames ( componentNames ),
    m_componentMolarWeight ( componentMolarWeight )
  {}

  ~ComponentProperties() = default;
  ComponentProperties( const ComponentProperties & ) = default;
  const ComponentProperties & operator=( const ComponentProperties & ) = delete;

  /**
   * @brief Get the number of components
   * @return The number of components
   */
  integer getNumberOfComponents() const { return m_componentNames.size(); }

  /**
   * Data accessors
   */
  arrayView1d< real64 > const & getComponentMolarWeight() const { return m_componentMolarWeight; }
  arrayView1d< real64 > const & getComponentCriticalPressure() const { return m_componentCriticalPressure; }
  arrayView1d< real64 > const & getComponentCriticalTemperature() const { return m_componentCriticalTemperature; }
  arrayView1d< real64 > const & getComponentAcentricFactor() const { return m_componentAcentricFactor; }
  arrayView1d< real64 > const & getComponentVolumeShift() const { return m_componentVolumeShift; }

  struct KernelWrapper
  {
    KernelWrapper( arrayView1d< real64 const > const & componentMolarWeight,
                   arrayView1d< real64 const > const & componentCriticalPressure,
                   arrayView1d< real64 const > const & componentCriticalTemperature,
                   arrayView1d< real64 const > const & componentAcentricFactor,
                   arrayView1d< real64 const > const & componentVolumeShift,
                   arrayView2d< real64 const > const & componentBinaryCoeff ):
      m_componentMolarWeight ( componentMolarWeight ),
      m_componentCriticalPressure ( componentCriticalPressure ),
      m_componentCriticalTemperature( componentCriticalTemperature ),
      m_componentAcentricFactor( componentAcentricFactor ),
      m_componentVolumeShift( componentVolumeShift ),
      m_componentBinaryCoeff( componentBinaryCoeff )
    {}

    /**
     * @brief Move the KernelWrapper to the given execution space, optionally touching it.
     * @param space the space to move the KernelWrapper to
     * @param touch whether the KernelWrapper should be touched in the new space or not
     * @note This function exists to enable holding KernelWrapper objects in an ArrayView
     *       and have their contents properly moved between memory spaces.
     */
    void move( LvArray::MemorySpace const space, bool const touch )
    {
      m_componentMolarWeight.move( space, touch );
      m_componentCriticalPressure.move( space, touch );
      m_componentCriticalTemperature.move( space, touch );
      m_componentAcentricFactor.move( space, touch );
      m_componentVolumeShift.move( space, touch );
      m_componentBinaryCoeff.move( space, touch );
    }

    // Standard compositional input
    arrayView1d< real64 const > m_componentMolarWeight;
    arrayView1d< real64 const > m_componentCriticalPressure;
    arrayView1d< real64 const > m_componentCriticalTemperature;
    arrayView1d< real64 const > m_componentAcentricFactor;
    arrayView1d< real64 const > m_componentVolumeShift;
    arrayView2d< real64 const > m_componentBinaryCoeff;
  };

  /**
   * @brief Function to create and return a KernelWrapper
   * @return the KernelWrapper object
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_componentMolarWeight,
                          m_componentCriticalPressure,
                          m_componentCriticalTemperature,
                          m_componentAcentricFactor,
                          m_componentVolumeShift,
                          m_componentBinaryCoeff );
  }

public:
  // Standard compositional input
  string_array const & m_componentNames;
  array1d< real64 > const & m_componentMolarWeight;
  array1d< real64 > m_componentCriticalPressure;
  array1d< real64 > m_componentCriticalTemperature;
  array1d< real64 > m_componentAcentricFactor;
  array1d< real64 > m_componentVolumeShift;
  array2d< real64 > m_componentBinaryCoeff;
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPONENTPROPERTIES_HPP_
