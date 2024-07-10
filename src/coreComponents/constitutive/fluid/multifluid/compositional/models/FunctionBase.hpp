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
 * @file FunctionBase.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_FUNCTIONBASE_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_FUNCTIONBASE_HPP_

#include "ComponentProperties.hpp"
#include "ModelParameters.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

enum class FunctionType : integer
{
  UNKNOWN,
  FLASH,
  DENSITY,
  VISCOSITY,
  ENTHALPY,
  INTERNAL_ENERGY
};

class FunctionBaseUpdate
{
public:
  GEOS_HOST_DEVICE FunctionBaseUpdate(){}
  GEOS_HOST_DEVICE FunctionBaseUpdate( FunctionBaseUpdate const & ){}

  /**
   * @brief Move the KernelWrapper to the given execution space, optionally touching it.
   * @param space the space to move the KernelWrapper to
   * @param touch whether the KernelWrapper should be touched in the new space or not
   * @note This function exists to enable holding KernelWrapper objects in an ArrayView
   *       and have their contents properly moved between memory spaces.
   */
  virtual void move( LvArray::MemorySpace const space, bool const touch )
  {
    GEOS_UNUSED_VAR( space, touch );
  }

protected:
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static void setZero( real64 & val ){ val = 0.0; }
};

class FunctionBase
{

public:

  FunctionBase( string const & name,
                ComponentProperties const & componentProperties ):
    m_functionName( name ),
    m_componentProperties( componentProperties )
  {}

  virtual ~FunctionBase() = default;

  string const & functionName() const { return m_functionName; }

  virtual FunctionType functionType() const = 0;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters )
  {
    return parameters;
  }

protected:
  /// Name of the PVT function
  string m_functionName;

  /// Compositional component properties
  ComponentProperties const & m_componentProperties;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_FUNCTIONBASE_HPP_
