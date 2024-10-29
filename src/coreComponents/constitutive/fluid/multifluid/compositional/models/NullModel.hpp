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
 * @file NullModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLMODEL_HPP_

#include "FunctionBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class NullModelUpdate final : public FunctionBaseUpdate
{
public:
  NullModelUpdate() = default;

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                arraySlice2d< real64 const, USD2 > const & dPhaseComposition,
                real64 & value,
                arraySlice1d< real64, USD3 > const & dValue,
                bool useMass ) const
  {
    GEOS_UNUSED_VAR( componentProperties,
                     pressure, temperature,
                     phaseComposition, dPhaseComposition,
                     value, dValue,
                     useMass );
  }
};

class NullModel : public FunctionBase
{
public:

  NullModel( string const & name,
             ComponentProperties const & componentProperties,
             integer const phaseIndex,
             ModelParameters const & modelParameters ):
    FunctionBase( name, componentProperties )
  {
    GEOS_UNUSED_VAR( phaseIndex, modelParameters );
  }

  NullModel( string const & name,
             ComponentProperties const & componentProperties,
             ModelParameters const & modelParameters ):
    FunctionBase( name, componentProperties )
  {
    GEOS_UNUSED_VAR( modelParameters );
  }

  virtual ~NullModel() override = default;

  static string catalogName() { return "NullPVTModel"; }

  static constexpr FunctionType function(){ return FunctionType::UNKNOWN; }

  virtual FunctionType functionType() const override
  {
    return function();
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = NullModelUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper();
  };

};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NULLMODEL_HPP_
