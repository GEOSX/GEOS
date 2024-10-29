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
 * @file ModelParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_MODELPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_MODELPARAMETERS_HPP_

namespace geos
{

namespace constitutive
{

class MultiFluidBase;

namespace compositional
{

class ComponentProperties;

class ModelParameters
{
public:
  ModelParameters( std::unique_ptr< ModelParameters > parameters = nullptr ): baseParameters( std::move( parameters ) ) {}
  virtual ~ModelParameters() = default;

  void registerParameters( MultiFluidBase * fluid )
  {
    registerParametersImpl( fluid );
    if( baseParameters )
    {
      baseParameters->registerParameters( fluid );
    }
  }

  void postInputInitialization( MultiFluidBase const * fluid, ComponentProperties const & componentProperties )
  {
    postInputInitializationImpl( fluid, componentProperties );
    if( baseParameters )
    {
      baseParameters->postInputInitialization( fluid, componentProperties );
    }
  }

  template< typename PARAMETERS >
  PARAMETERS const * get() const
  {
    PARAMETERS const * parameters = dynamic_cast< PARAMETERS const * >(this);
    if( parameters == nullptr && baseParameters )
    {
      return baseParameters->get< PARAMETERS >();
    }
    return parameters;
  }

private:
  virtual void registerParametersImpl( MultiFluidBase * fluid )
  {
    GEOS_UNUSED_VAR( fluid );
  }

  virtual void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties )
  {
    GEOS_UNUSED_VAR( fluid );
    GEOS_UNUSED_VAR( componentProperties );
  }

private:
  std::unique_ptr< ModelParameters > baseParameters{};
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_MODELPARAMETERS_HPP_
