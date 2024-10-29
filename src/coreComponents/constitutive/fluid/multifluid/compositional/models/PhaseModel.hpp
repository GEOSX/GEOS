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
 * @file PhaseModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHASEMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHASEMODEL_HPP_

#include "NullModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ModelParameters;

/**
 * @brief Struct storing the submodels describing the fluid phase behavior.
 * @tparam DENSITY Class describing the density model
 * @tparam VISCOSITY Class describing the viscosity model
 * @tparam ENTHALPY Class describing the enthalpy model
 */
template< typename DENSITY, typename VISCOSITY, typename ENTHALPY >
struct PhaseModel
{
  using Density = DENSITY;
  using Viscosity = VISCOSITY;
  using Enthalpy = ENTHALPY;

  /**
   * @brief Constructor of the struct, in charge of the instantiation of the submodels
   * @param[in] phaseModelName name of the phase model, used only in the instantiation of the submodels
   * @param[in] componentNames names of the components
   * @param[in] componentMolarWeight molar weights of the components
   * @param[in] componentProperties EOS parameters for components
   */
  PhaseModel( string const & phaseModelName,
              ComponentProperties const & componentProperties,
              integer const phaseIndex,
              ModelParameters const & modelParameters ):
    density( phaseModelName + "_" + Density::catalogName(),
             componentProperties,
             phaseIndex,
             modelParameters ),
    viscosity( phaseModelName + "_" + Viscosity::catalogName(),
               componentProperties,
               phaseIndex,
               modelParameters ),
    enthalpy( phaseModelName + "_" + Enthalpy::catalogName(),
              componentProperties,
              phaseIndex,
              modelParameters )
  {}

  /// The phase density model
  Density density;

  /// The phase viscosity model
  Viscosity viscosity;

  /// The phase enthalpy model (can be NoOp for non-thermal models)
  Enthalpy enthalpy;

  /**
   * @brief Struct storing the submodels wrappers used for in-kernel fluid updates
   */
  struct KernelWrapper
  {

    /**
     * @brief Constructor for the kernel wrapper
     * @param[in] dens the density model
     * @param[in] visc the viscosity model
     * @param[in] enth the enthalpy model
     */
    KernelWrapper( Density const & dens,
                   Viscosity const & visc,
                   Enthalpy const & enth ):
      density( dens.createKernelWrapper() ),
      viscosity( visc.createKernelWrapper() ),
      enthalpy( enth.createKernelWrapper() )
    {}

    /**
     * @brief Moves the submodel data to the GPU if the phaseModel is stored in ArrayView
     * @param[in] space the target memory space
     * @param[in] touch the flag to register touch or not
     */
    void move( LvArray::MemorySpace const space, bool const touch )
    {
      density.move( space, touch );
      viscosity.move( space, touch );
      enthalpy.move( space, touch );
    }

    /// Kernel wrapper for density updates
    typename Density::KernelWrapper density;

    /// Kernel wrapper for viscosity updates
    typename Viscosity::KernelWrapper viscosity;

    /// Kernel wrapper for enthalpy updates
    typename Enthalpy::KernelWrapper enthalpy;

  };

  /**
   * @brief Function to create and return a KernelWrapper
   * @return the KernelWrapper object
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( density,
                          viscosity,
                          enthalpy );
  }

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters )
  {
    std::unique_ptr< ModelParameters > phaseParameters = std::move( parameters );
    phaseParameters = Density::createParameters( std::move( phaseParameters ) );
    phaseParameters = Viscosity::createParameters( std::move( phaseParameters ) );
    phaseParameters = Enthalpy::createParameters( std::move( phaseParameters ) );
    return phaseParameters;
  }

};

// A no-op phase model
using NullPhaseModel = PhaseModel< NullModel, NullModel, NullModel >;

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHASEMODEL_HPP_
