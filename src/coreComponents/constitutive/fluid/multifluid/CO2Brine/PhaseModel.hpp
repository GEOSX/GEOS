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
 * @file PhaseModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_PHASEMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_CO2BRINE_PHASEMODEL_HPP_

namespace geos
{

namespace constitutive
{

/**
 * @brief Struct storing the submodels describing the fluid phase behavior.
 *        For now, this is used only in the CO2BrineFluid class
 * @tparam DENS Class describing the density model
 * @tparam VISC Class describing the viscosity model
 * @tparam ENTH Class describing the enthalpy model
 */
template< typename DENS, typename VISC, typename ENTH >
struct PhaseModel
{
  using Density = DENS;
  using Viscosity = VISC;
  using Enthalpy = ENTH;

  /// Enum used in the constructor to make the distinction between submodel params
  enum InputParamOrder : integer
  {
    DENSITY        = 0, ///< position of the density params
    VISCOSITY      = 1, ///< position of the viscosity params
    ENTHALPY       = 2  ///< position of the enthalpy params
  };

  /**
   * @brief Constructor of the struct, in charge of the instantiation of the submodels
   * @param[in] phaseModelName name of the phase model, used only in the instantiation of the submodels
   * @param[in] inputParams input parameters read from files
   * @param[in] componentNames names of the components
   * @param[in] componentMolarWeight molar weights of the components
   */
  PhaseModel( string const & phaseModelName,
              array1d< array1d< string > > const & inputParams,
              string_array const & componentNames,
              array1d< real64 > const & componentMolarWeight,
              bool const printTable )
    : density( phaseModelName + "_" + Density::catalogName(),
               inputParams[InputParamOrder::DENSITY],
               componentNames,
               componentMolarWeight,
               printTable ),
    viscosity( phaseModelName + "_" + Viscosity::catalogName(),
               inputParams[InputParamOrder::VISCOSITY],
               componentNames,
               componentMolarWeight,
               printTable ),
    enthalpy( phaseModelName + "_" + Enthalpy::catalogName(),
              inputParams[InputParamOrder::ENTHALPY],
              componentNames,
              componentMolarWeight,
              printTable )
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
                   Enthalpy const & enth )
      : density( dens.createKernelWrapper() ),
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
};

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_PHASEMODEL_HPP_
