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
 * @file multiFluidSelector.hpp
 */
#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/fluid/multifluid/blackOil/DeadOilFluid.hpp"
#include "constitutive/fluid/multifluid/blackOil/BlackOilFluid.hpp"
#include "constitutive/fluid/multifluid/CO2Brine/CO2BrineFluid.hpp"

#include "common/TypeDispatch.hpp"
#ifdef GEOS_USE_PVTPackage
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluidPVTPackage.hpp"
#endif
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"

namespace geos
{

namespace constitutive
{

template< typename LAMBDA >
void constitutiveUpdatePassThru( constitutive::MultiFluidBase const & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
#ifdef GEOS_USE_PVTPackage
                               CompositionalMultiphaseFluidPVTPackage,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid,
// Including these in a CUDA build will lead to compiler segfault.
// Need to split compilation units for all the options
#if !defined(GEOS_DEVICE_COMPILE)
                               CO2BrineEzrokhiThermalFluid,
                               CompositionalTwoPhaseLohrenzBrayClarkViscosity,
#endif
                               CompositionalTwoPhaseConstantViscosity
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( constitutive::MultiFluidBase & fluid,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< DeadOilFluid,
                               BlackOilFluid,
#ifdef GEOS_USE_PVTPackage
                               CompositionalMultiphaseFluidPVTPackage,
#endif
                               CO2BrinePhillipsFluid,
                               CO2BrineEzrokhiFluid,
                               CO2BrinePhillipsThermalFluid,
// Including these in a CUDA build will lead to compiler segfault.
// Need to split compilation units for all the options"
#if !defined(GEOS_DEVICE_COMPILE)
                               CO2BrineEzrokhiThermalFluid,
                               CompositionalTwoPhaseLohrenzBrayClarkViscosity,
#endif
                               CompositionalTwoPhaseConstantViscosity
                               >::execute( fluid, std::forward< LAMBDA >( lambda ) );
}

namespace detail
{
/**
 * @brief Structure to execute a lambda on a fluid model depending on the number of components
 * @tparam COMPONENT_LIST Type listing the components to expand over
 */
template< typename COMPONENT_LIST >
struct ComponentSelector;

template< camp::idx_t ... Is >
struct ComponentSelector< camp::idx_seq< Is ... > >
{
  template< typename FluidType, typename LAMBDA >
  static void execute( int numComps, FluidType & fluid, LAMBDA && lambda )
  {
// With gcc8, the fold expression below issues a spurious warning
// warning: suggest parentheses around '&&' within '||' [-Wparentheses]
// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94505
#if (defined(__GNUC__) && (__GNUC__ < 10))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#endif

    bool const supported = ( ((numComps == Is) && (lambda( fluid, std::integral_constant< integer, Is >() ), true)) || ...);

#if (defined(__GNUC__) && (__GNUC__ < 10))
#pragma GCC diagnostic pop
#endif
    GEOS_THROW_IF( !supported,
                   "Unsupported number of components: " << numComps << " for fluid " << FluidType::catalogName(),
                   InputError );
  }
};

}

template< bool THERMAL = false, typename LAMBDA = NoOpFunc >
void constitutiveComponentUpdatePassThru( constitutive::MultiFluidBase & fluidBase,
                                          integer const numComps,
                                          LAMBDA && lambda )
{
  constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
  {
    using FluidType = TYPEOFREF( fluid );
    static_assert( FluidType::min_n_components <= FluidType::max_n_components );
    using Components = typename types::IntegerSequence< FluidType::min_n_components, FluidType::max_n_components >::type;
    if constexpr (!THERMAL || FluidType::isThermalType())
    {
      detail::ComponentSelector< Components >::execute( numComps, fluid, std::forward< LAMBDA >( lambda ));
    }
    else
    {
      GEOS_THROW( "Unsupported thermal call for fluid " << FluidType::catalogName(),
                  InputError );
    }
  } );
}

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUIDSELECTOR_HPP_
