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
 * @file KernelLaunchSelector.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_KERNELLAUNCHSELECTOR_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_KERNELLAUNCHSELECTOR_HPP

#include "physicsSolvers/KernelLaunchSelectors.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** Kernel launch machinery ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: type should be integral" );

  switch( value )
  {
    case 1:
    { lambda( std::integral_constant< T, 1 >() ); return; }
    case 2:
    { lambda( std::integral_constant< T, 2 >() ); return; }
    case 3:
    { lambda( std::integral_constant< T, 3 >() ); return; }
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return; }
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return; }
    default:
    { GEOS_ERROR( "Unsupported number of components: " << value ); }
  }
}

} // namespace internal

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelectorCompTherm( integer const numComp, bool const isThermal, ARGS && ... args )
{
  geos::internal::kernelLaunchSelectorCompThermSwitch( numComp, isThermal, [&] ( auto NC, auto ISTHERMAL )
  {
    KERNELWRAPPER::template launch< NC(), ISTHERMAL() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector1( integer const numComp, ARGS && ... args )
{
  internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    KERNELWRAPPER::template launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector2( integer const numComp, integer const numPhase, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( numPhase == 2 )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 2 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else if( numPhase == 3 )
  {
    internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
    {
      KERNELWRAPPER::template launch< NC(), 3 >( std::forward< ARGS >( args ) ... );
    } );
  }
  else
  {
    GEOS_ERROR( "Unsupported number of phases: " << numPhase );
  }
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector_NC_NP_THERM( integer const numComp, integer const numPhase, integer const isThermal, ARGS && ... args )
{
  // Ideally this would be inside the dispatch, but it breaks on Summit with GCC 9.1.0 and CUDA 11.0.3.
  if( isThermal )
  {
    if( numPhase == 2 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        KERNELWRAPPER::template launch< NC(), 2, 1 >( std::forward< ARGS >( args ) ... );
      } );
    }
    else if( numPhase == 3 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        KERNELWRAPPER::template launch< NC(), 3, 1 >( std::forward< ARGS >( args ) ... );
      } );
    }
    else
    {
      GEOS_ERROR( "Unsupported number of phases: " << numPhase );
    }
  }
  else
  {
    if( numPhase == 2 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        KERNELWRAPPER::template launch< NC(), 2, 0 >( std::forward< ARGS >( args ) ... );
      } );
    }
    else if( numPhase == 3 )
    {
      internal::kernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
      {
        KERNELWRAPPER::template launch< NC(), 3, 0 >( std::forward< ARGS >( args ) ... );
      } );
    }
    else
    {
      GEOS_ERROR( "Unsupported number of phases: " << numPhase );
    }
  }
}

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_KERNELLAUNCHSELECTOR_HPP
