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
 * @file KernelLaunchSelectors.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_KERNELLAUNCHSELECTORS_HPP
#define GEOS_PHYSICSSOLVERS_KERNELLAUNCHSELECTORS_HPP

namespace geos
{
namespace internal
{

template< typename S, typename T, typename LAMBDA >
void invokePhaseDispatchLambda ( S val, T numPhases, LAMBDA && lambda )
{
  if( numPhases == 1 )
  {
    lambda( val, std::integral_constant< T, 1 >());
    return;
  }
  else if( numPhases == 2 )
  {
    lambda( val, std::integral_constant< T, 2 >());
    return;
  }
  else if( numPhases == 3 )
  {
    lambda( val, std::integral_constant< T, 3 >());
    return;
  }
  else
  {
    GEOS_ERROR( "Unsupported state: " << numPhases );
  }
}

template< typename S, typename T, typename LAMBDA >
void invokeThermalDispatchLambda ( S val, T isThermal, LAMBDA && lambda )
{
  if( isThermal == 1 )
  {
    lambda( val, std::integral_constant< T, 1 >());
    return;
  }
  else if( isThermal == 0 )
  {
    lambda( val, std::integral_constant< T, 0 >());
    return;
  }
  else
  {
    GEOS_ERROR( "Unsupported state: " << isThermal );
  }
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorThermalSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorThermalSwitch: type should be integral" );

  switch( value )
  {
    case 0:
    {
      lambda( std::integral_constant< T, 0 >() );
      return;
    }
    case 1:
    {
      lambda( std::integral_constant< T, 1 >() );
      return;
    }
    default:
    {
      GEOS_ERROR( "Unsupported thermal state: " << value );
    }
  }
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompThermSwitch( T value, bool const isThermal, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: value type should be integral" );


  switch( value )
  {
    case 1:
    {
      invokeThermalDispatchLambda( std::integral_constant< T, 1 >(), isThermal, lambda );  return;
    }
    case 2:
    {
      invokeThermalDispatchLambda( std::integral_constant< T, 2 >(), isThermal, lambda );
      return;
    }
    case 3:
    {
      invokeThermalDispatchLambda( std::integral_constant< T, 3 >(), isThermal, lambda );
      return;
    }
    case 4:
    {
      invokeThermalDispatchLambda( std::integral_constant< T, 4 >(), isThermal, lambda );
      return;
    }
    case 5:
    {
      invokeThermalDispatchLambda( std::integral_constant< T, 5 >(), isThermal, lambda );
      return;
    }
    default:
    {
      GEOS_ERROR( "Unsupported number of components: " << value );
    }
  }
}

template< typename T, typename LAMBDA >
void kernelLaunchSelectorCompPhaseSwitch( T value, T n_phase, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorCompSwitch: value type should be integral" );
  switch( value )
  {
    case 1:
    {
      invokePhaseDispatchLambda( std::integral_constant< T, 1 >(), n_phase, lambda );
      return;
    }
    case 2:
    {
      invokePhaseDispatchLambda( std::integral_constant< T, 2 >(), n_phase, lambda );
      return;
    }
    case 3:
    {
      invokePhaseDispatchLambda( std::integral_constant< T, 3 >(), n_phase, lambda );
      return;
    }
    case 4:
    {
      invokePhaseDispatchLambda( std::integral_constant< T, 4 >(), n_phase, lambda );
      return;
    }
    case 5:
    {
      invokePhaseDispatchLambda( std::integral_constant< T, 5 >(), n_phase, lambda );
      return;
    }
    default:
    { GEOS_ERROR( "Unsupported number of components: " << value ); }
  }
}


} // end namspace internal
} // end namespace geos


#endif // GEOS_PHYSICSSOLVERS_KERNELLAUNCHSELECTORS_HPP
