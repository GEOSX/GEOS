/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLE_KERNELLAUNCHSELECTOR_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLE_KERNELLAUNCHSELECTOR_HPP


namespace geos
{

namespace immiscibleMultiphaseKernels
{

/******************************** Kernel launch machinery ********************************/


template< typename T, typename LAMBDA >
void kernelLaunchSelectorPhaseSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "kernelLaunchSelectorPhaseSwitch: type should be integral" );

  switch( value )
  {
    case 2:
    { lambda( std::integral_constant< T, 2 >() ); return; }
    default:
    { GEOS_ERROR( "Unsupported number of phases: " << value ); }
  }

}

} // namespace immiscibleMultiphaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLE_KERNELLAUNCHSELECTOR_HPP
