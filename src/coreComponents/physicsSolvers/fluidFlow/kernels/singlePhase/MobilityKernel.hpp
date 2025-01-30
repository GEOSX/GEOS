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
 * @file MobilityKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MOBILITYKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MOBILITYKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** MobilityKernel ********************************/

struct MobilityKernel
{

// Value-only (no derivatives) version
  GEOS_HOST_DEVICE
  inline
  static void
  compute( real64 const & dens,
           real64 const & visc,
           real64 & mob )
  {
    mob = dens / visc;
  }

// derivatives-only  version
  GEOS_HOST_DEVICE
  inline
  static void
  compute_derivative( real64 const & mobility,
                      real64 const & dDensity,
                      real64 const & viscosity,
                      real64 const & dViscosity,
                      real64 & dMobility )
  {
    dMobility = dDensity/viscosity - mobility/viscosity*dViscosity;
  }

  // Computes mobility and derivtives
  template< typename POLICY, integer NUMDOF >
  static void compute_value_and_derivatives( localIndex const size,
                                             arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & density,
                                             arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dDensity,
                                             arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & viscosity,
                                             arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > const & dViscosity,
                                             arrayView1d< real64 > const & mobility,
                                             arrayView2d< real64, constitutive::singlefluid::USD_FLUID > const & dMobility )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      compute( density[a][0], viscosity[a][0], mobility[a] );
      for( int i=0; i<NUMDOF; i++ )
      {
        compute_derivative( mobility[a], dDensity[a][0][i], viscosity[a][0], dViscosity[a][0][i], dMobility[a][i] );
      }

    } );
  }

  // Value-only (no derivatives) version
  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & dens,
                      arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const & visc,
                      arrayView1d< real64 > const & mob )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      compute( dens[a][0],
               visc[a][0],
               mob[a] );
    } );
  }
};

} // namespace singlePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MOBILITYKERNEL_HPP
