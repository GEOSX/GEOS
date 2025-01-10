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
 * @file MobilityKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MOBILITYKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_MOBILITYKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** MobilityKernel ********************************/

struct MobilityKernel
{
  // Isothermal version
  GEOS_HOST_DEVICE
  inline
  static void
  compute( real64 const & dens,
           real64 const & dDens_dP,
           real64 const & visc,
           real64 const & dVisc_dP,  //tjb
           real64 & mob,
           real64 & dMob_dPres )
  {
    mob = dens / visc;
    dMob_dPres = dDens_dP / visc - mob / visc * dVisc_dP;  // tjb keep
  }

// Thermal version
  GEOS_HOST_DEVICE
  inline
  static void
  old_compute( real64 const & dens,
               real64 const & dDens_dP, // tjb
               real64 const & dDens_dT, // tjb
               real64 const & dDens_dTemp,
               real64 const & visc,
               real64 const & dVisc_dP, // tjb
               real64 const & dVisc_dT, // tjb
               real64 const & dVisc_dTemp,
               real64 & mob,
               real64 & dMob_dPres,
               real64 & dMob_dTemp )
  {
    mob = dens / visc;
    dMob_dPres = dDens_dP / visc - mob / visc * dVisc_dP;
    dMob_dTemp = dDens_dTemp / visc - mob / visc * dVisc_dTemp;
    dMob_dTemp = dDens_dT / visc - mob / visc * dVisc_dT;  // tjb keep
    assert( fabs( dDens_dT - dDens_dTemp ) < FLT_EPSILON );
    assert( fabs( dVisc_dT - dVisc_dTemp ) < FLT_EPSILON );
  }

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

  // Isothermal version
  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView3d< real64 const > const & dDens,
                      arrayView2d< real64 const > const & visc,
                      arrayView3d< real64 const > const & dVisc,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres )
  {
    using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 0 >;
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      compute( dens[a][0],
               dDens[a][0][DerivOffset::dP],
               visc[a][0],
               dVisc[a][0][DerivOffset::dP],   // tjb use deriv::dp
               mob[a],
               dMob_dPres[a] );
    } );
  }

  // Generic version
  template< typename POLICY, integer NUMDOF >
  static void compute_value_and_derivatives( localIndex const size,
                                             arrayView2d< real64 const > const & density,
                                             arrayView3d< real64 const > const & dDensity,
                                             arrayView2d< real64 const > const & viscosity,
                                             arrayView3d< real64 const > const & dViscosity,
                                             arrayView1d< real64 > const & mobility,
                                             arrayView2d< real64 > const & dMobility )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      mobility[a] = density[a][0] / viscosity[a][0];
      for( int i=0; i<NUMDOF; i++ )
      {
        dMobility[a][i] = dDensity[a][0][i]/viscosity[a][0] - mobility[a]/viscosity[a][0]*dViscosity[a][0][i];
      }

    } );
  }

  // Thermal version
  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView3d< real64 const > const & dDens, // tjb
                      arrayView2d< real64 const > const & dDens_dTemp,
                      arrayView2d< real64 const > const & visc,
                      arrayView3d< real64 const > const & dVisc, // tjb
                      arrayView2d< real64 const > const & dVisc_dTemp,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres,
                      arrayView1d< real64 > const & dMob_dTemp )
  {
    using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      old_compute( dens[a][0],
                   dDens[a][0][DerivOffset::dP], // tjb use deriv::dp
                   dDens[a][0][DerivOffset::dT], // tjb use deriv::dt
                   dDens_dTemp[a][0],
                   visc[a][0],
                   dVisc[a][0][DerivOffset::dP], // tjb use deriv::dp
                   dVisc[a][0][DerivOffset::dT], // tjb use deriv::dt
                   dVisc_dTemp[a][0],
                   mob[a],
                   dMob_dPres[a],
                   dMob_dTemp[a] );
      /*
         compute( dens[a][0],
         dDens[a][0],
         visc[a][0],
         dVisc[a][0],
         mob[a],
         dMob_dPres[a],
         dMob_dTemp[a] );
       */
    } );
  }

// Value-only (no derivatives) version
  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & visc,
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
