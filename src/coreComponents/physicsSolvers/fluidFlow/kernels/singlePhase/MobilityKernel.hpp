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
           real64 const & dDens_dPres,
           real64 const & visc,
           real64 const & dVisc_dPres,
           real64 & mob,
           real64 & dMob_dPres )
  {
    mob = dens / visc;
    dMob_dPres = dDens_dPres / visc - mob / visc * dVisc_dPres;
  }

// Thermal version
  GEOS_HOST_DEVICE
  inline
  static void
  compute( real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & dDens_dTemp,
           real64 const & visc,
           real64 const & dVisc_dPres,
           real64 const & dVisc_dTemp,
           real64 & mob,
           real64 & dMob_dPres,
           real64 & dMob_dTemp )
  {
    mob = dens / visc;
    dMob_dPres = dDens_dPres / visc - mob / visc * dVisc_dPres;
    dMob_dTemp = dDens_dTemp / visc - mob / visc * dVisc_dTemp;
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
                      arrayView2d< real64 const > const & dDens_dPres,
                      arrayView2d< real64 const > const & visc,
                      arrayView2d< real64 const > const & dVisc_dPres,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      compute( dens[a][0],
               dDens_dPres[a][0],
               visc[a][0],
               dVisc_dPres[a][0],
               mob[a],
               dMob_dPres[a] );
    } );
  }

  // Thermal version
  template< typename POLICY >
  static void launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & dDens_dPres,
                      arrayView2d< real64 const > const & dDens_dTemp,
                      arrayView2d< real64 const > const & visc,
                      arrayView2d< real64 const > const & dVisc_dPres,
                      arrayView2d< real64 const > const & dVisc_dTemp,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres,
                      arrayView1d< real64 > const & dMob_dTemp )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      compute( dens[a][0],
               dDens_dPres[a][0],
               dDens_dTemp[a][0],
               visc[a][0],
               dVisc_dPres[a][0],
               dVisc_dTemp[a][0],
               mob[a],
               dMob_dPres[a],
               dMob_dTemp[a] );
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
