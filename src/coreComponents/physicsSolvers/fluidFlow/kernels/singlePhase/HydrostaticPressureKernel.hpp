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
 * @file HydrostaticPressureKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_HYDROSTATICPRESSUREKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_HYDROSTATICPRESSUREKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureKernel
{

  template< typename FLUID_WRAPPER >
  static bool
  computeHydrostaticPressure( integer const maxNumEquilIterations,
                              real64 const & equilTolerance,
                              real64 const (&gravVector)[ 3 ],
                              FLUID_WRAPPER fluidWrapper,
                              real64 const & refElevation,
                              real64 const & refPres,
                              real64 const & refDens,
                              real64 const & newElevation,
                              real64 & newPres,
                              real64 & newDens )
  {
    // Step 1: guess the pressure with the refDensity

    real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
    real64 pres0 = refPres - refDens * gravCoef;
    real64 pres1 = 0.0;

    // Step 2: compute the mass density at this elevation using the guess, and update pressure

    real64 dens = 0.0;
    real64 visc = 0.0;
    constitutive::SingleFluidBaseUpdate::computeValues( fluidWrapper,
                                                        pres0,
                                                        dens,
                                                        visc );
    pres1 = refPres - 0.5 * ( refDens + dens ) * gravCoef;

    // Step 3: fixed-point iteration until convergence

    bool equilHasConverged = false;
    for( localIndex eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
    {

      // check convergence
      equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
      pres0 = pres1;

      // if converged, move on
      if( equilHasConverged )
      {
        break;
      }

      // compute the density at this elevation using the previous pressure, and compute the new pressure
      constitutive::SingleFluidBaseUpdate::computeValues( fluidWrapper,
                                                          pres0,
                                                          dens,
                                                          visc );
      pres1 = refPres - 0.5 * ( refDens + dens ) * gravCoef;
    }

    // Step 4: save the hydrostatic pressure and the corresponding density

    newPres = pres1;
    newDens = dens;

    return equilHasConverged;
  }


  template< typename FLUID_WRAPPER >
  static bool
  launch( localIndex const size,
          integer const maxNumEquilIterations,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & minElevation,
          real64 const & elevationIncrement,
          real64 const & datumElevation,
          real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues,
          arrayView1d< real64 > pressureValues )
  {
    bool hasConverged = true;

    // Step 1: compute the mass density at the datum elevation

    real64 datumDens = 0.0;
    real64 datumVisc = 0.0;

    constitutive::SingleFluidBaseUpdate::computeValues( fluidWrapper,
                                                        datumPres,
                                                        datumDens,
                                                        datumVisc );

    // Step 2: find the closest elevation to datumElevation

    forAll< parallelHostPolicy >( size, [=] ( localIndex const i )
    {
      real64 const elevation = minElevation + i * elevationIncrement;
      elevationValues[0][i] = elevation;
    } );
    integer const iRef = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                 elevationValues[0].size(),
                                                                 datumElevation );


    // Step 3: compute the mass density and pressure at the reference elevation

    array1d< real64 > dens( pressureValues.size() );

    bool const hasConvergedRef =
      computeHydrostaticPressure( maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  datumElevation,
                                  datumPres,
                                  datumDens,
                                  elevationValues[0][iRef],
                                  pressureValues[iRef],
                                  dens[iRef] );
    if( !hasConvergedRef )
    {
      return false;
    }


    // Step 4: for each elevation above the reference elevation, compute the pressure

    localIndex const numEntriesAboveRef = size - iRef - 1;
    forAll< serialPolicy >( numEntriesAboveRef, [=, &hasConverged] ( localIndex const i )
    {
      bool const hasConvergedAboveRef =
        computeHydrostaticPressure( maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    elevationValues[0][iRef+i],
                                    pressureValues[iRef+i],
                                    dens[iRef+i],
                                    elevationValues[0][iRef+i+1],
                                    pressureValues[iRef+i+1],
                                    dens[iRef+i+1] );
      if( !hasConvergedAboveRef )
      {
        hasConverged = false;
      }


    } );

    // Step 5: for each elevation below the reference elevation, compute the pressure

    localIndex const numEntriesBelowRef = iRef;
    forAll< serialPolicy >( numEntriesBelowRef, [=, &hasConverged] ( localIndex const i )
    {
      bool const hasConvergedBelowRef =
        computeHydrostaticPressure( maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    elevationValues[0][iRef-i],
                                    pressureValues[iRef-i],
                                    dens[iRef-i],
                                    elevationValues[0][iRef-i-1],
                                    pressureValues[iRef-i-1],
                                    dens[iRef-i-1] );
      if( !hasConvergedBelowRef )
      {
        hasConverged = false;
      }
    } );

    return hasConverged;
  }
};

} // namespace singlePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_HYDROSTATICPRESSUREKERNEL_HPP
