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
 * @file TaperKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_TAPERKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_TAPERKERNEL_HPP_

#include "WaveSolverUtils.hpp"

namespace geos
{

using wsCoordType = WaveSolverUtils::wsCoordType;

struct TaperKernel
{

  /**
   * @brief Compute coefficients for the taper layers
   * @tparam EXEC_POLICY the execution policy
   * @param[in] size the number of nodes
   * @param[in] sizeT Taper thickness
   * @param[in] dt time-step
   * @param[in] vMin Min wavespeed (P-wavespeed for acoustic, S-wavespeed for elastic)
   * @param[in] r desired reflectivity of the Taper
   * @param[out] taperCoeff array which contains the taper coefficient on each node (which will be equal to 1 when we are outside of the
   * taper layers)
   */
  template< typename EXEC_POLICY >
  static void
  computeTaperCoeff( localIndex const size,
                     arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                     real32 const sizeT,
                     real32 const dt,
                     real32 const vMax,
                     real32 const r,
                     arrayView1d< real32 > const taperCoeff )
  {

    ///Seek the global maximum and minimum of the domain
    RAJA::ReduceMin< parallelDeviceReduce, real32 > xMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > yMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > zMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > xMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > yMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > zMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > xMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > yMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > zMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > xMaxInterior( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > yMaxInterior( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > zMaxInterior( -LvArray::NumericLimits< real32 >::max );

    real32 xGlobalMin[3]{};
    real32 xGlobalMax[3]{};

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      xMinGlobal.min( nodeCoords[a][0] );
      yMinGlobal.min( nodeCoords[a][1] );
      zMinGlobal.min( nodeCoords[a][2] );
      xMaxGlobal.max( nodeCoords[a][0] );
      yMaxGlobal.max( nodeCoords[a][1] );
      zMaxGlobal.max( nodeCoords[a][2] );
    } );

    xGlobalMin[0] = xMinGlobal.get();
    xGlobalMin[1] = yMinGlobal.get();
    xGlobalMin[2] = zMinGlobal.get();
    xGlobalMax[0] = xMaxGlobal.get();
    xGlobalMax[1] = yMaxGlobal.get();
    xGlobalMax[2] = zMaxGlobal.get();


    for( integer i=0; i<3; ++i )
    {
      xGlobalMin[i] = MpiWrapper::min( xGlobalMin[i] );
      xGlobalMax[i] = MpiWrapper::max( xGlobalMax[i] );
    }

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      real32 dist=0;

      real32 distXmin = LvArray::math::abs( nodeCoords[a][0]-xGlobalMin[0] );
      real32 distXmax = LvArray::math::abs( nodeCoords[a][0]-xGlobalMax[0] );
      real32 distYmin = LvArray::math::abs( nodeCoords[a][1]-xGlobalMin[1] );
      real32 distYmax = LvArray::math::abs( nodeCoords[a][1]-xGlobalMax[1] );
      real32 distZmax = LvArray::math::abs( nodeCoords[a][2]-xGlobalMax[2] );

      dist=LvArray::math::min( distXmin, (LvArray::math::min( distXmax, LvArray::math::min( distYmin, (LvArray::math::min( distYmax, distZmax ))))));

      taperCoeff[a] = dist;


    } );

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {

      real32 dist = taperCoeff[a];

      if( dist<sizeT )
      {
        taperCoeff[a] = LvArray::math::exp((((3*vMax)/(2*sizeT))*log( r )*pow((sizeT-dist)/sizeT, 2 ))*dt );
      }
      else
      {
        taperCoeff[a] = 1.0;
      }

    } );

  }


  /**
   * @brief Multiply an array with the taper coefficients
   * @tparam EXEC_POLICY the execution policy
   * @param[in] size the number of nodes
   * @param[in] taperCoeff array which contains the taper coefficient on each node (which will be equal to 1 when we are outside of the
   * taper layers)
   * @param[inout] vector array which is multiplied by the taper array
   */
  template< typename EXEC_POLICY >
  static void
  multiplyByTaperCoeff( localIndex const size,
                        arrayView1d< real32 const > const taperCoeff,
                        arrayView1d< real32 > const vector )
  {
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      vector[a] *= taperCoeff[a];
    } );

  }


};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_TAPERERNEL_HPP_
