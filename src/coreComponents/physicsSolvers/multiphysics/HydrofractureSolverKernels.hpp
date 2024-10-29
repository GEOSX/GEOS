/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrofractureSolverKernels.hpp
 *
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_

#include "HydrofractureSolverKernels.hpp"

namespace geos
{

namespace hydrofractureSolverKernels
{

struct DeformationUpdateKernel
{

  template< typename POLICY, typename CONTACT_WRAPPER, typename POROUS_WRAPPER >
  static std::tuple< double, double, double, double, double, double >
  launch( localIndex const size,
          CONTACT_WRAPPER const & contactWrapper,
          POROUS_WRAPPER const & porousMaterialWrapper,
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u,
          arrayView2d< real64 const > const & faceNormal,
          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
          ArrayOfArraysView< localIndex const > const & elemsToFaces,
          arrayView1d< real64 const > const & area,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 > const & deltaVolume,
          arrayView1d< real64 > const & aperture,
          arrayView1d< real64 > const & hydraulicAperture
#ifdef GEOS_USE_SEPARATION_COEFFICIENT
          ,
          arrayView1d< real64 const > const & apertureAtFailure,
          arrayView1d< real64 > const & separationCoeff,
          arrayView1d< real64 > const & dSeparationCoeff_dAper,
          arrayView1d< real64 const > const & separationCoeff0
#endif
          )
  {

    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxApertureChange( 0.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxHydraulicApertureChange( 0.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minAperture( 1e10 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxAperture( -1e10 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > minHydraulicAperture( 1e10 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxHydraulicAperture( -1e10 );

    forAll< POLICY >( size,
                      [=] GEOS_HOST_DEVICE ( localIndex const kfe ) mutable
    {
      if( elemsToFaces.sizeOfArray( kfe ) != 2 )
      { return; }

      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const kf1 = elemsToFaces[kfe][1];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      real64 temp[ 3 ] = { 0 };
      for( localIndex a=0; a<numNodesPerFace; ++a )
      {
        LvArray::tensorOps::add< 3 >( temp, u[ faceToNodeMap( kf0, a ) ] );
        LvArray::tensorOps::subtract< 3 >( temp, u[ faceToNodeMap( kf1, a ) ] );
      }

      // TODO this needs a proper contact based strategy for aperture
      real64 const normalJump = -LvArray::tensorOps::AiBi< 3 >( temp, faceNormal[kf0] ) / numNodesPerFace;
      maxApertureChange.max( std::fabs( normalJump - aperture[kfe] ));
      aperture[kfe] = normalJump;
      minAperture.min( aperture[kfe] );
      maxAperture.max( aperture[kfe] );

      real64 dHydraulicAperture_dNormalJump = 0;
      real64 const newHydraulicAperture = contactWrapper.computeHydraulicAperture( aperture[kfe], dHydraulicAperture_dNormalJump );
      maxHydraulicApertureChange.max( std::fabs( newHydraulicAperture - hydraulicAperture[kfe] ));
      real64 const oldHydraulicAperture = hydraulicAperture[kfe];
      hydraulicAperture[kfe] = newHydraulicAperture;
      minHydraulicAperture.min( hydraulicAperture[kfe] );
      maxHydraulicAperture.max( hydraulicAperture[kfe] );

      real64 const jump[3] = { normalJump, 0.0, 0.0 };
      real64 const traction[3] = {0.0, 0.0, 0.0};

      porousMaterialWrapper.updateStateFromPressureApertureJumpAndTraction( kfe, 0, 0.0,
                                                                            oldHydraulicAperture, newHydraulicAperture,
                                                                            dHydraulicAperture_dNormalJump,
                                                                            jump, traction );

#ifdef GEOS_USE_SEPARATION_COEFFICIENT
      real64 const s = aperture[kfe] / apertureAtFailure[kfe];
      if( separationCoeff0[kfe]<1.0 && s>separationCoeff0[kfe] )
      {
        if( s >= 1.0 )
        {
          separationCoeff[kfe] = 1.0;
          dSeparationCoeff_dAper[kfe] = 0.0;
        }
        else
        {
          separationCoeff[kfe] = s;
          dSeparationCoeff_dAper[kfe] = 1.0/apertureAtFailure[kfe];
        }
      }
#endif
      deltaVolume[kfe] = hydraulicAperture[kfe] * area[kfe] - volume[kfe];
    } );

    return std::make_tuple( maxApertureChange.get(), maxHydraulicApertureChange.get(), minAperture.get(), maxAperture.get(), minHydraulicAperture.get(), maxHydraulicAperture.get() );
  }
};

struct FluidMassResidualDerivativeAssemblyKernel
{
  template< typename CONTACT_WRAPPER >
  GEOS_HOST_DEVICE
  inline
  static void
  computeAccumulationDerivative( CONTACT_WRAPPER const & contactWrapper,
                                 localIndex const numNodesPerFace,
                                 arraySlice1d< localIndex const > const elemsToFaces,
                                 ArrayOfArraysView< localIndex const > const faceToNodeMap,
                                 arrayView1d< globalIndex const > const dispDofNumber,
                                 real64 const (&Nbar)[ 3 ],
                                 real64 const & area,
                                 real64 const & aperture,
                                 real64 const & dens,
                                 globalIndex (& nodeDOF)[8 * 3],
                                 arraySlice1d< real64 > const dRdU )
  {
    real64 dHydraulicAperture_dNormalJump = 0;
    real64 const hydraulicAperture = contactWrapper.computeHydraulicAperture( aperture, dHydraulicAperture_dNormalJump );
    GEOS_UNUSED_VAR( hydraulicAperture );

    constexpr integer kfSign[2] = { -1, 1 };
    for( localIndex kf = 0; kf < 2; ++kf )
    {
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        for( int i = 0; i < 3; ++i )
        {
          nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kf], a )] + i;

          real64 const dNormalJump_dDisplacement = kfSign[kf] * Nbar[i] / numNodesPerFace;
          real64 const dHydraulicAperture_dDisplacement = dHydraulicAperture_dNormalJump * dNormalJump_dDisplacement;
          real64 const dVolume_dDisplacement = area * dHydraulicAperture_dDisplacement;

          dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dens * dVolume_dDisplacement;
        }
      }
    }
  }

  GEOS_HOST_DEVICE
  inline
  static void
  computeFluxDerivative( localIndex const kfe2,
                         localIndex const numNodesPerFace,
                         arraySlice1d< localIndex const > const & columns,
                         arraySlice1d< real64 const > const & values,
                         ArrayOfArraysView< localIndex const > const elemsToFaces,
                         ArrayOfArraysView< localIndex const > const faceToNodeMap,
                         arrayView1d< globalIndex const > const dispDofNumber,
                         real64 const (&Nbar)[ 3 ],
                         globalIndex (& nodeDOF)[8 * 3],
                         arraySlice1d< real64 > const dRdU )
  {
    constexpr integer kfSign[2] = { -1, 1 };

    real64 const dR_dNormalJump = values[kfe2];
    localIndex const ei2 = columns[kfe2];

    for( localIndex kf = 0; kf < 2; ++kf )
    {
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[ei2][kf], a )] + i;

          real64 const dNormalJump_dDisplacement = kfSign[kf] * Nbar[i] / numNodesPerFace;

          dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dR_dNormalJump * dNormalJump_dDisplacement;
        }
      }
    }
  }

  template< typename POLICY, typename CONTACT_WRAPPER >
  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          CONTACT_WRAPPER const & contactWrapper,
          integer const useQuasiNewton,
          ArrayOfArraysView< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const faceToNodeMap,
          arrayView2d< real64 const > const faceNormal,
          arrayView1d< real64 const > const area,
          arrayView1d< real64 const > const aperture,
          arrayView1d< globalIndex const > const presDofNumber,
          arrayView1d< globalIndex const > const dispDofNumber,
          arrayView2d< real64 const > const dens,
          CRSMatrixView< real64 const, localIndex const > const dFluxResidual_dNormalJump,
          CRSMatrixView< real64, globalIndex const > const & localMatrix )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex ei )
    {
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[ei][0] );

      real64 Nbar[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[elemsToFaces[ei][0]] );
      LvArray::tensorOps::subtract< 3 >( Nbar, faceNormal[elemsToFaces[ei][1]] );
      LvArray::tensorOps::normalize< 3 >( Nbar );

      globalIndex const rowNumber = presDofNumber[ei] - rankOffset;
      globalIndex nodeDOF[8 * 3];
      stackArray1d< real64, 24 > dRdU( 2 * numNodesPerFace * 3 );
//
      computeAccumulationDerivative( contactWrapper,
                                     numNodesPerFace,
                                     elemsToFaces[ei],
                                     faceToNodeMap,
                                     dispDofNumber,
                                     Nbar,
                                     area[ei],
                                     aperture[ei],
                                     dens[ei][0],
                                     nodeDOF,
                                     dRdU );

      if( rowNumber >= 0  && rowNumber < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                          nodeDOF,
                                                                          dRdU.data(),
                                                                          2 * numNodesPerFace * 3 );
      }
//
      if( useQuasiNewton == 0 ) // when Quasi Newton is not enabled - add flux derivatives
      {
        localIndex const numColumns = dFluxResidual_dNormalJump.numNonZeros( ei );
        arraySlice1d< localIndex const > const & columns = dFluxResidual_dNormalJump.getColumns( ei );
        arraySlice1d< real64 const > const & values = dFluxResidual_dNormalJump.getEntries( ei );

        for( localIndex kfe2 = 0; kfe2 < numColumns; ++kfe2 )
        {
          computeFluxDerivative( kfe2,
                                 numNodesPerFace,
                                 columns,
                                 values,
                                 elemsToFaces,
                                 faceToNodeMap,
                                 dispDofNumber,
                                 Nbar,
                                 nodeDOF,
                                 dRdU );

          if( rowNumber >= 0 && rowNumber < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                              nodeDOF,
                                                                              dRdU.data(),
                                                                              2 * numNodesPerFace * 3 );
          }
        }
      }
    } );
  }
};

} /* namespace hydrofractureSolverKernels */

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_
