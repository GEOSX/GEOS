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
 * @file WaveSolverUtils.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_

#include "mesh/utilities/ComputationalGeometry.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

struct WaveSolverUtils
{
  static constexpr real64 epsilonLoc = 1e-8;

  using EXEC_POLICY = parallelDevicePolicy< >;
  using wsCoordType = real32;

  enum class DASType : integer
  {
    none,               ///< deactivate DAS computation
    dipole,             ///< use dipole formulation for DAS
    strainIntegration,  ///< use strain integration for DAS
  };

  enum class AttenuationType : integer
  {
    none,               ///< deactivate attenuation (default)
    sls,                ///< istandard-linear-solid description [Fichtner 2014]
  };


  GEOS_HOST_DEVICE
  static real32 evaluateRicker( real64 const time_n, real32 const f0, real32 const t0, localIndex const order )
  {
    real32 const delay = t0 > 0 ? t0 : 1 / f0;
    real32 pulse = 0.0;
    real32 const alpha = -pow( f0 * M_PI, 2 );
    real32 const time_d = time_n - delay;
    real32 const gaussian = exp( alpha * pow( time_d, 2 ));
    localIndex const sgn = pow( -1, order + 1 );

    switch( order )
    {
      case 0:
        pulse = sgn * gaussian;
        break;
      case 1:
        pulse = sgn * (2 * alpha * time_d) * gaussian;
        break;
      case 2:
        pulse = sgn * (2 * alpha + 4 * pow( alpha, 2 ) * pow( time_d, 2 )) * gaussian;
        break;
      case 3:
        pulse = sgn * (12 * pow( alpha, 2 ) * time_d + 8 * pow( alpha, 3 ) * pow( time_d, 3 )) * gaussian;
        break;
      case 4:
        pulse = sgn * (12 * pow( alpha, 2 ) + 48 * pow( alpha, 3 ) * pow( time_d, 2 ) + 16 * pow( alpha, 4 ) * pow( time_d, 4 )) * gaussian;
        break;
      default:
        GEOS_ERROR( "This option is not supported yet, rickerOrder must be in range {0:4}" );
    }

    return pulse;
  }

  /**
   * @brief Initialize (clear) the trace file.
   * @param[in] prefix Prefix of the output file
   * @param[in] name Name of the solver on which you write the seismo trace
   * @param[in] outputSeismoTrace Boolean equals to 1 if you want to output the seismotrace on a txt file 0 either
   * @param[in] nReceivers Number of receivers
   * @param[in] receiverIsLocal Array to check if the receiver is local to the MPI partition
   */
  static void initTrace( char const * prefix,
                         string const & name,
                         bool const outputSeismoTrace,
                         localIndex const nReceivers,
                         arrayView1d< localIndex const > const receiverIsLocal )
  {
    if( !outputSeismoTrace ) return;

    string const outputDir = OutputBase::getOutputDirectory();
    RAJA::ReduceSum< ReducePolicy< serialPolicy >, localIndex > count( 0 );

    forAll< serialPolicy >( nReceivers, [=] ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        count += 1;
        string const fn = joinPath( outputDir, GEOS_FMT( "{}_{}_{:03}.txt", prefix, name, ircv ) );
        std::ofstream f( fn, std::ios::out | std::ios::trunc );
      }
    } );

    localIndex const total = MpiWrapper::sum( count.get() );
    GEOS_ERROR_IF( nReceivers != total, GEOS_FMT( ": Invalid distribution of receivers: nReceivers={} != MPI::sum={}.", nReceivers, total ) );
  }

  /**
   * @brief Convenient helper for 3D vectors calling 3 times the scalar version with only the sampled variable argument changed.
   * @param[in] prefix Prefix of the output file
   * @param[in] name Name of the solver on which you write the seismo trace
   * @param[in] outputSeismoTrace Boolean equals to 1 if you want to output the seismotrace on a txt file 0 either
   * @param[in] nReceivers Number of receivers
   * @param[in] receiverIsLocal Array to check if the receiver is local to the MPI partition
   * @param[in] nsamplesSeismoTrace Number of samples per seismo trace
   * @param[out] varAtReceiversx Array containing the variable (x-direction) computed at the receivers
   * @param[out] varAtReceiversy Array containing the variable (y-direction) computed at the receivers
   * @param[out] varAtReceiversz Array containing the variable (z-direction) computed at the receivers
   */
  static void writeSeismoTraceVector( char const * prefix,
                                      string const & name,
                                      bool const outputSeismoTrace,
                                      localIndex const nReceivers,
                                      arrayView1d< localIndex const > const receiverIsLocal,
                                      localIndex const nsamplesSeismoTrace,
                                      arrayView2d< real32 const > const varAtReceiversx,
                                      arrayView2d< real32 const > const varAtReceiversy,
                                      arrayView2d< real32 const > const varAtReceiversz )
  {
    writeSeismoTrace( prefix, name, outputSeismoTrace, nReceivers, receiverIsLocal, nsamplesSeismoTrace, varAtReceiversx );
    writeSeismoTrace( prefix, name, outputSeismoTrace, nReceivers, receiverIsLocal, nsamplesSeismoTrace, varAtReceiversy );
    writeSeismoTrace( prefix, name, outputSeismoTrace, nReceivers, receiverIsLocal, nsamplesSeismoTrace, varAtReceiversz );
  }

  /**
   * @brief Write the seismo traces to a file.
   * @param[in] prefix Prefix of the output file
   * @param[in] name Name of the solver on which you write the seismo trace
   * @param[in] outputSeismoTrace Boolean equals to 1 if you want to output the seismotrace on a txt file 0 either
   * @param[in] nReceivers Number of receivers
   * @param[in] receiverIsLocal Array to check if the receiver is local to the MPI partition
   * @param[in] nsamplesSeismoTrace Number of samples per seismo trace
   * @param[in] varAtReceivers Array containing the the variable computed at the receivers
   */
  static void writeSeismoTrace( char const * prefix,
                                string const & name,
                                bool const outputSeismoTrace,
                                localIndex const nReceivers,
                                arrayView1d< localIndex const > const receiverIsLocal,
                                localIndex const nsamplesSeismoTrace,
                                arrayView2d< real32 const > const varAtReceivers )
  {
    if( !outputSeismoTrace ) return;

    string const outputDir = OutputBase::getOutputDirectory();
    forAll< serialPolicy >( nReceivers, [=] ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        string const fn = joinPath( outputDir, GEOS_FMT( "{}_{}_{:03}.txt", prefix, name, ircv ) );
        std::ofstream f( fn, std::ios::app );
        if( f )
        {
          GEOS_LOG_RANK( GEOS_FMT( "Append to seismo trace file {}", fn ) );
          for( localIndex iSample = 0; iSample < nsamplesSeismoTrace; ++iSample )
          {
            // index - time - value
            f << iSample << " " << varAtReceivers[iSample][nReceivers] << " " << varAtReceivers[iSample][ircv] << std::endl;
          }
          f.close();
        }
        else
        {
          GEOS_WARNING( GEOS_FMT( "Failed to open output file {}", fn ) );
        }
      }
    } );
  }

  /**
   * @brief Compute the seismo traces.
   * @param[in] time_n Current time iteration
   * @param[in] dt time-step
   * @param[in] timeSeismo time when the seismo is computed
   * @param[in] iSeismo i-th seismo trace
   * @param[in] receiverNodeIds indices of the nodes of the element where the receiver is located
   * @param[in] receiverConstants constant part of the receiver term
   * @param[in] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[in] var_np1 Array containing the variable at time n+1
   * @param[in] var_n Array containing the variable at time n
   * @param[out] varAtReceivers Array containing the the variable computed at the receivers
   * @param[in] coeffs Coefficients array for receivers
   * @param[in] add Boolean to say if you want to add the value of interpolation to the same receiver coefficient or not
   */
  static void computeSeismoTrace( real64 const time_n,
                                  real64 const dt,
                                  real64 const timeSeismo,
                                  localIndex const iSeismo,
                                  arrayView2d< localIndex const > const receiverNodeIds,
                                  arrayView2d< real64 const > const receiverConstants,
                                  arrayView1d< localIndex const > const receiverIsLocal,
                                  arrayView1d< real32 const > const var_np1,
                                  arrayView1d< real32 const > const var_n,
                                  arrayView2d< real32 > varAtReceivers,
                                  arrayView1d< real32 > coeffs = {},
                                  bool add = false )
  {
    real64 const time_np1 = time_n + dt;

    real32 const a1 = LvArray::math::abs( dt ) < epsilonLoc ? 1.0 : (time_np1 - timeSeismo) / dt;
    real32 const a2 = 1.0 - a1;

    localIndex const nReceivers = receiverConstants.size( 0 );
    forAll< EXEC_POLICY >( nReceivers, [=] GEOS_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] > 0 )
      {
        real32 vtmp_np1 = 0.0, vtmp_n = 0.0;
        for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
        {
          if( receiverNodeIds( ircv, inode ) >= 0 )
          {
            vtmp_np1 += var_np1[receiverNodeIds( ircv, inode )] * receiverConstants( ircv, inode );
            vtmp_n += var_n[receiverNodeIds( ircv, inode )] * receiverConstants( ircv, inode );
          }
        }
        // linear interpolation between the pressure value at time_n and time_{n+1}
        real32 receiverCoeff = coeffs.size( 0 ) == 0 ? 1.0 : coeffs( ircv );
        if( add )
        {
          varAtReceivers( iSeismo, ircv ) += receiverCoeff * ( a1 * vtmp_n + a2 * vtmp_np1 );
        }
        else
        {
          varAtReceivers( iSeismo, ircv ) = receiverCoeff * ( a1 * vtmp_n + a2 * vtmp_np1 );
        }
        // NOTE: varAtReceivers has size(1) = numReceiversGlobal + 1, this does not OOB
        // left in the forAll loop for sync issues since the following does not depend on `ircv`
        varAtReceivers( iSeismo, nReceivers ) = a1 * time_n + a2 * time_np1;
      }
    } );
  }

  /**
   * @brief Compute the seismo traces for 2d arrays
   * @param[in] time_n Current time iteration
   * @param[in] dt time-step
   * @param[in] regionIndex Index of the current region
   * @param[in] receiverRegion Array containing the region in which the receiver is located
   * @param[in] timeSeismo time when the seismo is computed
   * @param[in] iSeismo i-th seismo trace
   * @param[in] receiverElem Array containing the element on which the receiver is located
   * @param[in] receiverConstants constant part of the receiver term
   * @param[in] receiverIsLocal flag indicating whether the receiver is local or not
   * @param[in] var_np1 Array containing the variable at time n+1
   * @param[in] var_n Array containing the variable at time n
   * @param[out] varAtReceivers Array containing the the variable computed at the receivers
   */
  static void compute2dVariableSeismoTrace( real64 const time_n,
                                            real64 const dt,
                                            localIndex const regionIndex,
                                            arrayView1d< localIndex const > const receiverRegion,
                                            real64 const timeSeismo,
                                            localIndex const iSeismo,
                                            arrayView1d< localIndex const > const receiverElem,
                                            arrayView2d< real64 const > const receiverConstants,
                                            arrayView1d< localIndex const > const receiverIsLocal,
                                            arrayView2d< real32 const > const var_np1,
                                            arrayView2d< real32 const > const var_n,
                                            arrayView2d< real32 > varAtReceivers )
  {
    real64 const time_np1 = time_n + dt;

    real32 const a1 = dt < epsilonLoc ? 1.0 : (time_np1 - timeSeismo) / dt;
    real32 const a2 = 1.0 - a1;

    localIndex const nReceivers = receiverConstants.size( 0 );

    forAll< EXEC_POLICY >( nReceivers, [=] GEOS_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        if( receiverRegion[ircv] == regionIndex )
        {
          real32 vtmp_np1 = 0.0, vtmp_n = 0.0;
          for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
          {
            vtmp_np1 += var_np1( receiverElem[ircv], inode ) * receiverConstants( ircv, inode );
            vtmp_n += var_n( receiverElem[ircv], inode ) * receiverConstants( ircv, inode );
          }
          // linear interpolation between the pressure value at time_n and time_{n+1}
          varAtReceivers( iSeismo, ircv ) = a1 * vtmp_n + a2 * vtmp_np1;
          // NOTE: varAtReceivers has size(1) = numReceiversGlobal + 1, this does not OOB
          // left in the forAll loop for sync issues since the following does not depend on `ircv`
          varAtReceivers( iSeismo, nReceivers ) = a1 * time_n + a2 * time_np1;
        }
      }
    } );
  }

  /**
   * @brief Check if the source point is inside an element or not
   * @param numFacesPerElem number of face on an element
   * @param elemCenter array containing the center of the elements
   * @param faceNormal array containing the normal of all faces
   * @param faceCenter array containing the center of all faces
   * @param elemsToFaces map to get the global faces from element index and local face index
   * @param coords coordinate of the point
   * @return true if coords is inside the element
   */
  GEOS_HOST_DEVICE
  static bool
  locateSourceElement( real64 const numFacesPerElem,
                       real64 const (&elemCenter)[3],
                       arrayView2d< real64 const > const faceNormal,
                       arrayView2d< real64 const > const faceCenter,
                       arraySlice1d< localIndex const > const elemsToFaces,
                       real64 const (&coords)[3] )
  {
    //Loop over the element faces
    real64 tmpVector[3]{};
    for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
    {

      localIndex const iface = elemsToFaces[kfe];
      real64 faceCenterOnFace[3] = {faceCenter[iface][0],
                                    faceCenter[iface][1],
                                    faceCenter[iface][2]};
      real64 faceNormalOnFace[3] = {faceNormal[iface][0],
                                    faceNormal[iface][1],
                                    faceNormal[iface][2]};

      //Test to make sure if the normal is outwardly directed
      LvArray::tensorOps::copy< 3 >( tmpVector, faceCenterOnFace );
      LvArray::tensorOps::subtract< 3 >( tmpVector, elemCenter );
      if( LvArray::tensorOps::AiBi< 3 >( tmpVector, faceNormalOnFace ) < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormalOnFace, -1 );
      }

      // compute the vector face center to query point
      LvArray::tensorOps::subtract< 3 >( faceCenterOnFace, coords );
      localIndex const s = computationalGeometry::sign( LvArray::tensorOps::AiBi< 3 >( faceNormalOnFace, faceCenterOnFace ));

      // all dot products should be non-negative (we enforce outward normals)
      if( s < 0 ) return false;

    }
    return true;
  }

  /**
   * @brief Convert a mesh element point coordinate into a coordinate on the reference element
   * @tparam FE_TYPE finite element type
   * @param[in] coords coordinate of the point
   * @param[in] elemsToNodes map to obtaint global nodes from element index
   * @param[in] nodeCoords array of mesh nodes coordinates
   * @param[out] coordsOnRefElem to contain the coordinate computed in the reference element
   */
  template< typename FE_TYPE >
  GEOS_HOST_DEVICE
  static void
  computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                        arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const elemsToNodes,
                                        arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                                        real64 (& coordsOnRefElem)[3] )
  {
    // only the eight corners of the mesh cell are needed to compute the Jacobian
    real64 xLocal[8][3]{};
    for( localIndex a = 0; a < 8; ++a )
    {
      LvArray::tensorOps::copy< 3 >( xLocal[a], nodeCoords[ elemsToNodes[ FE_TYPE::meshIndexToLinearIndex3D( a )] ] );
    }
    // coordsOnRefElem = invJ*(coords-coordsNode_0)
    real64 invJ[3][3]{};
    FE_TYPE::invJacobianTransformation( 0, xLocal, invJ );
    for( localIndex i = 0; i < 3; ++i )
    {
      // init at (-1,-1,-1) as the origin of the referential elem
      coordsOnRefElem[i] = -1.0;
      for( localIndex j = 0; j < 3; ++j )
      {
        coordsOnRefElem[i] += invJ[i][j] * (coords[j] - xLocal[0][j]);
      }
    }
  }

/**
 * @brief Converts the DAS direction from dip/azimuth to a 3D unit vector
 * @param[in] dip the dip of the linear DAS
 * @param[in] azimuth the azimuth of the linear DAS
 * @param[out] a unit vector pointing in the DAS direction
 */
  GEOS_HOST_DEVICE
  static
  R1Tensor computeDASVector( real64 const dip, real64 const azimuth )
  {
    real64 cd = cos( dip );
    real64 v1 = cd * cos( azimuth );
    real64 v2 = cd * sin( azimuth );
    real64 v3 = sin( dip );
    R1Tensor dasVector = { v1, v2, v3 };
    return dasVector;
  }

};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( WaveSolverUtils::DASType,
              "none",
              "dipole",
              "strainIntegration" );

ENUM_STRINGS( WaveSolverUtils::AttenuationType,
              "none",
              "sls" );

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_ */
