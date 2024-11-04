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
 * @file MultilinearInterpolatorAdaptiveKernels.hpp
 */

#ifndef GEOS_FUNCTIONS_MULTILINEARINTERPOLATORADAPTIVEKERNELS_HPP_
#define GEOS_FUNCTIONS_MULTILINEARINTERPOLATORADAPTIVEKERNELS_HPP_

#include "MultilinearInterpolatorBaseKernels.hpp"

namespace geos
{
/**
 * @class KernelWrapper
 *
 *
 */

/**
 * @class MultilinearInterpolatorAdaptiveKernel
 *
 * A class for multilinear piecewise interpolation with access to exact function evaluator
 * All functions are interpolated using the same uniformly discretized space
 *
 * @tparam NUM_DIMS number of dimensions (inputs)
 * @tparam NUM_OPS number of interpolated functions (outputs)
 */
template< integer NUM_DIMS, integer NUM_OPS, typename INDEX_T = __uint128_t >
class MultilinearInterpolatorAdaptiveKernel : public MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >
{
public:
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::numDims;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::numOps;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::numVerts;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::m_axisMinimums;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::m_axisMaximums;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::m_axisPoints;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::m_axisSteps;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::m_axisStepInvs;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::m_axisHypercubeMults;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::m_axisPointMults;
  using MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::m_coordinates;
  using typename MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS, INDEX_T >::longIndex;

  /**
   * @brief Construct a new Multilinear Interpolator Adaptive Kernel object
   * @param[in] evalFunc evaluator of exact function values
   */
  MultilinearInterpolatorAdaptiveKernel( PythonFunction< longIndex > * evalFunc ):
    MultilinearInterpolatorBaseKernel< NUM_DIMS, NUM_OPS >( evalFunc->getAxisMinimums(),
                                                            evalFunc->getAxisMaximums(),
                                                            evalFunc->getAxisPoints(),
                                                            evalFunc->getAxisSteps(),
                                                            evalFunc->getAxisStepInvs(),
                                                            evalFunc->getAxisHypercubeMults() ),
    m_evalFunc( evalFunc )
  {};
protected:
  /**
   * @brief Get values of operators at a given point
   * Provide a reference to correct location in the adaptive point storage.
   * If the point is not found, compute it first, and then return the reference.
   *
   * @param[in] pointIndex index of point
   * @return operator values at given point
   */
  GEOS_HOST_DEVICE
  array1d< real64 > const &
  getPointData( longIndex const pointIndex ) const
  {
    auto & m_pointData = m_evalFunc->getPointData();
    auto item = m_pointData.find( pointIndex );
    if( item == m_pointData.end())
    {
      array1d< real64 > newPoint ( numOps );
      this->getPointCoordinates( pointIndex, m_coordinates );
      m_evalFunc->evaluate( m_coordinates, newPoint );
      m_pointData[pointIndex] = newPoint;
      return m_pointData.at( pointIndex );
    }
    else
      return item->second;
  }
  /**
   * @brief Get indexes of all vertices for given hypercube
   *
   * @param[in] hypercubeIndex index of the hyporcube
   * @param[out] hypercubePoints indexes of all vertices of hypercube
   */
  GEOS_HOST_DEVICE
  inline
  void
  getHypercubePoints( longIndex const hypercubeIndex,
                      array1d< longIndex > & hypercubePoints ) const
  {
    longIndex axisIdx, remainderIdx = hypercubeIndex;
    integer pwr = numVerts;

    for( integer i = 0; i < numVerts; ++i )
      hypercubePoints[i] = 0;

    for( integer i = 0; i < numDims; ++i )
    {
      axisIdx = remainderIdx / m_axisHypercubeMults[i];
      remainderIdx = remainderIdx % m_axisHypercubeMults[i];
      pwr /= 2;

      for( integer j = 0; j < numVerts; ++j )
      {
        integer zeroOrOne = (j / pwr) % 2;
        hypercubePoints[j] += (axisIdx + zeroOrOne) * m_axisPointMults[i];
      }
    }
  }
  /**
   * @brief Get pointer to hypercube data
   *
   * @param[in] hypercubeIndex
   * @return pointer to hypercube data
   */
  GEOS_HOST_DEVICE
  inline
  real64 const *
  getHypercubeData( longIndex const hypercubeIndex ) const override
  {
    auto & m_hypercubeData = m_evalFunc->getHypercubeData();
    auto item = m_hypercubeData.find( hypercubeIndex );
    if( item == m_hypercubeData.end())
    {
      array1d< longIndex > points ( numVerts );
      array1d< real64 > newHypercube ( numVerts * numOps );

      this->getHypercubePoints( hypercubeIndex, points );

      for( integer i = 0; i < numVerts; ++i )
      {
        // obtain point data and copy it to hypercube data
        auto const & pData = this->getPointData( points[i] );
        for( integer op = 0; op < numOps; op++ )
        {
          newHypercube[i * numOps + op] = pData[op];
        }
      }
      m_hypercubeData[hypercubeIndex] = newHypercube;
      return m_hypercubeData[hypercubeIndex].data();
    }
    else
      return item->second.data();
  }
  /**
   * @brief wrapper providing interface to Python function
   */
  PythonFunction< longIndex > * m_evalFunc;
};

} /* namespace geos */

#endif /* GEOS_FUNCTIONS_MULTILINEARINTERPOLATORADAPTIVEKERNELS_HPP_ */
