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
 * @file QuadratureWeights.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_WEIGHTS_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_WEIGHTS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/TeamKernelInterface/TeamKernelFunctions/common.hpp"

namespace geosx
{

namespace stackVariables
{

template < localIndex num_quads_1d >
struct StackQuadratureWeights;

template < >
struct StackQuadratureWeights< 2 >
{
  // real64 weights[num_quads_1d];

  GEOSX_HOST_DEVICE
  StackQuadratureWeights( RAJA::LaunchContext & ctx )
  {
    // Initialize quadrature weights
    // TODO generalize/use threads
    // weights[0] = 1.0;
    // weights[1] = 1.0;
  }

  GEOSX_HOST_DEVICE
  real64 operator()( TensorIndex const & quad_index )
  {
    return 1.0; //(*weights)[ quad_index.x ] * (*weights)[ quad_index.y ] * (*weights)[ quad_index.z ];
  }
};

template < localIndex num_quads_1d >
struct SharedQuadratureWeights
{
  real64 ( * weights )[num_quads_1d];

  GEOSX_HOST_DEVICE
  SharedQuadratureWeights( RAJA::LaunchContext & ctx )
  {
    // Initialize quadrature weights
    GEOSX_STATIC_SHARED real64 s_weights[num_quads_1d];
    weights = &s_weights;
    // TODO generalize/use threads
    s_weights[0] = 1.0;
    s_weights[1] = 1.0;
  }

  GEOSX_HOST_DEVICE
  real64 operator()( TensorIndex const & quad_index )
  {
    return (*weights)[ quad_index.x ] * (*weights)[ quad_index.y ] * (*weights)[ quad_index.z ];
  }
};

} // namespace stackVariables

} // namespace geosx

#endif // GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_WEIGHTS_HPP_
