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
 * @file applyTestFunctions/distributed_3d.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_3D_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_3D_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "tensor/tensor_types.hpp"

namespace geosx
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

namespace impl
{

template < localIndex... Sizes >
using SharedTensor = tensor::StaticPointerDTensor< Sizes... >;

// 3D Distributed

/**
 * @brief "Apply" test functions.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] q_values Values of the "qFunction" at quadrature points.
 * @return Contribution of the q_values to the degrees of freedom.
 */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyTestFunctions( StackVariables & stack,
                         real64 const (& basis)[num_dofs_1d][num_quads_1d],
                         tensor::Static3dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d > const & q_values,
                         tensor::Static3dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d > & dofs )
{
  RAJA::LaunchContext & ctx = stack.ctx;

  SharedTensor< num_quads_1d, num_quads_1d, num_quads_1d > Du( stack.shared_mem[0] );
  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex quad_x, localIndex quad_y, localIndex quad_z){
    Du( quad_x, quad_y, quad_z ) = q_values( quad_x, quad_y, quad_z );
  } );

  // Load in registers values for (quad_x, quad_y, quad_z)
  // FIXME: not interesting for low order and if basis already in stack mem
  real64 Bx[num_quads_1d], By[num_quads_1d], Bz[num_quads_1d];
  localIndex dof_x = stack.tidx;
  if ( dof_x < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
    {
      Bx[quad_x] = basis[dof_x][quad_x];
    }
  }
  localIndex dof_y = stack.tidy;
  if ( dof_y < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
    {
      By[quad_y] = basis[dof_y][quad_y];
    }
  }
  localIndex dof_z = stack.tidz;
  if ( dof_z < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      Bz[quad_z] = basis[dof_z][quad_z];
    }
  }

  ctx.teamSync();

  loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
          [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z){
    real64 bbbu = 0.0;
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      real64 const bz = Bz[ quad_z ];
      #pragma unroll
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const by = By[ quad_y ];
        #pragma unroll
        for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
        {
          real64 const bx = Bx[ quad_x ];
          real64 const val = Du( quad_x, quad_y, quad_z );
          bbbu += bx * by * bz * val;
        }
      }
    }
    q_values( dof_x, dof_y, dof_z ) = bbbu;
  } );
}

/**
 * @brief "Apply" test functions.
 * 
 * @param[in] basis Operator returning the shape functions values at the desired
 *                  quadrature points.
 *                  `basis ( dof, quad ) = phi_dof ( x_quad )`
 * @param[in] q_values Values of the "qFunction" at quadrature points.
 * @return Contribution of the q_values to the degrees of freedom.
 */
template < typename StackVariables,
           localIndex num_dofs_1d,
           localIndex num_quads_1d,
           localIndex num_comp >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void applyTestFunctions(
  StackVariables & stack,
  real64 const (& basis)[num_dofs_1d][num_quads_1d],
  tensor::Static3dThreadDTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp > const & q_values,
  tensor::Static3dThreadDTensor< num_dofs_1d, num_dofs_1d, num_dofs_1d, num_comp > & dofs )
{
  RAJA::LaunchContext & ctx = stack.ctx;

  // Load in registers values for (quad_x, quad_y, quad_z)
  // FIXME: not interesting for low order and if basis already in stack mem
  real64 Bx[num_quads_1d], By[num_quads_1d], Bz[num_quads_1d];
  localIndex dof_x = stack.tidx;
  if ( dof_x < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
    {
      Bx[quad_x] = basis[dof_x][quad_x];
    }
  }
  localIndex dof_y = stack.tidy;
  if ( dof_y < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
    {
      By[quad_y] = basis[dof_y][quad_y];
    }
  }
  localIndex dof_z = stack.tidz;
  if ( dof_z < num_dofs_1d )
  {
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      Bz[quad_z] = basis[dof_z][quad_z];
    }
  }

  SharedTensor< num_quads_1d, num_quads_1d, num_quads_1d, num_comp > Du( stack.shared_mem[0] );
  loop3D( stack, num_quads_1d, num_quads_1d, num_quads_1d,
          [&]( localIndex quad_x, localIndex quad_y, localIndex quad_z ){
    #pragma unroll
    for (size_t comp = 0; comp < num_comp; comp++)
    {
      Du( quad_x, quad_y, quad_z, comp ) = q_values( quad_x, quad_y, quad_z, comp );
    }
  } );

  ctx.teamSync();

  loop3D( stack, num_dofs_1d, num_dofs_1d, num_dofs_1d,
          [&]( localIndex dof_x, localIndex dof_y, localIndex dof_z ){
    real64 bbbu[ num_comp ];
    #pragma unroll
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      bbbu[ comp ] = 0.0;
    }
    #pragma unroll
    for (localIndex quad_z = 0; quad_z < num_quads_1d; quad_z++)
    {
      real64 const bz = Bz[quad_z];
      #pragma unroll
      for (localIndex quad_y = 0; quad_y < num_quads_1d; quad_y++)
      {
        real64 const by = By[quad_y];
        #pragma unroll
        for (localIndex quad_x = 0; quad_x < num_quads_1d; quad_x++)
        {
          real64 const bx = Bx[quad_x];
          real64 const b = bx * by * bz;
          #pragma unroll
          for (localIndex comp = 0; comp < num_comp; comp++)
          {
            real64 const val = Du( quad_x, quad_y, quad_z, comp );
            bbbu[ comp ] += b * val;
          }
        }
      }
    }
    #pragma unroll
    for (localIndex comp = 0; comp < num_comp; comp++)
    {
      dofs( dof_x, dof_y, dof_z, comp ) = bbbu[ comp ];
    }
  } );
}

} // namespace impl

} // namespace finiteElement

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_APPLY_TEST_3D_HPP_ */
