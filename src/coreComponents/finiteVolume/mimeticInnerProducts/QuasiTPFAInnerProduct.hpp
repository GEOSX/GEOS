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
 * @file QuasiTPFAInnerProduct.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_QUASITPFAINNERPRODUCT_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_QUASITPFAINNERPRODUCT_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @class QuasiTPFAInnerProduct
 *
 * Provides an implementation of a quasi-TPFA inner product in the hybrid FVM solvers
 */
class QuasiTPFAInnerProduct : public MimeticInnerProductBase
{
public:

  /**
   * @brief In a given element, recompute the transmissibility matrix using the quasi TPFA inner product.
   * @param[in] nodePosition the position of the nodes
   * @param[in] transMultiplier the transmissibility multipliers at the mesh faces
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   *
   * @details Reference: K-A Lie, An Introduction to Reservoir Simulation Using MATLAB/GNU Octave (2019)
   */
  template< localIndex NF >
  GEOS_HOST_DEVICE
  static void
  compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           arrayView1d< real64 const > const & transMultiplier,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const & elemToFaces,
           arraySlice1d< real64 const > const & elemCenter,
           real64 const & elemVolume,
           real64 const (&elemPerm)[ 3 ],
           real64 const & lengthTolerance,
           arraySlice2d< real64 > const & transMatrix );

};

template< localIndex NF >
GEOS_HOST_DEVICE
void
QuasiTPFAInnerProduct::compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                arrayView1d< real64 const > const & transMultiplier,
                                ArrayOfArraysView< localIndex const > const & faceToNodes,
                                arraySlice1d< localIndex const > const & elemToFaces,
                                arraySlice1d< real64 const > const & elemCenter,
                                real64 const & elemVolume,
                                real64 const (&elemPerm)[ 3 ],
                                real64 const & lengthTolerance,
                                arraySlice2d< real64 > const & transMatrix )
{
  MimeticInnerProductBase::computeParametricInnerProduct< NF >( nodePosition,
                                                                transMultiplier,
                                                                faceToNodes,
                                                                elemToFaces,
                                                                elemCenter,
                                                                elemVolume,
                                                                elemPerm,
                                                                2.0,
                                                                lengthTolerance,
                                                                transMatrix );
}

} // end namespace mimeticInnerProduct

} // end namespace geos

#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_QUASITPFAINNERPRODUCT_HPP_
