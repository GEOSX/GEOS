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
 * @file BoundaryStencil.hpp
 */

#ifndef GEOS_FINITEVOLUME_BOUNDARYSTENCIL_HPP_
#define GEOS_FINITEVOLUME_BOUNDARYSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geos
{

/**
 * Provides access to the boundary stencil that may be called from a kernel function.
 */
class BoundaryStencilWrapper : public StencilWrapperBase< TwoPointStencilTraits >
{
public:

  /// Coefficient view accessory type
  template< typename VIEWTYPE >
  using CoefficientAccessor = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief Constructor
   * @param elementRegionIndices The container for the element region indices for each point in each stencil
   * @param elementSubRegionIndices The container for the element sub region indices for each point in each stencil
   * @param elementIndices The container for the element indices for each point in each stencil
   * @param weights The container for the weights for each point in each stencil
   * @param faceNormal Face normal vector
   * @param cellToFaceVec Cell center to face center vector
   * @param weightMultiplier Weight multiplier
   */
  BoundaryStencilWrapper( IndexContainerType const & elementRegionIndices,
                          IndexContainerType const & elementSubRegionIndices,
                          IndexContainerType const & elementIndices,
                          WeightContainerType const & weights,
                          arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & faceNormal,
                          arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & cellToFaceVec,
                          arrayView1d< real64 > const & weightMultiplier );

  /**
   * @brief Compute weights and derivatives w.r.t to the coefficient.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[out] weight view weights
   * @param[out] dWeight_dCoef derivative of the weights w.r.t to the coefficient
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                       real64 & weight,
                       real64 ( &dWeight_dCoef )[3] ) const;

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex size() const
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr localIndex stencilSize( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }

  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  constexpr localIndex numPointsInFlux( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxNumPointsInFlux;
  }

  void getFaceNormal( localIndex const iconn, real64 (& faceNormal)[3] ) const
  {
    GEOS_UNUSED_VAR( iconn, faceNormal );
  }

private:

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > m_faceNormal;
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > m_cellToFaceVec;
  arrayView1d< real64 > m_weightMultiplier;
};

/**
 * @brief Provides management of the boundary stencil points
 * (stencils used to prescribe boundary conditions on domain boundaries, i.e. faces)
 */
class BoundaryStencil final : public StencilBase< TwoPointStencilTraits, BoundaryStencil >
{
public:

  /**
   * @brief Constructor.
   */
  BoundaryStencil();

  /**
   * @brief Defines the order of element/face in the stencil.
   */
  struct Order
  {
    static constexpr localIndex ELEM = 0; ///< Order of element index in stencil
    static constexpr localIndex FACE = 1; ///< Order of face index in stencil
  };

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override;

  /**
   * @brief Add the vectors need to compute the transmissiblity to the Stencil.
   * @param[in] transMultiplier the transmissibility multiplier
   * @param[in] faceNormal the normal to the face
   * @param[in] cellToFaceVec distance vector between the cell center and the face
   */
  void addVectors( real64 const & transMultiplier,
                   real64 const (&faceNormal)[3],
                   real64 const (&cellToFaceVec)[3] );

  /**
   * @copydoc StencilBase<BoundaryStencilTraits,BoundaryStencil>::size
   */
  virtual localIndex size() const override
  {
    return m_elementRegionIndices.size( 0 );
  }

  /**
   * @brief Gives the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BoundaryStencilWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

private:

  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_faceNormal;
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_cellToFaceVec;
  array1d< real64 > m_weightMultiplier;

};

GEOS_HOST_DEVICE
inline void
BoundaryStencilWrapper::
  computeWeights( localIndex const iconn,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                  real64 & weight,
                  real64 (& dWeight_dCoef)[3] ) const
{
  localIndex const er  = m_elementRegionIndices[iconn][BoundaryStencil::Order::ELEM];
  localIndex const esr = m_elementSubRegionIndices[iconn][BoundaryStencil::Order::ELEM];
  localIndex const ei  = m_elementIndices[iconn][BoundaryStencil::Order::ELEM];

  // Preload data onto the stack
  real64 faceNormal[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_faceNormal[iconn] );
  real64 const cellToFace[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_cellToFaceVec[iconn] );
  arraySlice1d< real64 const > const coef = coefficient[er][esr][ei][0];

  // Compute the face-cell transmissibility
  auto const compute = [&]
  {
    real64 faceConormal[3];
    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, coef, faceNormal );
    weight = LvArray::tensorOps::AiBi< 3 >( cellToFace, faceConormal );
    LvArray::tensorOps::hadamardProduct< 3 >( dWeight_dCoef, cellToFace, faceNormal );
  };

  compute();
  if( weight < 0.0 )
  {
    // Correct negative weight issue arising from non-K-orthogonal grids
    // by replacing face normal with c2f vector
    LvArray::tensorOps::copy< 3 >( faceNormal, cellToFace );
    compute();
  }

  // Scale by face/distance and trans multiplier
  real64 const mult = m_weights[iconn][BoundaryStencil::Order::ELEM] * m_weightMultiplier[iconn];
  weight *= mult;
  LvArray::tensorOps::scale< 3 >( dWeight_dCoef, mult );
}

} // namespace geos

#endif // GEOS_FINITEVOLUME_BOUNDARYSTENCIL_HPP_
