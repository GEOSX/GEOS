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
 * @file CellElementStencilTPFA.hpp
 */

#ifndef GEOS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_
#define GEOS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_

#include "StencilBase.hpp"

namespace geos
{

/**
 * Provides access to the cellElement stencil that may be called from a kernel function.
 */
class CellElementStencilTPFAWrapper : public StencilWrapperBase< TwoPointStencilTraits >
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
   * @param transMultiplier Transmissibility multiplier
   * @param geometricStabilizationCoef Geometric coefficient for pressure jump stabilization
   */
  CellElementStencilTPFAWrapper( IndexContainerType const & elementRegionIndices,
                                 IndexContainerType const & elementSubRegionIndices,
                                 IndexContainerType const & elementIndices,
                                 WeightContainerType const & weights,
                                 arrayView2d< real64 > const & faceNormal,
                                 arrayView3d< real64 > const & cellToFaceVec,
                                 arrayView1d< real64 > const & transMultiplier,
                                 arrayView1d< real64 > const & geometricStabilizationCoef );

  /**
   * @brief Compute weights and derivatives w.r.t to one variable.
   * @param[in] iconn connection index
   * @param[in] coefficient view accessor to the coefficient used to compute the weights
   * @param[in] dCoeff_dVar view accessor to the derivative of the coefficient w.r.t to the variable
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weights w.r.t to the variable
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex const iconn,
                       CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                       CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;

  /**
   * @brief Compute weights and derivatives w.r.t to one variable without coefficient
   * Used in ReactiveCompositionalMultiphaseOBL solver for thermal transmissibility computation:
   * here, conductivity is a part of operator and cannot be used directly as a coefficient
   * @param[in] iconn connection index
   * @param[out] weight view weights
   * @param[out] dWeight_dVar derivative of the weights w.r.t to the variable
   */
  GEOS_HOST_DEVICE
  void computeWeights( localIndex iconn,
                       real64 ( &weight )[1][2],
                       real64 ( &dWeight_dVar )[1][2] ) const;

  /**
   * @brief Compute the stabilization weights
   * @param[in] iconn connection index
   * @param[out] stabilizationWeight view weights
   */
  GEOS_HOST_DEVICE
  void computeStabilizationWeights( localIndex iconn,
                                    real64 ( &stabilizationWeight )[1][2] ) const;

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  localIndex size() const
  {
    return m_elementRegionIndices.size( 0 );
  }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  localIndex stencilSize( localIndex const index ) const
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
  localIndex numPointsInFlux( localIndex const index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxNumPointsInFlux;
  }

  void getFaceNormal( localIndex const iconn, real64 (& faceNormal)[3] ) const
  {
    LvArray::tensorOps::copy< 3 >( faceNormal, m_faceNormal[iconn] );
  }
private:

  arrayView2d< real64 > m_faceNormal;
  arrayView3d< real64 > m_cellToFaceVec;
  arrayView1d< real64 > m_transMultiplier;
  arrayView1d< real64 > m_geometricStabilizationCoef;
};


/**
 * @class CellElementStencilTPFA
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
class CellElementStencilTPFA final : public StencilBase< TwoPointStencilTraits, CellElementStencilTPFA >
{
public:

  /**
   * @brief Default constructor.
   */
  CellElementStencilTPFA();

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override;

  /**
   * @brief Add the vectors need to compute the transmissiblity to the Stencil.
   * @param[in] transMultiplier the transmissibility multiplier
   * @param[in] geometricStabilizationCoef the stabilization weight
   * @param[in] faceNormal the normal to the face
   * @param[in] cellToFaceVec distance vector between the cell center and the face
   */
  void addVectors( real64 const & transMultiplier,
                   real64 const & geometricStabilizationCoef,
                   real64 const (&faceNormal)[3],
                   real64 const (&cellToFaceVec)[2][3] );

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Reserve the size of the stencil
   * @param[in] size the size of the stencil to reserve
   */
  virtual void reserve( localIndex const size ) override;

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxStencilSize;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CellElementStencilTPFAWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;



private:

  array2d< real64 > m_faceNormal;
  array3d< real64 > m_cellToFaceVec;
  array1d< real64 > m_transMultiplier;
  array1d< real64 > m_geometricStabilizationCoef;
};

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::
  computeWeights( localIndex const iconn,
                  CoefficientAccessor< arrayView3d< real64 const > > const & coefficient,
                  CoefficientAccessor< arrayView3d< real64 const > > const & dCoeff_dVar,
                  real64 (& weight)[1][2],
                  real64 (& dWeight_dVar )[1][2] ) const
{
  real64 halfWeight[2];
  real64 dHalfWeight_dVar[2];

  // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex i = 0; i < 2; ++i )
  {
    localIndex const er  = m_elementRegionIndices[iconn][i];
    localIndex const esr = m_elementSubRegionIndices[iconn][i];
    localIndex const ei  = m_elementIndices[iconn][i];

    halfWeight[i] = m_weights[iconn][i];
    dHalfWeight_dVar[i] = m_weights[iconn][i];

    // Proper computation
    real64 faceNormal[3];
    LvArray::tensorOps::copy< 3 >( faceNormal, m_faceNormal[iconn] );
    if( LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    real64 faceConormal[3];
    real64 dFaceConormal_dVar[3];
    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, coefficient[er][esr][ei][0], faceNormal );
    LvArray::tensorOps::hadamardProduct< 3 >( dFaceConormal_dVar, dCoeff_dVar[er][esr][ei][0], faceNormal );
    halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceConormal );
    dHalfWeight_dVar[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], dFaceConormal_dVar );

    // correct negative weight issue arising from non-K-orthogonal grids
    if( halfWeight[i] < 0.0 )
    {
      LvArray::tensorOps::hadamardProduct< 3 >( faceConormal,
                                                coefficient[er][esr][ei][0],
                                                m_cellToFaceVec[iconn][i] );
      LvArray::tensorOps::hadamardProduct< 3 >( dFaceConormal_dVar,
                                                dCoeff_dVar[er][esr][ei][0],
                                                m_cellToFaceVec[iconn][i] );
      halfWeight[i] = m_weights[iconn][i];
      dHalfWeight_dVar[i] = m_weights[iconn][i];
      halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceConormal );
      dHalfWeight_dVar[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], dFaceConormal_dVar );
    }
  }

  // Do harmonic and arithmetic averaging
  real64 const product = halfWeight[0]*halfWeight[1];
  real64 const sum = halfWeight[0]+halfWeight[1];

  real64 const harmonicWeight   = sum > 0 ? product / sum : 0.0;
  real64 const arithmeticWeight = sum / 2;

  real64 dHarmonicWeight_dVar[2];
  real64 dArithmeticWeight_dVar[2];

  dHarmonicWeight_dVar[0] = sum > 0 ? (dHalfWeight_dVar[0]*sum*halfWeight[1] - dHalfWeight_dVar[0]*halfWeight[0]*halfWeight[1]) / ( sum*sum ) : 0.0;
  dHarmonicWeight_dVar[1] = sum > 0 ? (dHalfWeight_dVar[1]*sum*halfWeight[0] - dHalfWeight_dVar[1]*halfWeight[1]*halfWeight[0]) / ( sum*sum ) : 0.0;

  dArithmeticWeight_dVar[0] = dHalfWeight_dVar[0] / 2;
  dArithmeticWeight_dVar[1] = dHalfWeight_dVar[1] / 2;

  real64 const meanPermCoeff = 1.0; //TODO make it a member if it is really necessary

  real64 const value = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    weight[0][ke] = m_transMultiplier[iconn] * value * (ke == 0 ? 1 : -1);

    real64 const dValue_dVar = meanPermCoeff * dHarmonicWeight_dVar[ke] + (1 - meanPermCoeff) * dArithmeticWeight_dVar[ke];
    dWeight_dVar[0][ke] = m_transMultiplier[iconn] * dValue_dVar;
  }
}

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::
  computeWeights( localIndex iconn,
                  real64 (& weight)[1][2],
                  real64 (& dWeight_dVar )[1][2] ) const
{
  real64 halfWeight[2];

  // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex i =0; i<2; i++ )
  {
    halfWeight[i] = m_weights[iconn][i];

    // Proper computation
    real64 faceNormal[3];

    LvArray::tensorOps::copy< 3 >( faceNormal, m_faceNormal[iconn] );

    if( LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceNormal );

    // correct negative weight issue arising from non-K-orthogonal grids
    if( halfWeight[i] < 0.0 )
    {
      halfWeight[i] = m_weights[iconn][i];
      halfWeight[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], m_cellToFaceVec[iconn][i] );
    }
  }

  // Do harmonic and arithmetic averaging
  real64 const product = halfWeight[0]*halfWeight[1];
  real64 const sum = halfWeight[0]+halfWeight[1];

  real64 const harmonicWeight   = sum > 0 ? product / sum : 0.0;
  real64 const arithmeticWeight = sum / 2;

  real64 const meanPermCoeff = 1.0; //TODO make it a member if it is really necessary

  real64 const value = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    weight[0][ke] = m_transMultiplier[iconn] * value * (ke == 0 ? 1 : -1);
  }

  dWeight_dVar[0][0] = 0.0;
  dWeight_dVar[0][1] = 0.0;
}

GEOS_HOST_DEVICE
inline void
CellElementStencilTPFAWrapper::
  computeStabilizationWeights( localIndex iconn,
                               real64 ( & stabilizationWeight )[1][2] ) const
{
  stabilizationWeight[0][0] = m_geometricStabilizationCoef[iconn];
  stabilizationWeight[0][1] = -m_geometricStabilizationCoef[iconn];
}


} /* namespace geos */

#endif /* GEOS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_ */
