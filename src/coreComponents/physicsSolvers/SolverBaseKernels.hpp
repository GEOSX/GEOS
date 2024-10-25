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
 * @file SolverBaseKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLVERBASEKERNELS_HPP
#define GEOS_PHYSICSSOLVERS_SOLVERBASEKERNELS_HPP

#include "codingUtilities/EnumStrings.hpp"
#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{

namespace solverBaseKernels
{

/******************************** ResidualNormKernelBase ********************************/

/**
 * @tparam NUM_NORM number of norms to compute (NUM_NORM is related to the number of equations to solve)
 * For instance, NUM_NORM=1 for isothermal simulations and NUM_NORM=2 for thermal simulations
 * @brief Define the base interface for the residual calculations
 */
template< integer NUM_NORM >
class ResidualNormKernelBase
{
public:


  /// Compile time value for the number of norms to compute
  static constexpr integer numNorm = NUM_NORM;

  ResidualNormKernelBase( globalIndex const rankOffset,
                          arrayView1d< real64 const > const & localResidual,
                          arrayView1d< globalIndex const > const & dofNumber,
                          arrayView1d< localIndex const > const & ghostRank,
                          real64 const minNormalizer ):
    m_rankOffset( rankOffset ),
    m_localResidual( localResidual ),
    m_dofNumber( dofNumber ),
    m_ghostRank( ghostRank ),
    m_minNormalizer( minNormalizer )
  {}

  /**
   * @struct LinfStackVariables
   * @brief Kernel variables located on the stack for Linf norm
   */
  struct LinfStackVariables
  {
    /// Index of the local row in the residual vector
    localIndex localRow;

    /// Normalized residual value for the element/node/face
    real64 localValue[numNorm]{};
  };

  /**
   * @struct L2StackVariables
   * @brief Kernel variables located on the stack for L2 norm
   */
  struct L2StackVariables : public LinfStackVariables
  {
    /// Normalizer value for the element/node/face
    real64 localNormalizer[numNorm]{};
  };


  /**
   * @brief Getter for the ghost rank
   * @param[in] i the looping index of the element/node/face
   * @return the ghost rank of the element/node/face
   */
  GEOS_HOST_DEVICE
  integer ghostRank( localIndex const i ) const
  { return m_ghostRank( i ); }

  /**
   * @brief Setup the residual Linf normal calculations
   * @param[in] i the element/node/face index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  virtual void setupLinf( localIndex const i,
                          LinfStackVariables & stack ) const
  {
    stack.localRow = m_dofNumber[i] - m_rankOffset;
  }

  /**
   * @brief Setup the residual L2 normal calculations
   * @param[in] i the element/node/face index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  virtual void setupL2( localIndex const i,
                        L2StackVariables & stack ) const
  {
    stack.localRow = m_dofNumber[i] - m_rankOffset;
  }


  /**
   * @brief Compute the local values for the Linf norm
   * @param[in] i the element/node/face index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const i,
                            LinfStackVariables & stack ) const = 0;

  /**
   * @brief Compute the local values and normalizer for the L2 norm
   * @param[in] i the element/node/face index
   * @param[inout] stack the stack variables
   */
  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const i,
                          L2StackVariables & stack ) const = 0;

  /**
   * @brief Performs the kernel launch for the L-\infty norm
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] size the number of elements/nodes/faces
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   * @param[inout] residualNorms the norms to compute
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launchLinf( localIndex const size,
              KERNEL_TYPE const & kernelComponent,
              real64 (& residualNorm)[numNorm] )
  {
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > localResidualNorm[numNorm]{};

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      if( kernelComponent.ghostRank( i ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::LinfStackVariables stack;
      kernelComponent.setupLinf( i, stack );
      kernelComponent.computeLinf( i, stack );

      for( integer j = 0; j < numNorm; ++j )
      {
        localResidualNorm[j].max( LvArray::math::abs( stack.localValue[j] ) );
      }
    } );

    for( integer j = 0; j < numNorm; ++j )
    {
      residualNorm[j] = localResidualNorm[j].get();
    }
  }

  /**
   * @brief Performs the kernel launch for the L2 norm
   * @tparam POLICY the policy used in the RAJA kernels
   * @tparam KERNEL_TYPE the kernel type
   * @param[in] size the number of elements/nodes/faces
   * @param[inout] kernelComponent the kernel component providing access to the compute function
   * @param[inout] residualNorm the norms to compute
   * @param[inout] residualNormalizer the norms to compute
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static void
  launchL2( localIndex const size,
            KERNEL_TYPE const & kernelComponent,
            real64 (& residualNorm)[numNorm],
            real64 (& residualNormalizer)[numNorm] )
  {
    RAJA::ReduceSum< ReducePolicy< POLICY >, real64 > localResidualNorm[numNorm]{};
    RAJA::ReduceSum< ReducePolicy< POLICY >, real64 > localResidualNormalizer[numNorm]{};

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      if( kernelComponent.ghostRank( i ) >= 0 )
      {
        return;
      }

      typename KERNEL_TYPE::L2StackVariables stack;
      kernelComponent.setupL2( i, stack );
      kernelComponent.computeL2( i, stack );

      for( integer j = 0; j < numNorm; ++j )
      {
        localResidualNorm[j] += stack.localValue[j];
        localResidualNormalizer[j] += stack.localNormalizer[j];
      }
    } );

    for( integer j = 0; j < numNorm; ++j )
    {
      residualNorm[j] = localResidualNorm[j].get();
      residualNormalizer[j] = localResidualNormalizer[j].get();
    }
  }


protected:

  /// Offset for my MPI rank
  globalIndex const m_rankOffset;

  /// View on the local residual
  arrayView1d< real64 const > const m_localResidual;

  /// View on the dof numbers
  arrayView1d< globalIndex const > const m_dofNumber;

  /// View on the ghost ranks
  arrayView1d< integer const > const m_ghostRank;

  /// Value used to make sure that normalizers are never zero
  real64 const m_minNormalizer;

};

/**
 * @class LinfResidualNormHelper
 * @brief Utility class to compute the global Linf residual norm
 */
class LinfResidualNormHelper
{
public:

  template< integer NUM_NORM >
  static void updateLocalNorm( real64 const (&subRegionResidualNorm)[NUM_NORM],
                               array1d< real64 > & localResidualNorm )
  {
    for( integer i = 0; i < NUM_NORM; ++i )
    {
      if( subRegionResidualNorm[i] > localResidualNorm[i] )
      {
        localResidualNorm[i] = subRegionResidualNorm[i];
      }
    }
  }

  static void computeGlobalNorm( real64 const & localResidualNorm,
                                 real64 & globalResidualNorm )
  {
    globalResidualNorm = MpiWrapper::max( localResidualNorm );
  }

  static void computeGlobalNorm( array1d< real64 > const & localResidualNorm,
                                 array1d< real64 > & globalResidualNorm )
  {
    MpiWrapper::allReduce( localResidualNorm.data(),
                           globalResidualNorm.data(),
                           localResidualNorm.size(),
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Max ),
                           MPI_COMM_GEOS );
  }
};

/**
 * @class L2ResidualNormHelper
 * @brief Utility class to compute the global L2 residual norm
 */
class L2ResidualNormHelper
{
public:

  template< integer NUM_NORM >
  static void updateLocalNorm( real64 const (&subRegionResidualNorm)[NUM_NORM],
                               real64 const (&subRegionResidualNormalizer)[NUM_NORM],
                               array1d< real64 > & localResidualNorm,
                               array1d< real64 > & localResidualNormalizer )
  {
    for( integer i = 0; i < NUM_NORM; ++i )
    {
      localResidualNorm[i] += subRegionResidualNorm[i];
      localResidualNormalizer[i] += subRegionResidualNormalizer[i];
    }
  }

  static void computeGlobalNorm( real64 const & localResidualNorm,
                                 real64 const & localResidualNormalizer,
                                 real64 & globalResidualNorm )
  {
    globalResidualNorm = sqrt( MpiWrapper::sum( localResidualNorm ) ) / sqrt( MpiWrapper::sum( localResidualNormalizer ) );
  }

  static void computeGlobalNorm( array1d< real64 > const & localResidualNorm,
                                 array1d< real64 > const & localResidualNormalizer,
                                 array1d< real64 > & globalResidualNorm )
  {
    array1d< real64 > sumLocalResidualNorm( localResidualNorm.size() );
    array1d< real64 > sumLocalResidualNormalizer( localResidualNormalizer.size() );
    MpiWrapper::allReduce( localResidualNorm.data(),
                           sumLocalResidualNorm.data(),
                           localResidualNorm.size(),
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                           MPI_COMM_GEOS );
    MpiWrapper::allReduce( localResidualNormalizer.data(),
                           sumLocalResidualNormalizer.data(),
                           localResidualNormalizer.size(),
                           MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
                           MPI_COMM_GEOS );
    for( integer i = 0; i < localResidualNorm.size(); ++i )
    {
      globalResidualNorm[i] = sqrt( sumLocalResidualNorm[i] ) / sqrt( sumLocalResidualNormalizer[i] );
    }
  }

};

/**
 * @brief Type of norm used to check convergence
 * TODO: find a way to put this inside the class
 */
enum class NormType : integer
{
  Linf,  /**< Linfinity norm */
  L2     /**< L2 */
};

ENUM_STRINGS( NormType,
              "Linfinity",
              "L2" );


} // namespace solverBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_SOLVERBASEKERNELS_HPP
