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

#include "ImplicitKernelBase.hpp"

/**
 * @file SparsityKernelBase.hpp
 */

#ifndef GEOS_FINITEELEMENT_SPARSITYKERNELBASE_HPP_
#define GEOS_FINITEELEMENT_SPARSITYKERNELBASE_HPP_



namespace geos
{

namespace finiteElement
{

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @class SparsityKernelBase
 * @brief Define the base interface for implicit finite element kernels.
 * @copydoc geos::finiteElement::KernelBase
 *
 * ### SparsityKernelBase Description
 * Provide common kernels for generation of CRS Sparsity patterns.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class SparsityKernelBase : public ImplicitKernelBase< SUBREGION_TYPE,
                                                      CONSTITUTIVE_TYPE,
                                                      FE_TYPE,
                                                      NUM_DOF_PER_TEST_SP,
                                                      NUM_DOF_PER_TRIAL_SP >
{
public:
  /// Alias for the base class. (i.e. #geos::finiteElement::ImplicitKernelBase)
  using Base = ImplicitKernelBase< SUBREGION_TYPE,
                                   CONSTITUTIVE_TYPE,
                                   FE_TYPE,
                                   NUM_DOF_PER_TEST_SP,
                                   NUM_DOF_PER_TRIAL_SP >;


  using typename Base::StackVariables;
  using Base::m_dofRankOffset;
  using Base::m_dt;

  using Base::setup;

  /**
   * @brief Constructor
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param targetRegionIndex Index of the region the subregion belongs to.
   * @param inputDofNumber The dof number for the primary field.
   * @param rankOffset dof index offset of current rank
   * @param inputSparsity The sparsity pattern to fill.
   * @param inputDt The timestep for the physics update.
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   */
  SparsityKernelBase( NodeManager const & nodeManager,
                      EdgeManager const & edgeManager,
                      FaceManager const & faceManager,
                      localIndex const targetRegionIndex,
                      SUBREGION_TYPE const & elementSubRegion,
                      FE_TYPE const & finiteElementSpace,
                      CONSTITUTIVE_TYPE & inputConstitutiveType,
                      arrayView1d< globalIndex const > const & inputDofNumber,
                      globalIndex const rankOffset,
                      real64 const inputDt,
                      SparsityPattern< globalIndex > & inputSparsity ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          CRSMatrixView< real64, globalIndex const >(),
          arrayView1d< real64 >(),
          inputDt ),
    m_sparsity( inputSparsity )
  {}


  /**
   * @copydoc geos::finiteElement::KernelBase::complete
   *
   * In this implementation, only the matrix values are inserted, making this
   * implementation appropriate for generating the sparsity pattern.
   */
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );

    for( localIndex r=0; r<stack.numRows; ++r )
    {
      localIndex const row = stack.localRowDofIndex[r] - m_dofRankOffset;
      if( row < 0 || row >= m_sparsity.numRows() ) continue;
      for( localIndex c=0; c<stack.numCols; ++c )
      {
        m_sparsity.insertNonZero( row,
                                  stack.localColDofIndex[c] );
      }
    }
    return 0;
  }



  /**
   * @brief Kernel Launcher.
   * @tparam POLICY The RAJA policy to use for the launch.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam KERNEL_TYPE The type of Kernel to execute.
   * @param numElems The number of elements to process in this launch.
   * @param kernelComponent The instantiation of KERNEL_TYPE to execute.
   * @return The maximum residual.
   *
   * This is a generic launching function for all of the finite element kernels
   * that follow the interface set by KernelBase.
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    // launch the kernel
    forAll< POLICY >( numElems,
                      [=] ( localIndex const k )
    {
      // allocate the stack variables
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );

      kernelComponent.complete( k, stack );

    } );
    return 0;
  }
private:
  SparsityPattern< globalIndex > & m_sparsity;
};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @brief Helper struct to define a specialization of
 *   #::geos::finiteElement::SparsityKernelBase that may be used to generate the sparsity pattern.
 * @tparam KERNEL_TEMPLATE Templated class that defines the physics kernel.
 *   Most likely derives from SparsityKernelBase.
 */
template< template< typename,
                    typename,
                    typename > class KERNEL_TEMPLATE >
class SparsityKernelFactory
{
public:

  /**
   * @brief Constructor.
   * @param inputDofNumber An array containing the input degree of freedom numbers.
   * @param rankOffset The global rank offset.
   * @param inputSparsityPattern The local sparsity pattern.
   */
  SparsityKernelFactory( arrayView1d< globalIndex const > const & inputDofNumber,
                         globalIndex const rankOffset,
                         SparsityPattern< globalIndex > & inputSparsityPattern ):
    m_inputDofNumber( inputDofNumber ),
    m_rankOffset( rankOffset ),
    m_inputSparsityPattern( inputSparsityPattern )
  {}

  /**
   * @brief Return a new instance of @c SparsityKernelBase specialized for @c KERNEL_TEMPLATE.
   * @tparam SUBREGION_TYPE The type of of @p elementSubRegion.
   * @tparam CONSTITUTIVE_TYPE The type of @p inputConstitutiveType.
   * @tparam FE_TYPE The type of @p finiteElementSpace.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param targetRegionIndex The target region index.
   * @param elementSubRegion The sub region on which to generate the sparsity.
   * @param finiteElementSpace The finite element space.
   * @param inputConstitutiveType The constitutive relation.
   * @return A new instance of @c SparsityKernelBase specialized for @c KERNEL_TEMPLATE.
   */
  template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
  auto createKernel( NodeManager const & nodeManager,
                     EdgeManager const & edgeManager,
                     FaceManager const & faceManager,
                     localIndex const targetRegionIndex,
                     SUBREGION_TYPE const & elementSubRegion,
                     FE_TYPE const & finiteElementSpace,
                     CONSTITUTIVE_TYPE & inputConstitutiveType )
  {
    using Kernel = KERNEL_TEMPLATE< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >;

    return SparsityKernelBase< SUBREGION_TYPE,
                               CONSTITUTIVE_TYPE,
                               FE_TYPE,
                               Kernel::numDofPerTestSupportPoint,
                               Kernel::numDofPerTrialSupportPoint >( nodeManager,
                                                                     edgeManager,
                                                                     faceManager,
                                                                     targetRegionIndex,
                                                                     elementSubRegion,
                                                                     finiteElementSpace,
                                                                     inputConstitutiveType,
                                                                     m_inputDofNumber,
                                                                     m_rankOffset,
                                                                     0.0, //dt but not needed
                                                                     m_inputSparsityPattern );
  }

private:
  /// The input degree of freedom numbers.
  arrayView1d< globalIndex const > const & m_inputDofNumber;
  /// The global rank offset.
  globalIndex const m_rankOffset;
  /// The local sparsity pattern.
  SparsityPattern< globalIndex > & m_inputSparsityPattern;
};

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @brief Fills matrix sparsity.
 * @tparam REGION_TYPE The type of region to loop over. TODO make this a
 *                     parameter pack?
 * @tparam KERNEL_TEMPLATE The type of template for the physics kernel, which
 *                         conforms to the interface specified by KernelBase.
 * @param mesh The MeshLevel object.
 * @param targetRegions The names of the target regions(of type @p REGION_TYPE)
 *                      to apply the @p KERNEL_TEMPLATE.
 * @param discretizationName The name of the finite element discretization.
 * @param inputDofNumber The global degree of freedom numbers.
 * @param rankOffset Offset of dof indices on curren rank.
 * @param inputSparsityPattern The local sparsity pattern to fill.
 * @return 0
 *
 * Fills matrix sparsity using information from physics specific implementation
 * of #geos::finiteElement::KernelBase interface using the
 * #geos::finiteElement::ImplicitKernelBase alias to specialize
 * #geos::finiteElement::SparsityKernelBase to conform with the template
 * pattern specified in the physics kernels.
 */
template< typename REGION_TYPE,
          template< typename SUBREGION_TYPE,
                    typename CONSTITUTIVE_TYPE,
                    typename FE_TYPE > class KERNEL_TEMPLATE >
static
real64 fillSparsity( MeshLevel & mesh,
                     arrayView1d< string const > const & targetRegions,
                     string const & discretizationName,
                     arrayView1d< globalIndex const > const & inputDofNumber,
                     globalIndex const rankOffset,
                     SparsityPattern< globalIndex > & inputSparsityPattern )
{
  GEOS_MARK_FUNCTION;

  SparsityKernelFactory< KERNEL_TEMPLATE > KernelFactory( inputDofNumber, rankOffset, inputSparsityPattern );

  regionBasedKernelApplication< serialPolicy,
                                constitutive::NullModel,
                                REGION_TYPE >( mesh,
                                               targetRegions,
                                               discretizationName,
                                               string(),
                                               KernelFactory );

  return 0;
}

}
}



#endif /* GEOS_FINITEELEMENT_SPARSITYKERNELBASE_HPP_ */
