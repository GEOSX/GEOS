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
 * @file KernelBase.hpp
 */

#ifndef GEOS_FINITEELEMENT_KERNELBASE_HPP_
#define GEOS_FINITEELEMENT_KERNELBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/MeshLevel.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

/**
 * @brief This macro allows solvers to select a subset of FE_TYPES on which the dispatch is done. If none are selected, by default all the
 * FE_TYPES apply.
 */
#ifndef SELECTED_FE_TYPES
#define SELECTED_FE_TYPES BASE_FE_TYPES
#endif

namespace geos
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

/**
 * @class KernelBase
 * @brief Define the base interface for finite element kernels.
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
 * @tparam CONSTITUTIVE_TYPE The type of constitutive model present in the
 *                           subregion.
 * @tparam NUM_TEST_SUPPORT_POINTS_PER_ELEM The number of test space support
 *                                          points per element.
 * @tparam NUM_TRIAL_SUPPORT_POINTS_PER_ELEM The number of trial space support
 *                                           points per element.
 * @tparam NUM_DOF_PER_TEST_SP The number of DOF per test support point.
 * @tparam NUM_DOF_PER_TRIAL_SP The number of DOF per trial support point.
 *
 * ### General KernelBase Description
 *
 * KernelBase defines an interface for implementing finite element kernels
 * that will be callable by the family of kernel launch functions. Specific
 * physics kernels may or may not derive from KernelBase, but must follow
 * the same interface in order to be callable from the generic launching
 * functions.
 *
 * The template parameters of KernelBase should be duplicated as part of the
 * interface, EXCEPT for @p NUM_DOF_PER_TEST_SP and @p NUM_DOF_PER_TRIAL_SP.
 * These values should be set internally by the physics solver since each
 * physics discretization will have a constant intrinsic value for these
 * quantities. For example, when solving or the heat equation with scalar
 * temperature as the primary variable at the support point, these will have
 * a value of 1. In contrast, when solving a solid mechanics problem, with
 * vector displacement as the primary variable at the support point, these
 * will have a value of 3. Note that the interface provided by
 * geos::finiteElement::RegionBasedKernelApplication will construct a
 * kernel assuming only the first 4 template arguments.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class KernelBase
{
public:
  /// Compile time value for the number of test function support points per
  /// element.
  static constexpr int maxNumTestSupportPointsPerElem  = FE_TYPE::maxSupportPoints;

  /// Compile time value for the number of trial function support points per
  /// element.
  static constexpr int maxNumTrialSupportPointsPerElem = FE_TYPE::maxSupportPoints;

  /// Compile time value for the number of degrees of freedom per test function
  /// support point.
  static constexpr int numDofPerTestSupportPoint    = NUM_DOF_PER_TEST_SP;

  /// Compile time value for the number of degrees of freedom per trial
  /// function support point.
  static constexpr int numDofPerTrialSupportPoint   = NUM_DOF_PER_TRIAL_SP;

  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  /**
   * @brief Constructor
   * @param elementSubRegion Reference to the SUBREGION_TYPE(class template
   *                         parameter) object.
   * @param inputConstitutiveType The constitutive object.
   * @param finiteElementSpace Placeholder for the finite element space object,
   *                           which currently doesn't do much.
   */
  KernelBase( SUBREGION_TYPE const & elementSubRegion,
              FE_TYPE const & finiteElementSpace,
              CONSTITUTIVE_TYPE & inputConstitutiveType ):
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_elemGhostRank( elementSubRegion.ghostRank() ),
    m_constitutiveUpdate( inputConstitutiveType.createKernelUpdates() ),
    m_finiteElementSpace( finiteElementSpace )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables allocated on the stack.
   *
   * ### ImplicitKernelBase::StackVariables Description
   *
   * Contains variables that will be allocated on the stack of the main kernel.
   * This will typically consist of local arrays to hold data mapped from the
   * global data arrays, and/or local storage for the residual and jacobian
   * contributions.
   */
  struct StackVariables
  {};

  /**
   * @brief Performs the setup phase for the kernel.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::setup() Description
   *
   * The operations typically found in setup are thing such as the collection
   * of global data into local stack storage.
   */
  GEOS_HOST_DEVICE
  inline
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( stack );
  }

  /**
   * @brief Performs a state update at a quadrature point.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param q The quadrature point index.
   * @param stack The StackVariable object that hold the stack variables.
   *
   * ### KernelBase::quadraturePointKernel() Description
   *
   * The operations found here are the mapping from the support points to the
   * quadrature point, calculation of gradients, etc. From this data the
   * state of the constitutive model is updated if required by the physics
   * package.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( q );
    GEOS_UNUSED_VAR( stack );
  }

  /**
   * @brief Performs the complete phase for the kernel.
   * @tparam STACK_VARIABLE_TYPE The type of StackVariable that holds the stack
   *                             variables. This is most likely a defined in a
   *                             type that derives from KernelBase.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   * @return The maximum contribution to the residual.
   *
   * ### KernelBase::complete() Description
   *
   * The operations typically found in complete are the mapping of the local
   * Jacobian and Residual into the global Jacobian and Residual.
   */
  GEOS_HOST_DEVICE
  inline
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( k );
    GEOS_UNUSED_VAR( stack );
    return 0;
  }


  /**
   * @brief Kernel Launcher.
   * @tparam POLICY The RAJA policy to use for the launch.
   * @tparam NUM_QUADRATURE_POINTS The number of quadrature points per element.
   * @tparam KERNEL_TYPE The type of Kernel to execute.
   * @param numElems The number of elements to process in this launch.
   * @param kernelComponent The instantiation of KERNEL_TYPE to execute.
   * @return The maximum residual contribution.
   *
   * This is a generic launching function for all of the finite element kernels
   * that follow the interface set by KernelBase.
   */
  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;

    // Define a RAJA reduction variable to get the maximum residual contribution.
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( numElems,
                      [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      // #pragma unroll
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      maxResidual.max( kernelComponent.complete( k, stack ) );
    } );
    return maxResidual.get();
  }
  //END_kernelLauncher

protected:
  /// The element to nodes map.
  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elemsToNodes;

  /// The element ghost rank array.
  arrayView1d< integer const > const m_elemGhostRank;

  /// The constitutive update object used to update the constitutive state,
  /// and extract constitutive data.
  typename CONSTITUTIVE_TYPE::KernelWrapper const m_constitutiveUpdate;

  /// The finite element space/discretization object for the element type in
  /// the SUBREGION_TYPE.
  FE_TYPE const & m_finiteElementSpace;
};

/**
 * @class KernelFactory
 * @brief Used to forward arguments to a class that implements the KernelBase interface.
 * @tparam KERNEL_TYPE The template class to construct, should implement the KernelBase interface.
 * @tparam ARGS The arguments used to construct a @p KERNEL_TYPE in addition to the standard arguments.
 */
template< template< typename SUBREGION_TYPE,
                    typename CONSTITUTIVE_TYPE,
                    typename FE_TYPE > class KERNEL_TYPE,
          typename ... ARGS >
class KernelFactory
{
public:

  /**
   * @brief Initialize the factory.
   * @param args The arguments used to construct a @p KERNEL_TYPE in addition to the standard arguments.
   */
  KernelFactory( ARGS ... args ):
    m_args( args ... )
  {}

  /**
   * @brief Create a new kernel with the given standard arguments.
   * @tparam SUBREGION_TYPE The type of @p elementSubRegion.
   * @tparam CONSTITUTIVE_TYPE The type of @p inputConstitutiveType.
   * @tparam FE_TYPE The type of @p finiteElementSpace.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param targetRegionIndex The target region index.
   * @param elementSubRegion The subregion to execute on.
   * @param finiteElementSpace The finite element space.
   * @param inputConstitutiveType The constitutive relation.
   * @return A new kernel constructed with the given arguments and @c ARGS.
   */
  template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
  KERNEL_TYPE< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE > createKernel(
    NodeManager & nodeManager,
    EdgeManager const & edgeManager,
    FaceManager const & faceManager,
    localIndex const targetRegionIndex,
    SUBREGION_TYPE const & elementSubRegion,
    FE_TYPE const & finiteElementSpace,
    CONSTITUTIVE_TYPE & inputConstitutiveType )
  {
    camp::tuple< NodeManager &,
                 EdgeManager const &,
                 FaceManager const &,
                 localIndex const,
                 SUBREGION_TYPE const &,
                 FE_TYPE const &,
                 CONSTITUTIVE_TYPE & > standardArgs { nodeManager,
                                                      edgeManager,
                                                      faceManager,
                                                      targetRegionIndex,
                                                      elementSubRegion,
                                                      finiteElementSpace,
                                                      inputConstitutiveType };

    auto allArgs = camp::tuple_cat_pair( standardArgs, m_args );
    return camp::make_from_tuple< KERNEL_TYPE< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE > >( allArgs );
  }

private:
  /// The arguments to append to the standard kernel constructor arguments.
  camp::tuple< ARGS ... > m_args;
};


//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//START_regionBasedKernelApplication
/**
 * @brief Performs a loop over specific regions (by type and name) and calls a kernel launch on the subregions
 *   with compile time knowledge of sub-loop bounds such as number of nodes and quadrature points per element.
 * @tparam POLICY The RAJA launch policy to pass to the kernel launch.
 * @tparam CONSTITUTIVE_BASE The common base class for constitutive pass-thru/dispatch which gives the kernel
 *   launch compile time knowledge of the constitutive model. This is achieved through a call to the
 *   ConstitutivePassThru function which should have a specialization for CONSTITUTIVE_BASE implemented in
 *   order to perform the compile time dispatch.
 * @tparam SUBREGION_TYPE The type of subregion to loop over. TODO make this a parameter pack?
 * @tparam KERNEL_FACTORY The type of @p kernelFactory, typically an instantiation of @c KernelFactory, and
 *   must adhere to that interface.
 * @param mesh The MeshLevel object.
 * @param targetRegions The names of the target regions(of type @p SUBREGION_TYPE) to apply the @p KERNEL_TEMPLATE.
 * @param finiteElementName The name of the finite element.
 * @param constitutiveStringName The key to the constitutive model name found on the Region.
 * @param kernelFactory The object used to construct the kernel.
 * @return The maximum contribution to the residual, which may be used to scale the residual.
 *
 * @details Loops over all regions Applies/Launches a kernel specified by the @p KERNEL_TEMPLATE through
 * #::geos::finiteElement::KernelBase::kernelLaunch().
 */
template< typename POLICY,
          typename CONSTITUTIVE_BASE,
          typename SUBREGION_TYPE,
          typename KERNEL_FACTORY >
static
real64 regionBasedKernelApplication( MeshLevel & mesh,
                                     arrayView1d< string const > const & targetRegions,
                                     string const & finiteElementName,
                                     string const & constitutiveStringName,
                                     KERNEL_FACTORY & kernelFactory )
{
  GEOS_MARK_FUNCTION;
  // save the maximum residual contribution for scaling residuals for convergence criteria.
  real64 maxResidualContribution = 0;

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elementRegionManager = mesh.getElemManager();

  // Loop over all sub-regions in regions of type SUBREGION_TYPE, that are listed in the targetRegions array.
  elementRegionManager.forElementSubRegions< SUBREGION_TYPE >( targetRegions,
                                                               [&constitutiveStringName,
                                                                &maxResidualContribution,
                                                                &nodeManager,
                                                                &edgeManager,
                                                                &faceManager,
                                                                &kernelFactory,
                                                                &finiteElementName]
                                                                 ( localIndex const targetRegionIndex, auto & elementSubRegion )
  {
    localIndex const numElems = elementSubRegion.size();

    // Get the constitutive model...and allocate a null constitutive model if required.

    constitutive::ConstitutiveBase * constitutiveRelation = nullptr;
    constitutive::NullModel * nullConstitutiveModel = nullptr;
    if( elementSubRegion.template hasWrapper< string >( constitutiveStringName ) )
    {
      string const & constitutiveName = elementSubRegion.template getReference< string >( constitutiveStringName );
      constitutiveRelation = &elementSubRegion.template getConstitutiveModel( constitutiveName );
    }
    else
    {
      nullConstitutiveModel = &elementSubRegion.template registerGroup< constitutive::NullModel >( "nullModelGroup" );
      constitutiveRelation = nullConstitutiveModel;
    }

    // Call the constitutive dispatch which converts the type of constitutive model into a compile time constant.
    constitutive::ConstitutivePassThru< CONSTITUTIVE_BASE >::execute( *constitutiveRelation,
                                                                      [&maxResidualContribution,
                                                                       &nodeManager,
                                                                       &edgeManager,
                                                                       &faceManager,
                                                                       targetRegionIndex,
                                                                       &kernelFactory,
                                                                       &elementSubRegion,
                                                                       &finiteElementName,
                                                                       numElems]
                                                                        ( auto & castedConstitutiveRelation )
    {
      FiniteElementBase &
      subRegionFE = elementSubRegion.template getReference< FiniteElementBase >( finiteElementName );

      finiteElement::FiniteElementDispatchHandler< SELECTED_FE_TYPES >::dispatch3D( subRegionFE,
                                                                                    [&maxResidualContribution,
                                                                                     &nodeManager,
                                                                                     &edgeManager,
                                                                                     &faceManager,
                                                                                     targetRegionIndex,
                                                                                     &kernelFactory,
                                                                                     &elementSubRegion,
                                                                                     numElems,
                                                                                     &castedConstitutiveRelation] ( auto const finiteElement )
      {
        auto kernel = kernelFactory.createKernel( nodeManager,
                                                  edgeManager,
                                                  faceManager,
                                                  targetRegionIndex,
                                                  elementSubRegion,
                                                  finiteElement,
                                                  castedConstitutiveRelation );

        using KERNEL_TYPE = decltype( kernel );

        // Call the kernelLaunch function, and store the maximum contribution to the residual.
        maxResidualContribution =
          std::max( maxResidualContribution,
                    KERNEL_TYPE::template kernelLaunch< POLICY, KERNEL_TYPE >( numElems, kernel ) );
      } );
    } );

    // Remove the null constitutive model (not required, but cleaner)
    if( nullConstitutiveModel )
    {
      elementSubRegion.deregisterGroup( "nullModelGroup" );
    }

  } );

  return maxResidualContribution;
}
//END_regionBasedKernelApplication

} // namespace finiteElement
} // namespace geos



#endif /* GEOS_FINITEELEMENT_KERNELBASE_HPP_ */
