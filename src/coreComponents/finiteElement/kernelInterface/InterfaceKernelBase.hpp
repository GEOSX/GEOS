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
 * @file InterfaceKernelBase.hpp
 */

#ifndef GEOS_FINITEELEMENT_INTERFACEKERNELBASE_HPP_
#define GEOS_FINITEELEMENT_INTERFACEKERNELBASE_HPP_

#include "ImplicitKernelBase.hpp"

/**
 * @brief This macro allows solvers to select a subset of FE_TYPES_2D on which the dispatch is done.
 * If none are selected, by default all the BASE_FE_TYPES_2D apply.
 */
#ifndef SELECTED_FE_TYPES_2D
#define SELECTED_FE_TYPES_2D BASE_FE_TYPES_2D
#endif

namespace geos
{

/**
 * @namespace finiteElement Contains the finite element implementation.
 */
namespace finiteElement
{

/**
 * @class InterfaceKernelBase
 * @brief Define the base class for interface finite element kernels.
 *   (2D finite elements belong to FaceElementSubRegion).
 * @tparam CONSTITUTIVE_TYPE The type of constitutive model present in the
 *                           FaceElementSubRegion.
 * @tparam FE_TYPE The type of finite element.
 * @tparam NUM_DOF_PER_TEST_SP The number of DOF per test support point.
 * @tparam NUM_DOF_PER_TRIAL_SP The number of DOF per trial support point.
 *
 */

template< typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class InterfaceKernelBase : public ImplicitKernelBase< FaceElementSubRegion,
                                                       CONSTITUTIVE_TYPE,
                                                       FE_TYPE,
                                                       NUM_DOF_PER_TEST_SP,
                                                       NUM_DOF_PER_TRIAL_SP >
{
public:

  /// Alias for the base class, i.e., geos::finiteElement::ImplicitKernelBase)
  using Base = ImplicitKernelBase< FaceElementSubRegion,
                                   CONSTITUTIVE_TYPE,
                                   FE_TYPE,
                                   NUM_DOF_PER_TEST_SP,
                                   NUM_DOF_PER_TRIAL_SP >;

  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   */
  InterfaceKernelBase( NodeManager const & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       localIndex const targetRegionIndex,
                       FaceElementSubRegion & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE & inputConstitutiveType,
                       arrayView1d< globalIndex const > const inputDofNumber,
                       globalIndex const rankOffset,
                       CRSMatrixView< real64, globalIndex const > const inputMatrix,
                       arrayView1d< real64 > const inputRhs,
                       real64 const inputDt ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt )
  {}

  //***************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables
  {};

};

/**
 * @class InterfaceKernelFactory
 * @brief Used to forward arguments to a class that implements the InterfaceKernelBase interface.
 * @tparam KERNEL_TYPE The template class to construct, should implement the InterfaceKernelBase interface.
 * @tparam ARGS The arguments used to construct a @p KERNEL_TYPE in addition to the standard arguments.
 */
template< template< typename CONSTITUTIVE_TYPE,
                    typename FE_TYPE > class KERNEL_TYPE,
          typename ... ARGS >
class InterfaceKernelFactory
{
public:

  /**
   * @brief Initialize the factory.
   * @param args The arguments used to construct a @p KERNEL_TYPE in addition to the standard arguments.
   */
  InterfaceKernelFactory( ARGS ... args ):
    m_args( args ... )
  {}

  /**
   * @brief Create a new kernel with the given standard arguments.
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
  template< typename CONSTITUTIVE_TYPE, typename FE_TYPE >
  KERNEL_TYPE< CONSTITUTIVE_TYPE, FE_TYPE > createKernel(
    NodeManager & nodeManager,
    EdgeManager const & edgeManager,
    FaceManager const & faceManager,
    localIndex const targetRegionIndex,
    FaceElementSubRegion & elementSubRegion,
    FE_TYPE const & finiteElementSpace,
    CONSTITUTIVE_TYPE & inputConstitutiveType )
  {
    camp::tuple< NodeManager &,
                 EdgeManager const &,
                 FaceManager const &,
                 localIndex const,
                 FaceElementSubRegion &,
                 FE_TYPE const &,
                 CONSTITUTIVE_TYPE & > standardArgs { nodeManager,
                                                      edgeManager,
                                                      faceManager,
                                                      targetRegionIndex,
                                                      elementSubRegion,
                                                      finiteElementSpace,
                                                      inputConstitutiveType };

    auto allArgs = camp::tuple_cat_pair( standardArgs, m_args );
    return camp::make_from_tuple< KERNEL_TYPE< CONSTITUTIVE_TYPE, FE_TYPE > >( allArgs );

  }

private:
  /// The arguments to append to the standard kernel constructor arguments.
  camp::tuple< ARGS ... > m_args;
};

//*****************************************************************************

//START_interfaceBasedKernelApplication
/**
 * @brief Performs a loop over FaceElementSubRegion and calls a kernel launch
 *   with compile time knowledge of sub-loop bounds such as number of nodes and
 *   quadrature points per element.
 * @tparam POLICY The RAJA launch policy to pass to the kernel launch.
 * @tparam CONSTITUTIVE_BASE The common base class for constitutive pass-thru/dispatch which gives the kernel
 *   launch compile time knowledge of the constitutive model.
 * @tparam KERNEL_FACTORY The type of @p interfaceKernelFactory, typically an instantiation of @c InterfaceKernelFactory, and
 *   must adhere to that interface.
 * @param mesh The MeshLevel object.
 * @param targetRegionName The names of the target regions to apply the @p KERNEL_TEMPLATE.
 * @param faceElementList List of element of the same type belongs to FaceElementSubRegion.
 * @param subRegionFE Finite element object.
 * @param constitutiveStringName Key string used to retrieve the constitutive model.
 * @param interfaceKernelFactory The object used to construct the kernel.
 * @return The maximum contribution to the residual, which may be used to scale the residual.
 *
 * @details Loops over all regions Applies/Launches a kernel specified by the @p KERNEL_TEMPLATE through
 * #::geos::finiteElement::KernelBase::kernelLaunch().
 */
template< typename POLICY,
          typename CONSTITUTIVE_BASE,
          typename KERNEL_FACTORY >
static
real64 interfaceBasedKernelApplication( MeshLevel & mesh,
                                        string const & targetRegionName,
                                        arrayView1d< localIndex const > const & faceElementList,
                                        FiniteElementBase const & subRegionFE,
                                        string const & constitutiveStringName,
                                        KERNEL_FACTORY & interfaceKernelFactory )
{
  GEOS_MARK_FUNCTION;


  // save the maximum residual contribution for scaling residuals for convergence criteria.
  real64 maxResidualContribution = 0;

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elementManager = mesh.getElemManager();

  SurfaceElementRegion & region = elementManager.getRegion< SurfaceElementRegion >( targetRegionName );
  FaceElementSubRegion & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
  localIndex const targetRegionIndex = 0;

  // Get the constitutive model...and allocate a null constitutive model if required.
  constitutive::ConstitutiveBase * constitutiveRelation = nullptr;
  constitutive::NullModel * nullConstitutiveModel = nullptr;
  if( subRegion.template hasWrapper< string >( constitutiveStringName ) )
  {
    string const & constitutiveName = subRegion.template getReference< string >( constitutiveStringName );
    constitutiveRelation = &subRegion.template getConstitutiveModel( constitutiveName );
  }
  else
  {
    nullConstitutiveModel = &subRegion.template registerGroup< constitutive::NullModel >( "nullModelGroup" );
    constitutiveRelation = nullConstitutiveModel;
  }

  localIndex const numElems = faceElementList.size();

  // Call the constitutive dispatch which converts the type of constitutive model into a compile time constant.
  constitutive::ConstitutivePassThru< CONSTITUTIVE_BASE >::execute( *constitutiveRelation,
                                                                    [&maxResidualContribution,
                                                                     &nodeManager,
                                                                     &edgeManager,
                                                                     &faceManager,
                                                                     targetRegionIndex,
                                                                     &interfaceKernelFactory,
                                                                     &subRegion,
                                                                     &subRegionFE,
                                                                     numElems]
                                                                      ( auto & castedConstitutiveRelation )
  {

    finiteElement::FiniteElementDispatchHandler< SELECTED_FE_TYPES_2D >::dispatch2D( subRegionFE,
                                                                                     [&maxResidualContribution,
                                                                                      &nodeManager,
                                                                                      &edgeManager,
                                                                                      &faceManager,
                                                                                      targetRegionIndex,
                                                                                      &interfaceKernelFactory,
                                                                                      &subRegion,
                                                                                      numElems,
                                                                                      &castedConstitutiveRelation] ( auto const finiteElement )

    {
      auto kernel = interfaceKernelFactory.createKernel( nodeManager,
                                                         edgeManager,
                                                         faceManager,
                                                         targetRegionIndex,
                                                         subRegion,
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
    subRegion.deregisterGroup( "nullModelGroup" );
  }

  return maxResidualContribution;
}
//END_interfaceBasedKernelApplication

} // namespace finiteElement
} // namespace geos

#endif /* GEOS_FINITEELEMENT_INTERFACEKERNELBASE_HPP_ */
