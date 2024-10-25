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
 * @file StrainHelper.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_STRAINHELPER_HPP_
#define GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_STRAINHELPER_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{
/**
 * @class AverageStrainOverQuadraturePoints
 * @tparam SUBREGION_TYPE the subRegion type
 * @tparam FE_TYPE the finite element type
 */
template< typename SUBREGION_TYPE,
          typename FE_TYPE >
class AverageStrainOverQuadraturePoints :
  public AverageOverQuadraturePointsBase< SUBREGION_TYPE,
                                          FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = AverageOverQuadraturePointsBase< SUBREGION_TYPE,
                                                FE_TYPE >;

  using Base::m_elementVolume;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;

  /**
   * @brief Constructor for the class
   * @param nodeManager the node manager
   * @param edgeManager the edge manager
   * @param faceManager the face manager
   * @param elementSubRegion the element subRegion
   * @param finiteElementSpace the finite element space
   * @param displacement the displacement solution field
   * @param avgStrain the strain averaged over quadrature points
   */
  AverageStrainOverQuadraturePoints( NodeManager & nodeManager,
                                     EdgeManager const & edgeManager,
                                     FaceManager const & faceManager,
                                     SUBREGION_TYPE const & elementSubRegion,
                                     FE_TYPE const & finiteElementSpace,
                                     fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const displacement,
                                     fields::solidMechanics::arrayView2dLayoutStrain const avgStrain ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace ),
    m_displacement( displacement ),
    m_avgStrain( avgStrain )
  {}

  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : Base::StackVariables
  {real64 uLocal[FE_TYPE::maxSupportPoints][3]; };

  /**
   * @brief Performs the setup phase for the kernel.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    Base::setup( k, stack );

    for( localIndex a = 0; a < FE_TYPE::maxSupportPoints; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );
      for( int i = 0; i < 3; ++i )
      {
        stack.uLocal[a][i] = m_displacement[localNodeIndex][i];
      }
    }

    for( int icomp = 0; icomp < 6; ++icomp )
    {
      m_avgStrain[k][icomp] = 0.0;
    }
  }

  /**
   * @brief Increment the average property with the contribution of the property at this quadrature point
   * @param k The element index
   * @param q The quadrature point index
   * @param stack The StackVariables object that hold the stack variables.
   */
  GEOS_HOST_DEVICE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    //real64 const weight = FE_TYPE::transformedQuadratureWeight( q, stack.xLocal, stack.feStack ) / m_elementVolume[k];

    real64 dNdX[ FE_TYPE::maxSupportPoints ][3];
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, stack.feStack, dNdX );
    real64 strain[6] = {0.0};
    FE_TYPE::symmetricGradient( dNdX, stack.uLocal, strain );

    for( int icomp = 0; icomp < 6; ++icomp )
    {
      m_avgStrain[k][icomp] += detJxW*strain[icomp]/m_elementVolume[k];
    }
  }

  /**
   * @brief Launch the kernel over the elements in the subRegion
   * @tparam POLICY the kernel policy
   * @tparam KERNEL_TYPE the type of kernel
   * @param numElems the number of elements in the subRegion
   * @param kernelComponent the kernel component
   */
  template< typename POLICY,
            typename KERNEL_TYPE >
  static void
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    forAll< POLICY >( numElems,
                      [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q = 0; q < FE_TYPE::numQuadraturePoints; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
    } );
  }

protected:

  /// The displacement solution
  fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const m_displacement;

  /// The average strain
  fields::solidMechanics::arrayView2dLayoutStrain const m_avgStrain;

};



/**
 * @class AverageStrainOverQuadraturePointsKernelFactory
 * @brief Class to create and launch the kernel
 */
class AverageStrainOverQuadraturePointsKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam SUBREGION_TYPE the subRegion type
   * @tparam FE_TYPE the finite element type
   * @tparam POLICY the kernel policy
   * @param nodeManager the node manager
   * @param edgeManager the edge manager
   * @param faceManager the face manager
   * @param elementSubRegion the element subRegion
   * @param finiteElementSpace the finite element space
   * @param property the property at quadrature points
   * @param averageProperty the property averaged over quadrature points
   */
  template< typename SUBREGION_TYPE,
            typename FE_TYPE,
            typename POLICY >
  static void
  createAndLaunch( NodeManager & nodeManager,
                   EdgeManager const & edgeManager,
                   FaceManager const & faceManager,
                   SUBREGION_TYPE const & elementSubRegion,
                   FE_TYPE const & finiteElementSpace,
                   fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const displacement,
                   fields::solidMechanics::arrayView2dLayoutStrain const avgStrain )
  {
    AverageStrainOverQuadraturePoints< SUBREGION_TYPE, FE_TYPE >
    kernel( nodeManager, edgeManager, faceManager, elementSubRegion, finiteElementSpace,
            displacement, avgStrain );

    AverageStrainOverQuadraturePoints< SUBREGION_TYPE, FE_TYPE >::template
    kernelLaunch< POLICY >( elementSubRegion.size(), kernel );
  }
};



}


#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_STRAINHELPER_HPP_ */
