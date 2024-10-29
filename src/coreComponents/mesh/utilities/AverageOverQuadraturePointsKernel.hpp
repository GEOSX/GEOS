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

/**
 * @file AverageOverQuadraturePointsKernel.hpp
 */

#ifndef GEOS_MESH_UTILITIES_AVERAGEOVERQUADRATUREPOINTSKERNEL_HPP_
#define GEOS_MESH_UTILITIES_AVERAGEOVERQUADRATUREPOINTSKERNEL_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/FiniteElementDispatch.hpp"
#include "mesh/CellElementSubRegion.hpp"

namespace geos
{

/**
 * @class AverageOverQuadraturePointsBase
 * @tparam SUBREGION_TYPE the subRegion type
 * @tparam FE_TYPE the finite element type
 */
template< typename SUBREGION_TYPE,
          typename FE_TYPE >
class AverageOverQuadraturePointsBase
{
public:

  /**
   * @brief Constructor for the class
   * @param nodeManager the node manager
   * @param edgeManager the edge manager
   * @param faceManager the face manager
   * @param elementSubRegion the element subRegion
   * @param finiteElementSpace the finite element space
   */
  AverageOverQuadraturePointsBase( NodeManager & nodeManager,
                                   EdgeManager const & edgeManager,
                                   FaceManager const & faceManager,
                                   SUBREGION_TYPE const & elementSubRegion,
                                   FE_TYPE const & finiteElementSpace ):
    m_finiteElementSpace( finiteElementSpace ),
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_X( nodeManager.referencePosition() ),
    m_elementVolume( elementSubRegion.getElementVolume() )
  {
    finiteElement::FiniteElementBase::
      initialize< FE_TYPE >( nodeManager,
                             edgeManager,
                             faceManager,
                             elementSubRegion,
                             m_meshData );
  }

  //*****************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables
  {
public:

    /**
     * Default constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      xLocal()
    {}

    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[ FE_TYPE::maxSupportPoints ][ 3 ];

    /// Stack variables needed for the underlying FEM type
    typename FE_TYPE::StackVariables feStack;
  };
  //***************************************************************************

  /**
   * @brief Performs the setup phase for the kernel.
   * @param k The element index.
   * @param stack The StackVariable object that hold the stack variables.
   */
  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );

    for( localIndex a = 0; a < FE_TYPE::maxSupportPoints; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( integer i = 0; i < 3; ++i )
      {
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
      }
    }
  }

protected:

  /// The finite element space/discretization object for the element type in
  /// the SUBREGION_TYPE.
  FE_TYPE const & m_finiteElementSpace;

  /// The element to nodes map.
  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elemsToNodes;

  /// The reference position of the nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The volume of the elements
  arrayView1d< real64 const > const m_elementVolume;

  /// Data structure containing mesh data used to setup the finite element
  typename FE_TYPE::template MeshData< SUBREGION_TYPE > m_meshData;
};

/**
 * @class AverageOverQuadraturePoints1D
 * @tparam SUBREGION_TYPE the subRegion type
 * @tparam FE_TYPE the finite element type
 */
template< typename SUBREGION_TYPE,
          typename FE_TYPE >
class AverageOverQuadraturePoints1D :
  public AverageOverQuadraturePointsBase< SUBREGION_TYPE,
                                          FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = AverageOverQuadraturePointsBase< SUBREGION_TYPE,
                                                FE_TYPE >;

  using Base::m_elementVolume;

  /**
   * @brief Constructor for the class
   * @param nodeManager the node manager
   * @param edgeManager the edge manager
   * @param faceManager the face manager
   * @param elementSubRegion the element subRegion
   * @param finiteElementSpace the finite element space
   * @param property the property at quadrature points
   * @param averageProperty the property averaged over quadrature points
   */
  AverageOverQuadraturePoints1D( NodeManager & nodeManager,
                                 EdgeManager const & edgeManager,
                                 FaceManager const & faceManager,
                                 SUBREGION_TYPE const & elementSubRegion,
                                 FE_TYPE const & finiteElementSpace,
                                 arrayView2d< real64 const > const property,
                                 arrayView1d< real64 > const averageProperty ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace ),
    m_property( property ),
    m_averageProperty( averageProperty )
  {}

  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : Base::StackVariables
  {};

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
    m_averageProperty[k] = 0.0;
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
    real64 const weight = FE_TYPE::transformedQuadratureWeight( q, stack.xLocal, stack.feStack ) / m_elementVolume[k];
    m_averageProperty[k] += weight * m_property[k][q];
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

  /// The property living on quadrature points
  arrayView2d< real64 const > const m_property;

  /// The average property
  arrayView1d< real64 > const m_averageProperty;

};


/**
 * @class AverageOverQuadraturePoints1DKernelFactory
 * @brief Class to create and launch the kernel
 */
class AverageOverQuadraturePoints1DKernelFactory
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
                   arrayView2d< real64 const > const property,
                   arrayView1d< real64 > const averageProperty )
  {
    AverageOverQuadraturePoints1D< SUBREGION_TYPE, FE_TYPE >
    kernel( nodeManager, edgeManager, faceManager, elementSubRegion, finiteElementSpace,
            property, averageProperty );

    AverageOverQuadraturePoints1D< SUBREGION_TYPE, FE_TYPE >::template
    kernelLaunch< POLICY >( elementSubRegion.size(), kernel );
  }
};

}

#endif /* GEOS_MESH_UTILITIES_AVERAGEOVERQUADRATUREPOINTSKERNEL_HPP_ */
