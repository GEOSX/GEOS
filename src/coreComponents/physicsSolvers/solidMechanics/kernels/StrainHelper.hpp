/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
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
#include "constitutive/ConstitutivePassThru.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{
/**
 * @class AverageStrainOverQuadraturePoints
 * @tparam SUBREGION_TYPE the subRegion type
 * @tparam FE_TYPE the finite element type
 * @tparam SOLID_TYPE the solid mechanics constitutuve type
 */
template< typename FE_TYPE,
          typename SOLID_TYPE >
class AverageStrainOverQuadraturePoints :
  public AverageOverQuadraturePointsBase< CellElementSubRegion,
                                          FE_TYPE >
{
public:

  /// Alias for the base class;
  using Base = AverageOverQuadraturePointsBase< CellElementSubRegion,
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
                                     CellElementSubRegion const & elementSubRegion,
                                     FE_TYPE const & finiteElementSpace,
                                     SOLID_TYPE const & solidModel,
                                     fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const displacement,
                                     fields::solidMechanics::arrayViewConst2dLayoutIncrDisplacement const displacementInc,
                                     fields::solidMechanics::arrayView2dLayoutStrain const avgStrain,
                                     fields::solidMechanics::arrayView2dLayoutStrain const avgPlasticStrain ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace ),
    m_solidUpdate( solidModel.createKernelUpdates()),
    m_displacement( displacement ),
    m_displacementInc( displacementInc ),
    m_avgStrain( avgStrain ),
    m_avgPlasticStrain( avgPlasticStrain )
  {}

  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : Base::StackVariables
  {real64 uLocal[FE_TYPE::maxSupportPoints][3];
   real64 uHatLocal[FE_TYPE::maxSupportPoints][3]; };

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
        stack.uHatLocal[a][i] = m_displacementInc[localNodeIndex][i];
      }
    }

    for( int icomp = 0; icomp < 6; ++icomp )
    {
      m_avgStrain[k][icomp] = 0.0;
      //m_avgPlasticStrain[k][icomp] = 0.0;
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
    real64 strainInc[6] = {0.0};
    FE_TYPE::symmetricGradient( dNdX, stack.uLocal, strain );
    FE_TYPE::symmetricGradient( dNdX, stack.uHatLocal, strainInc );

    real64 elasticStrainInc[6] = {0.0};
    m_solidUpdate.getElasticStrainInc( k, q, elasticStrainInc );

    real64 conversionFactor[6] = {1.0, 1.0, 1.0, 0.5, 0.5, 0.5}; // used for converting from engineering shear to tensor shear

    for( int icomp = 0; icomp < 6; ++icomp )
    {
      m_avgStrain[k][icomp] += conversionFactor[icomp]*detJxW*strain[icomp]/m_elementVolume[k];

      // This is a hack to handle boundary conditions such as those seen in plane-strain wellbore problems
      // Essentially, if bcs are constraining the strain (and thus total displacement), we do not accumulate any plastic strain (regardless
      // of stresses in material law)
      if( std::abs( strainInc[icomp] ) > 1.0e-8 )
      {
        m_avgPlasticStrain[k][icomp] += conversionFactor[icomp]*detJxW*(strainInc[icomp] - elasticStrainInc[icomp])/m_elementVolume[k];
      }
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

  /// The material
  typename SOLID_TYPE::KernelWrapper const m_solidUpdate;

  /// The displacement solution
  fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const m_displacement;

  /// The displacement increment
  fields::solidMechanics::arrayViewConst2dLayoutIncrDisplacement const m_displacementInc;

  /// The average strain
  fields::solidMechanics::arrayView2dLayoutStrain const m_avgStrain;

  /// The average plastic strain
  fields::solidMechanics::arrayView2dLayoutStrain const m_avgPlasticStrain;

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
   * @tparam SOLID_TYPE the constitutive type
   * @tparam POLICY the kernel policy
   * @param nodeManager the node manager
   * @param edgeManager the edge manager
   * @param faceManager the face manager
   * @param elementSubRegion the element subRegion
   * @param finiteElementSpace the finite element space
   * @param property the property at quadrature points
   * @param averageProperty the property averaged over quadrature points
   */
  template< typename FE_TYPE,
            typename SOLID_TYPE,
            typename POLICY >
  static void
  createAndLaunch( NodeManager & nodeManager,
                   EdgeManager const & edgeManager,
                   FaceManager const & faceManager,
                   CellElementSubRegion const & elementSubRegion,
                   FE_TYPE const & finiteElementSpace,
                   SOLID_TYPE const & solidModel,
                   fields::solidMechanics::arrayViewConst2dLayoutTotalDisplacement const displacement,
                   fields::solidMechanics::arrayViewConst2dLayoutIncrDisplacement const displacementInc,
                   fields::solidMechanics::arrayView2dLayoutStrain const avgStrain,
                   fields::solidMechanics::arrayView2dLayoutStrain const avgPlasticStrain )
  {
    AverageStrainOverQuadraturePoints< FE_TYPE, SOLID_TYPE >
    kernel( nodeManager, edgeManager, faceManager, elementSubRegion, finiteElementSpace,
            solidModel, displacement, displacementInc, avgStrain, avgPlasticStrain );

    AverageStrainOverQuadraturePoints< FE_TYPE, SOLID_TYPE >::template
    kernelLaunch< POLICY >( elementSubRegion.size(), kernel );
  }
};



}


#endif /* GEOS_PHYSICSSOLVERS_SOLIDMECHANICS_KERNELS_STRAINHELPER_HPP_ */
