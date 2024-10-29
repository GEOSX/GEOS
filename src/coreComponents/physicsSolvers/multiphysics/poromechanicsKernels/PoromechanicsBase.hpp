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
 * @file PoromechanicsBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_POROMECHANICSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_POROMECHANICSBASE_HPP_

#include "finiteElement/BilinearFormUtilities.hpp"
#include "finiteElement/LinearFormUtilities.hpp"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{

namespace poromechanicsKernels
{

/**
 * @brief Defines the kernel structure for solving quasi-static poromechanics.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### Poromechanics Description
 * Defines the KernelBase interface functions required for solving the
 * quasi-static poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class PoromechanicsBase :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_meshData;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param gravityVector The gravity vector.
   */
  PoromechanicsBase( NodeManager const & nodeManager,
                     EdgeManager const & edgeManager,
                     FaceManager const & faceManager,
                     localIndex const targetRegionIndex,
                     SUBREGION_TYPE const & elementSubRegion,
                     FE_TYPE const & finiteElementSpace,
                     CONSTITUTIVE_TYPE & inputConstitutiveType,
                     arrayView1d< globalIndex const > const inputDispDofNumber,
                     globalIndex const rankOffset,
                     CRSMatrixView< real64, globalIndex const > const inputMatrix,
                     arrayView1d< real64 > const inputRhs,
                     real64 const inputDt,
                     real64 const (&gravityVector)[3],
                     string const inputFlowDofKey,
                     string const fluidModelKey ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDispDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputDt ),
    m_X( nodeManager.referencePosition() ),
    m_disp( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
    m_uhat( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
    m_gravityVector{ gravityVector[0], gravityVector[1], gravityVector[2] },
    m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( gravityVector ) ),
    m_solidDensity( inputConstitutiveType.getDensity() ),
    m_flowDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey ) ),
    m_pressure_n( elementSubRegion.template getField< fields::flow::pressure_n >() ),
    m_pressure( elementSubRegion.template getField< fields::flow::pressure >() )
  {
    GEOS_UNUSED_VAR( fluidModelKey );
  }

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geos::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
public:

    static constexpr int numDispDofPerElem =  Base::StackVariables::maxNumRows;

    /// Constructor.
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            u_local(),
            uhat_local(),
            localResidualMomentum( Base::StackVariables::localResidual ),
      dLocalResidualMomentum_dDisplacement( Base::StackVariables::localJacobian ),
      dLocalResidualMomentum_dPressure{ {0.0} },
      localPressureDofIndex{ 0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    // Storage for displacements

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint]{};
    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint]{};


    // Storage for helper variables used in the quadrature point kernel

    /// Strain increment
    real64 strainIncrement[6]{};

    /// Total stress
    real64 totalStress[6]{};
    /// Derivative of total stress wrt displacement
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
    /// Derivative of total stress wrt pressure
    real64 dTotalStress_dPressure[6]{};
    /// Derivative of total stress wrt temperature
    real64 dTotalStress_dTemperature[6]{};

    /// Body force
    real64 bodyForce[3]{};
    /// Derivative of body force wrt volumetric strain increment
    real64 dBodyForce_dVolStrainIncrement[3]{};
    /// Derivative of body force wrt pressure
    real64 dBodyForce_dPressure[3]{};

    /// Temperature
    real64 temperature{}; // for stress computation
    /// Delta temperature since last time step
    real64 deltaTemperatureFromLastStep{}; // for porosity update

    // Storage for residual and degrees of freedom

    /// Linear momentum balance residual
    real64 ( &localResidualMomentum )[numDispDofPerElem];
    /// Derivative of linear momentum balance residual wrt displacement
    real64 ( &dLocalResidualMomentum_dDisplacement )[numDispDofPerElem][numDispDofPerElem];
    /// Derivative of linear momentum balance residual wrt pressure
    real64 dLocalResidualMomentum_dPressure[numDispDofPerElem][1]{};

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localPressureDofIndex{};

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geos::finiteElement::ImplicitKernelBase::setup
   *
   * For the PoromechanicsBase implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    integer constexpr numDims = 3;

    m_finiteElementSpace.template setup< FE_TYPE >( k, m_meshData, stack.feStack );
    localIndex const numSupportPoints =
      m_finiteElementSpace.template numSupportPoints< FE_TYPE >( stack.feStack );

    for( localIndex a=0; a<numSupportPoints; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( integer i = 0; i < numDims; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
#endif
        stack.u_local[a][i] = m_disp[localNodeIndex][i];
        stack.uhat_local[a][i] = m_uhat[localNodeIndex][i];
        stack.localRowDofIndex[a*numDims+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*numDims+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

    stack.localPressureDofIndex = m_flowDofNumber[k];

    // Add stabilization to block diagonal parts of the local dResidualMomentum_dDisplacement (this
    // is a no-operation with FEM classes)
    real64 const stabilizationScaling = computeStabilizationScaling( k );
    m_finiteElementSpace.template addGradGradStabilizationMatrix
    < FE_TYPE, numDofPerTrialSupportPoint, false >( stack.feStack,
                                                    stack.dLocalResidualMomentum_dDisplacement,
                                                    -stabilizationScaling );
    m_finiteElementSpace.template
    addEvaluatedGradGradStabilizationVector< FE_TYPE,
                                             numDofPerTrialSupportPoint >
      ( stack.feStack,
      stack.uhat_local,
      reinterpret_cast< real64 (&)[numNodesPerElem][numDofPerTestSupportPoint] >(stack.localResidualMomentum),
      -stabilizationScaling );
  }

  /**
   * @brief Get a parameter representative of the stiffness, used as physical scaling for the
   * stabilization matrix.
   * @param[in] k Element index.
   * @return A parameter representative of the stiffness matrix dstress/dstrain
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  real64 computeStabilizationScaling( localIndex const k ) const
  {
    // TODO: generalize this to other constitutive models (currently we assume linear elasticity).
    return 2.0 * m_constitutiveUpdate.getShearModulus( k );
  }

protected:

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3]{};
  /// The L2-norm of the gravity vector
  real64 const m_gravityAcceleration;

  /// The rank global density
  arrayView2d< real64 const > m_solidDensity;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global fluid pressure at the previous converged time step
  arrayView1d< real64 const > const m_pressure_n;
  /// The rank-global fluid pressure
  arrayView1d< real64 const > const m_pressure;

};


} // namespace poromechanicsKernels

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_POROMECHANICSBASE_HPP_
