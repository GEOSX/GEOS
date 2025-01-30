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
 * @file PhaseFieldPressurizedDamageKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDPRESSURIZEDDAMAGEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDPRESSURIZEDDAMAGEKERNELS_HPP_

#include "physicsSolvers/simplePDE/PhaseFieldDamageFEMKernels.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"

namespace geos
{

/**
 * @brief Implements kernels for solving the Damage(or phase-field) equation
 * in a phase-field pressurized fracture problem.
 * @copydoc geos::finiteElement::ImplicitKernelBase
 *
 * ### PhaseFieldPressurizedDamageKernel Description
 * Implements the KernelBase interface functions required for solving the
 * the Damage(or phase-field) equation in a phase-field pressurized fracture
 * problem using one of the "finite element kernel application" functions such as
 * geos::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class PhaseFieldPressurizedDamageKernel :
  public PhaseFieldDamageKernel< SUBREGION_TYPE,
                                 CONSTITUTIVE_TYPE,
                                 FE_TYPE >
{
public:
  /// An alias for the base class.
  using Base = PhaseFieldDamageKernel< SUBREGION_TYPE,
                                       CONSTITUTIVE_TYPE,
                                       FE_TYPE >;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;
  using Base::m_dt;
  using Base::m_X;
  using Base::m_nodalDamage;
  using Base::m_quadDamage;
  using Base::m_quadExtDrivingForce;
  using Base::m_localDissipationOption;
  using LocalDissipation = PhaseFieldDamageKernelLocalDissipation;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  PhaseFieldPressurizedDamageKernel( NodeManager const & nodeManager,
                                     EdgeManager const & edgeManager,
                                     FaceManager const & faceManager,
                                     localIndex const targetRegionIndex,
                                     SUBREGION_TYPE const & elementSubRegion,
                                     FE_TYPE const & finiteElementSpace,
                                     CONSTITUTIVE_TYPE & inputConstitutiveType,
                                     arrayView1d< globalIndex const > const inputDofNumber,
                                     globalIndex const rankOffset,
                                     CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                     arrayView1d< real64 > const inputRhs,
                                     real64 const inputDt,
                                     string const fieldName,
                                     LocalDissipation localDissipationOption ):
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
          inputDt,
          fieldName,
          localDissipationOption ),
    m_disp( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
    m_fluidPressure( elementSubRegion.template getField< fields::flow::pressure >()  ),
    m_fluidPressureGradient( elementSubRegion.template getReference< array2d< real64 > >( "pressureGradient" ) )
  {}

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the primary field.
   */
  struct StackVariables : Base::StackVariables
  {
public:

    /**
     * @brief Constructor
     */
    GEOS_HOST_DEVICE
    StackVariables():
      Base::StackVariables()
    {}

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][3];
  };


  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geosx::finiteElement::ImplicitKernelBase::setup
   */
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    Base::setup( k, stack );

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );
      LvArray::tensorOps::copy< 3 >( stack.u_local[ a ], m_disp[ localNodeIndex ] );
    }
  }

  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    Base::quadraturePointKernel( k, q, stack );

    real64 const ell = m_constitutiveUpdate.getRegularizationLength();
    real64 const Gc = m_constitutiveUpdate.getCriticalFractureEnergy( k );
    real64 const volStrain = m_constitutiveUpdate.getVolStrain( k, q );
    real64 const biotCoeff = m_constitutiveUpdate.getBiotCoefficient( k );

    //Interpolate d and grad_d
    real64 N[ numNodesPerElem ];
    real64 dNdX[ numNodesPerElem ][ 3 ];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );
    FE_TYPE::calcN( q, N );

    real64 qp_damage = 0.0;
    FE_TYPE::value( N, stack.nodalDamageLocal, qp_damage );

    real64 qp_disp[3] = {0, 0, 0};
    FE_TYPE::value( N, stack.u_local, qp_disp );

    real64 elemPresGradient[3] = {0, 0, 0};
    for( integer i=0; i<3; ++i )
    {
      elemPresGradient[i] = m_fluidPressureGradient( k, i );
    }

    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      /// Add pressure effects
      stack.localResidual[ a ] -= detJ * 0.5 * ell/Gc * ( ( 1.0 - biotCoeff ) * volStrain * m_fluidPressure( k ) * m_constitutiveUpdate.pressureDamageFunctionDerivative( qp_damage ) * N[a]
                                                          + LvArray::tensorOps::AiBi< 3 >( qp_disp, elemPresGradient ) * m_constitutiveUpdate.pressureDamageFunctionDerivative( qp_damage ) * N[a] );

      for( localIndex b = 0; b < numNodesPerElem; ++b )
      {
        stack.localJacobian[ a ][ b ] -= detJ * 0.5 * ell/Gc *
                                         ( ( 1.0 - biotCoeff ) * volStrain * m_fluidPressure( k ) * m_constitutiveUpdate.pressureDamageFunctionSecondDerivative( qp_damage ) * N[a] * N[b]
                                           + LvArray::tensorOps::AiBi< 3 >( qp_disp, elemPresGradient ) * m_constitutiveUpdate.pressureDamageFunctionSecondDerivative( qp_damage ) * N[a] * N[b] );
      }
    }
  }

protected:

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  arrayView1d< real64 const > const m_fluidPressure;

  arrayView2d< real64 const > const m_fluidPressureGradient;

};

using PhaseFieldPressurizedDamageKernelFactory = finiteElement::KernelFactory< PhaseFieldPressurizedDamageKernel,
                                                                               arrayView1d< globalIndex const > const,
                                                                               globalIndex,
                                                                               CRSMatrixView< real64, globalIndex const > const,
                                                                               arrayView1d< real64 > const,
                                                                               real64 const,
                                                                               string const,
                                                                               PhaseFieldDamageKernelLocalDissipation >;

} // namespace geos

#endif // GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDPRESSURIZEDDAMAGEKERNELS_HPP_
