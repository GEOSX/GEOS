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
 * @file SinglePhasePoromechanicsEFEM_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_

#include "physicsSolvers/contact/ContactFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsEFEM.hpp"
#include "physicsSolvers/contact/SolidMechanicsEFEMKernelsHelper.hpp"

namespace geos
{

namespace poromechanicsEFEMKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
SinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
SinglePhasePoromechanicsEFEM( NodeManager const & nodeManager,
                              EdgeManager const & edgeManager,
                              FaceManager const & faceManager,
                              localIndex const targetRegionIndex,
                              SUBREGION_TYPE const & elementSubRegion,
                              FE_TYPE const & finiteElementSpace,
                              CONSTITUTIVE_TYPE & inputConstitutiveType,
                              EmbeddedSurfaceSubRegion const & embeddedSurfSubRegion,
                              arrayView1d< globalIndex const > const dispDofNumber,
                              arrayView1d< globalIndex const > const jumpDofNumber,
                              string const inputFlowDofKey,
                              globalIndex const rankOffset,
                              CRSMatrixView< real64, globalIndex const > const inputMatrix,
                              arrayView1d< real64 > const inputRhs,
                              real64 const inputDt,
                              real64 const (&inputGravityVector)[3],
                              string const fluidModelKey ):
  Base( nodeManager,
        edgeManager,
        faceManager,
        targetRegionIndex,
        elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType,
        dispDofNumber,
        rankOffset,
        inputMatrix,
        inputRhs,
        inputDt ),
  m_X( nodeManager.referencePosition()),
  m_disp( nodeManager.getField< fields::solidMechanics::totalDisplacement >() ),
  m_deltaDisp( nodeManager.getField< fields::solidMechanics::incrementalDisplacement >() ),
  m_w( embeddedSurfSubRegion.getField< fields::contact::dispJump >() ),
  m_effStress(inputConstitutiveType.getEffectiveStress()),
  m_matrixPresDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey ) ),
  m_fracturePresDofNumber( embeddedSurfSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey ) ),
  m_wDofNumber( jumpDofNumber ),
  m_solidDensity( inputConstitutiveType.getDensity() ),
  m_fluidDensity( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density() ),
  m_fluidDensity_n( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density_n() ),
  m_dFluidDensity_dPressure( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >(
                                                                                                                     fluidModelKey ) ).dDensity_dPressure() ),
  m_matrixPressure( elementSubRegion.template getField< fields::flow::pressure >() ),
  m_fracturePressure( embeddedSurfSubRegion.template getField< fields::flow::pressure >() ),
  m_porosity_n( inputConstitutiveType.getPorosity_n() ),
  m_tractionVec( embeddedSurfSubRegion.getField< fields::contact::traction >() ),
  m_dTraction_dJump( embeddedSurfSubRegion.getField< fields::contact::dTraction_dJump >() ),
  m_dTraction_dPressure( embeddedSurfSubRegion.getField< fields::contact::dTraction_dPressure >() ),
  m_nVec( embeddedSurfSubRegion.getNormalVector() ),
  m_tVec1( embeddedSurfSubRegion.getTangentVector1() ),
  m_tVec2( embeddedSurfSubRegion.getTangentVector2() ),
  m_surfaceCenter( embeddedSurfSubRegion.getElementCenter() ),
  m_surfaceArea( embeddedSurfSubRegion.getElementArea() ),
  m_elementVolume( elementSubRegion.getElementVolume() ),
  m_deltaVolume( elementSubRegion.template getField< fields::flow::deltaVolume >() ),
  m_fracturedElems( elementSubRegion.fracturedElementsList() ),
  m_cellsToEmbeddedSurfaces( elementSubRegion.embeddedSurfacesList().toViewConst() ),
  m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
  m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) )
{}


//START_kernelLauncher
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
SinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
kernelLaunch( localIndex const numElems,
              KERNEL_TYPE const & kernelComponent )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( numElems );

  // Define a RAJA reduction variable to get the maximum residual contribution.
  RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

  forAll< POLICY >( kernelComponent.m_fracturedElems.size(),
                    [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localIndex k = kernelComponent.m_fracturedElems[i];
    typename KERNEL_TYPE::StackVariables stack;

    kernelComponent.setup( k, stack );
    for( integer q=0; q<numQuadraturePointsPerElem; ++q )
    {
      kernelComponent.quadraturePointKernel( k, q, stack );
    }
    maxResidual.max( kernelComponent.complete( k, stack ) );
  } );

  return maxResidual.get();
}
//END_kernelLauncher


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

  stack.hInv = m_surfaceArea[embSurfIndex] / m_elementVolume[k];
  for( localIndex a=0; a<numNodesPerElem; ++a )
  {
    localIndex const localNodeIndex = m_elemsToNodes( k, a );

    for( int i=0; i<3; ++i )
    {
      stack.dispEqnRowIndices[a*3+i] = m_dofNumber[localNodeIndex]+i-m_dofRankOffset;
      stack.dispColIndices[a*3+i]    = m_dofNumber[localNodeIndex]+i;
      stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
      stack.dispLocal[ a*3 + i ] = m_disp[ localNodeIndex ][ i ];
      stack.deltaDispLocal[ a ][ i ] = m_deltaDisp[ localNodeIndex ][ i ];
    }
  }

  for( int i=0; i<3; ++i )
  {
    // need to grab the index.
    stack.jumpEqnRowIndices[i] = m_wDofNumber[embSurfIndex] + i - m_dofRankOffset;
    stack.jumpColIndices[i]    = m_wDofNumber[embSurfIndex] + i;
    stack.wLocal[ i ] = m_w[ embSurfIndex ][i];
    stack.tractionVec[ i ] = m_tractionVec[ embSurfIndex ][i] * m_surfaceArea[embSurfIndex];
    for( int ii=0; ii < 3; ++ii )
    {
      stack.dTractiondw[ i ][ ii ] = m_dTraction_dJump[embSurfIndex][i][ii] * m_surfaceArea[embSurfIndex];
    }
  }
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename FUNC >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void SinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack,
                       FUNC && kernelOp ) const
{

  // The quarature kernal deals with the fracture force balance (eq. 29 in https://onlinelibrary.wiley.com/doi/epdf/10.1002/nag.3168)
  // The total stress in matrix: sigma_tau = simga_eff_old + sigma_eff_incr - biot * p_m = sigma_eff_new - biot * p_m
  // We can either use formulation: sigma_eff_old + stiffness_matrix * incremental_strain - biot * p_m or directly effective_stress_current - biot * p_m
  // The latter one is adopted here.

  localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

  // Get displacement: (i) basis functions (N), (ii) basis function
  // derivatives (dNdX), and (iii) determinant of the Jacobian transformation
  // matrix times the quadrature weight (detJxW)
  real64 dNdX[numNodesPerElem][3];
  real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

  // EFEM part starts here
  constexpr int nUdof = numNodesPerElem*3;

  // Gauss contribution to Kww, Kwu and Kuw blocks
  real64 Kww_gauss[3][3]{}, Kwu_gauss[3][nUdof]{}, Kuw_gauss[nUdof][3]{}, Kwpm_gauss[3]{};

 // Gauss contirbution to eqMStress which is EqMatrix*effStress, all stresses are in Voigt notation
  real64 eqMStress_gauss[3]{};

  //  Compatibility, equilibrium and strain operators. The compatibility operator is constructed as
  //  a 3 x 6 because it is more convenient for construction purposes (reduces number of local var).
  real64 compMatrix[3][6]{}, strainMatrix[6][nUdof]{}, eqMatrix[3][6]{};
  real64 matBD[nUdof][6]{}, matED[3][6]{};
  real64 biotCoefficient{};
  int Heaviside[ numNodesPerElem ]{};

  m_constitutiveUpdate.getBiotCoefficient( k, biotCoefficient );


  // TODO: asking for the stiffness here will only work for elastic models.  most other models
  //       need to know the strain increment to compute the current stiffness value.
  m_constitutiveUpdate.getElasticStiffness( k, q, stack.constitutiveStiffness );

  solidMechanicsEFEMKernelsHelper::computeHeavisideFunction< numNodesPerElem >( Heaviside,
                                                                                stack.xLocal,
                                                                                m_nVec[embSurfIndex],
                                                                                m_surfaceCenter[embSurfIndex] );


  solidMechanicsEFEMKernelsHelper::assembleEquilibriumOperator( eqMatrix,
                                                                m_nVec[embSurfIndex],
                                                                m_tVec1[embSurfIndex],
                                                                m_tVec2[embSurfIndex],
                                                                stack.hInv );

  solidMechanicsEFEMKernelsHelper::assembleCompatibilityOperator< numNodesPerElem >( compMatrix,
                                                                                     m_nVec[embSurfIndex],
                                                                                     m_tVec1[embSurfIndex],
                                                                                     m_tVec2[embSurfIndex],
                                                                                     Heaviside,
                                                                                     dNdX );

  solidMechanicsEFEMKernelsHelper::assembleStrainOperator< 6, nUdof, numNodesPerElem >( strainMatrix, dNdX );

  // transp(B)D
  LvArray::tensorOps::Rij_eq_AkiBkj< nUdof, 6, 6 >( matBD, strainMatrix, stack.constitutiveStiffness );
  // ED
  LvArray::tensorOps::Rij_eq_AikBkj< 3, 6, 6 >( matED, eqMatrix, stack.constitutiveStiffness );
  // EDC
  LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 6 >( Kww_gauss, matED, compMatrix );
  // EDB
  LvArray::tensorOps::Rij_eq_AikBkj< 3, nUdof, 6 >( Kwu_gauss, matED, strainMatrix );
  // transp(B)DB
  LvArray::tensorOps::Rij_eq_AikBjk< nUdof, 3, 6 >( Kuw_gauss, matBD, compMatrix );
  // EqMatrix * effStress
  LvArray::tensorOps::Ri_eq_AijBj<3, 6> (eqMStress_gauss, eqMatrix, m_effStress[k][q]);

  LvArray::tensorOps::fill< 3 >( Kwpm_gauss, 0 );
  for( int i=0; i < 3; ++i )
  {
    Kwpm_gauss[0] += eqMatrix[0][i];
    Kwpm_gauss[1] += eqMatrix[1][i];
    Kwpm_gauss[2] += eqMatrix[2][i];
  }

  // multiply by determinant and add to element matrix
  LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, Kww_gauss, -detJ );
  LvArray::tensorOps::scaledAdd< 3, nUdof >( stack.localKwu, Kwu_gauss, -detJ );
  LvArray::tensorOps::scaledAdd< nUdof, 3 >( stack.localKuw, Kuw_gauss, -detJ );
  LvArray::tensorOps::scaledAdd< 3 > (stack.localEqMStress, eqMStress_gauss, -detJ);

  /// TODO: should this be negative???
  // I had No neg coz the total stress = effective stress - porePressure
  // and all signs are flipped here.
  LvArray::tensorOps::scaledAdd< 3 >( stack.localKwpm, Kwpm_gauss, detJ*biotCoefficient );

  kernelOp( eqMatrix, detJ );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 SinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  real64 maxForce = 0;
  constexpr int nUdof = numNodesPerElem*3;

  globalIndex matrixPressureColIndex = m_matrixPresDofNumber[k];

  // Compute the local residuals
  LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( stack.localJumpResidual, stack.localKww, stack.wLocal );
  LvArray::tensorOps::Ri_add_AijBj< nUdof, 3 >( stack.localDispResidual, stack.localKuw, stack.wLocal );
  // add EqM * effStress into the residual of enrichment nodes
  LvArray::tensorOps::add< 3 >( stack.localJumpResidual, stack.localEqMStress);

  // add pore pressure contribution
  LvArray::tensorOps::scaledAdd< 3 >( stack.localJumpResidual, stack.localKwpm, m_matrixPressure[ k ] );

  localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

  // Add total traction contribution from penalty force and fracture pressure  
  // total traction is T_total = -k * dispJump + pf (where dispJump < 0)
  // -1 is because k*dispJump was saved in tractionVec
  LvArray::tensorOps::scaledAdd< 3 >( stack.localJumpResidual, stack.tractionVec, -1 );
  LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, stack.dTractiondw, -1 );
  
  // fracture pressure only affects normal direction
  stack.localJumpResidual[0] += m_fracturePressure[embSurfIndex] * m_surfaceArea[embSurfIndex];
  // fracture force balance residual w.r.t. fracture pressure 
  real64 const localJumpFracPressureJacobian = m_surfaceArea[embSurfIndex];

  // Mass balance accumulation
  real64 const newVolume = m_elementVolume( embSurfIndex ) + m_deltaVolume( embSurfIndex );
  real64 const newMass =  m_fluidDensity( embSurfIndex, 0 ) * newVolume;
  real64 const oldMass =  m_fluidDensity_n( embSurfIndex, 0 ) * m_elementVolume( embSurfIndex );
  real64 const localFlowResidual = ( newMass - oldMass );
  real64 const localFlowJumpJacobian = m_fluidDensity( embSurfIndex, 0 ) * m_surfaceArea[ embSurfIndex ];
  real64 const localFlowFlowJacobian = m_dFluidDensity_dPressure( embSurfIndex, 0 ) * newVolume;

  for( localIndex i = 0; i < nUdof; ++i )
  {
    localIndex const uDof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );
    if( uDof < 0 || uDof >= m_matrix.numRows() )
      continue;

    RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[uDof], stack.localDispResidual[i] );

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( uDof,
                                                                            stack.jumpColIndices,
                                                                            stack.localKuw[i],
                                                                            3 );

  }

  for( localIndex i=0; i < 3; ++i )
  {
    localIndex const dof = LvArray::integerConversion< localIndex >( stack.jumpEqnRowIndices[ i ] );

    if( dof < 0 || dof >= m_matrix.numRows() )
      continue;

    RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localJumpResidual[i] );

    // fill in matrix
    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                            stack.jumpColIndices,
                                                                            stack.localKww[i],
                                                                            3 );
    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                            stack.dispColIndices,
                                                                            stack.localKwu[i],
                                                                            numNodesPerElem*3 );

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                            &matrixPressureColIndex,
                                                                            &stack.localKwpm[i],
                                                                            1 );
  }

//    // it only affects the normal jump

  if( stack.jumpEqnRowIndices[0] >= 0 && stack.jumpEqnRowIndices[0] < m_matrix.numRows() )
  {

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( stack.jumpEqnRowIndices[0],
                                                                            &m_fracturePresDofNumber[ embSurfIndex ],
                                                                            &localJumpFracPressureJacobian,
                                                                            1 );
  }

  localIndex const fracturePressureDof = m_fracturePresDofNumber[ embSurfIndex ] - m_dofRankOffset;
  if( fracturePressureDof >= 0 && fracturePressureDof < m_matrix.numRows() )
  {

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( fracturePressureDof,
                                                                            &stack.jumpColIndices[0],
                                                                            &localFlowJumpJacobian,
                                                                            1 );

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( fracturePressureDof,
                                                                            &m_fracturePresDofNumber[ embSurfIndex ],
                                                                            &localFlowFlowJacobian,
                                                                            1 );

    RAJA::atomicAdd< serialAtomic >( &m_rhs[ fracturePressureDof ], localFlowResidual );
  }

  return maxForce;
}


} // namespace poromechanicsEFEMKernels

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_
