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
  m_stressOld(inputConstitutiveType.getOldStress()),
  m_matrixPresDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey ) ),
  m_fracturePresDofNumber( embeddedSurfSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey ) ),
  m_wDofNumber( jumpDofNumber ),
  m_solidDensity( inputConstitutiveType.getDensity() ),
  m_fluidDensity( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density() ),
  m_fluidDensity_n( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).density_n() ),
  m_dFluidDensity_dPressure( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >(
                                                                                                                     fluidModelKey ) ).dDensity_dPressure() ),
  m_matrixPressure( elementSubRegion.template getField< fields::flow::pressure >() ),
  m_matrixPressure_n( elementSubRegion.template getField< fields::flow::pressure_n >() ),
  m_porosity_n( inputConstitutiveType.getPorosity_n() ),
  m_tractionVec( embeddedSurfSubRegion.getField< fields::contact::traction >() ),
  m_dTraction_dJump( embeddedSurfSubRegion.getField< fields::contact::dTraction_dJump >() ),
  m_dTraction_dPressure( embeddedSurfSubRegion.getField< fields::contact::dTraction_dPressure >() ),
  m_nVec( embeddedSurfSubRegion.getNormalVector() ),
  m_tVec1( embeddedSurfSubRegion.getTangentVector1() ),
  m_tVec2( embeddedSurfSubRegion.getTangentVector2() ),
  m_surfaceCenter( embeddedSurfSubRegion.getElementCenter() ),
  m_surfaceArea( embeddedSurfSubRegion.getElementArea() ),
  m_elementVolumeCell( elementSubRegion.getElementVolume() ),
  m_elementVolumeFrac( embeddedSurfSubRegion.getElementVolume() ),
  // m_deltaVolume( elementSubRegion.template getField< fields::flow::deltaVolume >() ),
  m_deltaVolume( embeddedSurfSubRegion.template getField< fields::flow::deltaVolume >() ),
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
      std::cout << "efem k = " << k << " "<< i <<  " " << q << std::endl;
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

  // std::cout << "size of m_elementVolumeCell = " << m_elementVolumeCell.size() << std::endl;
  // std::cout << "size of m_elementVolumeFrac = " << m_elementVolumeFrac.size() << std::endl;

  stack.hInv = m_surfaceArea[embSurfIndex] / m_elementVolumeCell[k];
  // stack.hInv = m_surfaceArea[embSurfIndex] / m_elementVolumeFrac[embSurfIndex];

  // std::cout << "hInv = " << stack.hInv << std::endl;
  // std::cout << "m_surfaceArea[embSurfIndex] = " << m_surfaceArea[embSurfIndex] << std::endl;
  // std::cout << "m_elementVolume[k] = " << m_elementVolumeCell[k] << std::endl;
  // std::cout << "m_elementVolume[embSurfIndex] = " << m_elementVolumeFrac[embSurfIndex] << std::endl;
  // std::cout << "embSurfIndex = " << embSurfIndex << std::endl;

  for( localIndex a=0; a<numNodesPerElem; ++a )
  {
    localIndex const localNodeIndex = m_elemsToNodes( k, a );

    for( int i=0; i<3; ++i )
    {
      stack.dispEqnRowIndices[a*3+i] = m_dofNumber[localNodeIndex]+i-m_dofRankOffset;
      stack.dispColIndices[a*3+i]    = m_dofNumber[localNodeIndex]+i;
      stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
      stack.dispLocal[ a*3 + i ] = m_disp[ localNodeIndex ][ i ];
      stack.deltaDispLocal[ a*3 + i ] = m_deltaDisp[ localNodeIndex ][ i ];
    }
  }

  for( int i=0; i<3; ++i )
  {
    // need to grab the index.
    stack.jumpEqnRowIndices[i] = m_wDofNumber[embSurfIndex] + i - m_dofRankOffset;
    stack.jumpColIndices[i]    = m_wDofNumber[embSurfIndex] + i;
    stack.wLocal[ i ] = m_w[ embSurfIndex ][i];
    stack.tractionVec[ i ] = m_tractionVec[ embSurfIndex ][i] * m_surfaceArea[embSurfIndex];
    std::cout << "stack.tractionVec " << m_tractionVec[embSurfIndex][i] << " " << m_surfaceArea[embSurfIndex] << std::endl;
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

  // Gauss contirbution to eqMStrOld = EqMatrix*Sigma_old (3x6) *(6*1) in Voigt notation
  real64 eqMStrOld_gauss[3]{};
  real64 oldStress[6] = {m_stressOld[k][q][0], 
                          m_stressOld[k][q][1], 
                          m_stressOld[k][q][2],
                          m_stressOld[k][q][3],
                          m_stressOld[k][q][4],
                          m_stressOld[k][q][5]};
  std::cout << "old stress: ";
  for (std::size_t i = 0; i < 6; ++i)
    std::cout << oldStress[i] << " ";
  std::cout << std::endl;

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

  // E*StressOld
  LvArray::tensorOps::Ri_eq_AijBj<3, 6> (eqMStrOld_gauss, eqMatrix, oldStress);

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
  LvArray::tensorOps::scaledAdd< 3 > (stack.localEqMStrold, eqMStrOld_gauss, -detJ);

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

  // rm later!!
   std::cout << "printing localJumpResidual: " << std::endl;
   for (int i = 0; i < 3; i++)
   {
     std::cout << "localJumpResidual[" << i << "] = " << stack.localJumpResidual[i] << std::endl;
   }

  // std::cout << "printing localKww: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   for (int j = 0; j < 3; j++)
  //   {
  //     std::cout << "localKww[" << i << "," << j << "] = " << stack.localKww[i][j] << std::endl;
  //   }
  // }

  // std::cout << "printing localKwu: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   for (int j = 0; j < 3; j++)
  //   {
  //     std::cout << "localKwu[" << i << "," << j << "] = " << stack.localKwu[i][j] << std::endl;
  //   }
  // }

  // std::cout << "printing localKuw: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   for (int j = 0; j < 3; j++)
  //   {
  //     std::cout << "localKuw[" << i << "," << j << "] = " << stack.localKuw[i][j] << std::endl;
  //   }
  // }

   std::cout << "printing wLocal: " << std::endl;
   for (int i = 0; i < 3; i++)
   {
     std::cout << "wLocal[" << i << "] = " << stack.wLocal[i] << std::endl;
   }

  // std::cout << "printing localDispResidual: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   std::cout << "localDispResidual[" << i << "] = " << stack.localDispResidual[i] << std::endl;
  // }

   std::cout << "printing delta dispLocal: " << std::endl;
   for (int i = 0; i < 3; i++)
   {
     std::cout << "dispLocal[" << i << "] = " << stack.deltaDispLocal[i] << std::endl;
   }

   std::cout << "printing localEqM stress old: " << std::endl;
   for (int i = 0; i < 3; i++)
   {
     std::cout << "dispLocal[" << i << "] = " << stack.localEqMStrold[i] << std::endl;
   }

  // Compute the local residuals
  LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( stack.localJumpResidual, stack.localKww, stack.wLocal );
  LvArray::tensorOps::Ri_add_AijBj< 3, nUdof >( stack.localJumpResidual, stack.localKwu, stack.deltaDispLocal );
  LvArray::tensorOps::Ri_add_AijBj< nUdof, 3 >( stack.localDispResidual, stack.localKuw, stack.wLocal );

  // add EqM * Stress_old into the residual of enrichment nodes
  LvArray::tensorOps::add< 3 >( stack.localJumpResidual, stack.localEqMStrold);

    // rm later!!
  // std::cout << "printing localJumpResidual: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   std::cout << "localJumpResidual[" << i << "] = " << stack.localJumpResidual[i] << std::endl;
  // }

  // std::cout << "printing localDispResidual: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   std::cout << "localDispResidual[" << i << "] = " << stack.localDispResidual[i] << std::endl;
  // }

  std::cout << "printing m_matrixPressure[ k ]: " << std::endl;
  std::cout << "m_matrixPressure[" << k << "] = " << m_matrixPressure[ k ] << " " << m_matrixPressure_n[ k ] << std::endl;
  
  // std::cout << "printing localKwpm: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   std::cout << "localKwpm[" << i << "] = " << stack.localKwpm[i] << std::endl;
  // }

  // add pore pressure contribution
  LvArray::tensorOps::scaledAdd< 3 >( stack.localJumpResidual, stack.localKwpm, m_matrixPressure[ k ] );

  // rm later!!
  // std::cout << "printing localJumpResidual after matrix pressure contribution added: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   std::cout << "localJumpResidual[" << i << "] = " << stack.localJumpResidual[i] << std::endl;
  // }

  localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

  // rm later!!
  std::cout << "printing tractionVec: " << std::endl;
  for (int i = 0; i < 3; i++)
  {
    std::cout << "tractionVec[" << i << "] = " << stack.tractionVec[i] << std::endl;
  }

  // std::cout << "printing dTractiondw: " << std::endl;
  // for (int i = 0; i < 3; i++)
  // {
  //   for (int j = 0; j < 3; j++)
  //   {
  //     std::cout << "dTractiondw[" << i << "," << j << "] = " << stack.dTractiondw[i][j] << std::endl;
  //   }
  // }

  // Add traction contribution tranction
  LvArray::tensorOps::scaledAdd< 3 >( stack.localJumpResidual, stack.tractionVec, -1 );
  LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, stack.dTractiondw, -1 );

  // rm later!!
   std::cout << "printing localJumpResidual after traction added added: " << std::endl;
   for (int i = 0; i < 3; i++)
   {
     std::cout << "localJumpResidual[" << i << "] = " << stack.localJumpResidual[i] << std::endl;
   }

  // JumpFractureFlowJacobian
  real64 const localJumpFracPressureJacobian = -m_dTraction_dPressure[embSurfIndex] * m_surfaceArea[embSurfIndex];

  // Mass balance accumulation
  real64 const newVolume = m_elementVolumeFrac( embSurfIndex ) + m_deltaVolume( embSurfIndex );
  real64 const newMass =  m_fluidDensity( embSurfIndex, 0 ) * newVolume;
  real64 const oldMass =  m_fluidDensity_n( embSurfIndex, 0 ) * m_elementVolumeFrac( embSurfIndex );
  real64 const localFlowResidual = ( newMass - oldMass );
  real64 const localFlowJumpJacobian = m_fluidDensity( embSurfIndex, 0 ) * m_surfaceArea[ embSurfIndex ];
  // real64 const localFlowJumpJacobian = 0.0;
  real64 const localFlowFlowJacobian = m_dFluidDensity_dPressure( embSurfIndex, 0 ) * newVolume;


  // std::cout << "size of m_fluidDensity = " << m_fluidDensity.size() << std::endl;
  // std::cout << "size of m_fluidDensity_n = " << m_fluidDensity_n.size() << std::endl;
  // std::cout << "size of m_deltaVolume = " << m_deltaVolume.size() << std::endl;
  // std::cout << "size of m_surfaceArea = " << m_surfaceArea.size() << std::endl;
  // std::cout << "size of m_dFluidDensity_dPressure = " << m_dFluidDensity_dPressure.size() << std::endl;
  // std::cout << "size of m_dTraction_dPressure = " << m_dTraction_dPressure.size() << std::endl;

  // std::cout << "embSurfIndex = " << embSurfIndex << std::endl;
  // std::cout << "dens_n = " << m_fluidDensity_n( embSurfIndex, 0 ) << std::endl;
  // std::cout << "dens = " << m_fluidDensity( embSurfIndex, 0 ) << std::endl;
  // std::cout << "vol = " << m_elementVolumeFrac( embSurfIndex ) << std::endl;
  // std::cout << "newVolume = " << newVolume << std::endl;
  // std::cout << "localFlowResidual = " << localFlowResidual << std::endl;
  // std::cout << "m_fluidDensity * m_surfaceArea = " << m_fluidDensity( embSurfIndex, 0 ) * m_surfaceArea[ embSurfIndex ] << std::endl;
  // std::cout << "localFlowJumpJacobian = " << localFlowJumpJacobian << std::endl;
  // std::cout << "localFlowFlowJacobian = " << localFlowFlowJacobian << std::endl;

  // std::cout << "nUdof = " << nUdof << std::endl;

  for( localIndex i = 0; i < nUdof; ++i )
  {
    localIndex const uDof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );
    if( uDof < 0 || uDof >= m_matrix.numRows() )
      continue;

    RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[uDof], stack.localDispResidual[i] );

    // std::cout << "i = " << i << ", uDof = " << uDof << std::endl;

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

    // std::cout << "i = " << i << ", dof = " << dof << std::endl;

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

    // std::cout << "fracturePressureDof = " << fracturePressureDof << ", m_fracturePresDofNumber[ embSurfIndex ] = " << m_fracturePresDofNumber[ embSurfIndex ] << std::endl;

    RAJA::atomicAdd< serialAtomic >( &m_rhs[ fracturePressureDof ], localFlowResidual );
  }

  return maxForce;
}


} // namespace poromechanicsEFEMKernels

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_SINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_
