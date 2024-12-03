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
 * @file ThermalSinglePhasePoromechanicsEFEM_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_

#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsEFEM.hpp"

namespace geos
{

namespace thermoPoromechanicsEFEMKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
ThermalSinglePhasePoromechanicsEFEM( NodeManager const & nodeManager,
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
        embeddedSurfSubRegion,
        dispDofNumber,
        jumpDofNumber,
        inputFlowDofKey,
        rankOffset,
        inputMatrix,
        inputRhs,
        inputDt,
        inputGravityVector,
        fluidModelKey ),
  m_dFluidDensity_dTemperature( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >(
                                                                                                                        fluidModelKey ) ).dDensity_dTemperature() ),
  m_fluidInternalEnergy_n( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).internalEnergy_n() ),
  m_fluidInternalEnergy( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >( fluidModelKey ) ).internalEnergy() ),
  m_dFluidInternalEnergy_dPressure( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >(
                                                                                                                            fluidModelKey ) ).dInternalEnergy_dPressure() ),
  m_dFluidInternalEnergy_dTemperature( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( elementSubRegion.template getReference< string >(
                                                                                                                               fluidModelKey ) ).dInternalEnergy_dTemperature() ),
  m_temperature_n( embeddedSurfSubRegion.template getField< fields::flow::temperature_n >() ),
  m_temperature( embeddedSurfSubRegion.template getField< fields::flow::temperature >() ),
  m_matrixTemperature( elementSubRegion.template getField< fields::flow::temperature >() )

{}


//START_kernelLauncher
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
kernelLaunch( localIndex const numElems,
              KERNEL_TYPE const & kernelComponent )
{
  return Base::template kernelLaunch< POLICY >( numElems, kernelComponent );
}
//END_kernelLauncher


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  Base::setup( k, stack );
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack ) const
{

  Base::quadraturePointKernel( k, q, stack, [&] ( real64 const (&eqMatrix)[3][6],
                                                  real64 const detJ )
  {
    real64 KwTm_gauss[3]{};
    real64 thermalExpansionCoefficient{};

    m_constitutiveUpdate.getThermalExpansionCoefficient( k, thermalExpansionCoefficient );

    // assemble KwTmLocal
    LvArray::tensorOps::fill< 3 >( KwTm_gauss, 0 );
    for( int i=0; i < 3; ++i )
    {
      KwTm_gauss[0] += eqMatrix[0][i];
      KwTm_gauss[1] += eqMatrix[1][i];
      KwTm_gauss[2] += eqMatrix[2][i];
    }
    LvArray::tensorOps::scaledAdd< 3 >( stack.localKwTm, KwTm_gauss, 3*detJ*thermalExpansionCoefficient );
  } );

}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  real64 const maxForce = Base::complete( k, stack );


  // add pore pressure contribution
  real64 localJumpResidual_tempContribution[3]{};
  LvArray::tensorOps::scaledAdd< 3 >( localJumpResidual_tempContribution, stack.localKwTm, m_matrixTemperature[ k ] );

  globalIndex const matrixTemperatureColIndex = m_matrixPresDofNumber[k] + 1;
  for( localIndex i=0; i < 3; ++i )
  {
    localIndex const dof = LvArray::integerConversion< localIndex >( stack.jumpEqnRowIndices[ i ] );

    if( dof < 0 || dof >= m_matrix.numRows() )
      continue;

    RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], localJumpResidual_tempContribution[i] );

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                            &matrixTemperatureColIndex,
                                                                            &stack.localKwTm[i],
                                                                            1 );
  }

  localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];
  // Energy balance accumulation
  real64 const volume        =  m_elementVolumeFrac( embSurfIndex ) + m_deltaVolume( embSurfIndex );
  real64 const volume_n      =  m_elementVolumeFrac( embSurfIndex );
  real64 const fluidEnergy   =  m_fluidDensity( embSurfIndex, 0 ) * m_fluidInternalEnergy( embSurfIndex, 0 ) * volume;
  real64 const fluidEnergy_n =  m_fluidDensity_n( embSurfIndex, 0 ) * m_fluidInternalEnergy_n( embSurfIndex, 0 ) * volume_n;

  stack.dFluidMassIncrement_dTemperature =  m_dFluidDensity_dTemperature( embSurfIndex, 0 ) * volume;

  stack.energyIncrement               = fluidEnergy - fluidEnergy_n;
  stack.dEnergyIncrement_dJump        = m_fluidDensity( embSurfIndex, 0 ) * m_fluidInternalEnergy( embSurfIndex, 0 ) * m_surfaceArea[ embSurfIndex ];
  stack.dEnergyIncrement_dPressure    = m_dFluidDensity_dPressure( embSurfIndex, 0 ) * m_fluidInternalEnergy( embSurfIndex, 0 ) * volume;
  stack.dEnergyIncrement_dTemperature = ( m_dFluidDensity_dTemperature( embSurfIndex, 0 ) * m_fluidInternalEnergy( embSurfIndex, 0 ) +
                                          m_fluidDensity( embSurfIndex, 0 ) * m_dFluidInternalEnergy_dTemperature( embSurfIndex, 0 )  ) * volume;

  globalIndex const fracturePressureDof        = m_fracturePresDofNumber[ embSurfIndex ];
  globalIndex const fractureTemperatureDof     = m_fracturePresDofNumber[ embSurfIndex ] + 1;
  localIndex const massBalanceEquationIndex   = fracturePressureDof - m_dofRankOffset;
  localIndex const energyBalanceEquationIndex = massBalanceEquationIndex + 1;

  if( massBalanceEquationIndex >= 0 && massBalanceEquationIndex < m_matrix.numRows() )
  {

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( massBalanceEquationIndex,
                                                                            &fractureTemperatureDof,
                                                                            &stack.dFluidMassIncrement_dTemperature,
                                                                            1 );
  }

  if( energyBalanceEquationIndex >= 0 && energyBalanceEquationIndex < m_matrix.numRows() )
  {

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( energyBalanceEquationIndex,
                                                                            &stack.jumpColIndices[0],
                                                                            &stack.dEnergyIncrement_dJump,
                                                                            1 );

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( energyBalanceEquationIndex,
                                                                            &fracturePressureDof,
                                                                            &stack.dEnergyIncrement_dPressure,
                                                                            1 );

    m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( energyBalanceEquationIndex,
                                                                            &fractureTemperatureDof,
                                                                            &stack.dEnergyIncrement_dTemperature,
                                                                            1 );

    RAJA::atomicAdd< serialAtomic >( &m_rhs[ energyBalanceEquationIndex ], stack.energyIncrement );
  }

  return maxForce;
}


} // namespace thermoPoromechanicsEFEMKernels

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_
